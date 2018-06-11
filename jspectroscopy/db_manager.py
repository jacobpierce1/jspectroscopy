# this script is made for managing the DB, including creating the initial schema,
# adding columns later for data analysis, wrappers for reading / writing to the 
# schema, and 
# 
# source: https://pymotw.com/3/sqlite3/
import copy 
import os
import sqlite3
import json
import numpy as np

from libjacob.jutils import isint, time_estimator
import sys
import jspectroscopy as spec 

import datetime

import dill # allow lambda function pickle

import libjacob.jmeas as meas 




# config 
DEBUG_DB = 0



# _current_abs_path = os.path.dirname( __file__ ) + '/'



schema_filename = 'initial_schema.sql'
_tablename = 'fits'


schema_cols = [ 'detnum', 'x', 'y', 'fit_id', 'successful_fit',
                'npeaks', 'last_attempt', 
                'params_guess', 'fit_bounds',
                'peak_guesses', 'model']





# class providing storage, access, and filtering capabilities
# to make the user's life less annoying.










# class providing complete set of operations for writing
# and specialized querying for the DBs.
# det_number is optional number of the detector, in case
# multiple detectors were used for the same dataset.
# sources is a list of the names of the sources used.

class spectrum_db( object ):

    
    def __init__( self, path, dets_used = None, dimensions = None,
                  peak_types = None, constrain_det_params = None,
                  name = None ) :

        create_new_db = ( peak_types is not None ) and ( dimensions is not None )
        
        self.path = path
        self.conn = None
        self.name = name

        # self.delete() 
        
        pathdir = os.path.dirname( path ) 
        if not os.path.exists( pathdir ) :
            os.mkdir( pathdir ) 
        
        if not os.path.exists( path ) :
            if create_new_db :
                self.is_empty = 1 
                self.create( dets_used, dimensions, peak_types, constrain_det_params )
            else :
                print( 'ERROR: attempted to load db, but the path does not exist' )
                sys.exit(0) 

        else : 
            self.is_empty = 0
            self.connect()
            self.read_metadata()

        # TODO remove this to fix bug
        
        return

    
        
    def create( self, dets_used, dimensions, peak_types, constrain_det_params ) :

        if dets_used is None :
            self.dets_used = [-1]
        else : 
            self.dets_used = dets_used 
        
        if constrain_det_params is None :
            constrain_det_params = { 'a' : 0, 'b' : 0, 'g' : 0 } 

        self.constrain_det_params = copy.deepcopy( constrain_det_params )  

        if self.name is None :
            name = ''

        self.name = name
        
        self.peak_types = peak_types
        self.dimensions = dimensions

        if not os.path.exists( self.path ) :
            print( 'INFO: creating DB and schema for ' + self.path + '...' )
            
        else:
            print( 'ERROR: database already exists, returning' )
            return 0

        self.is_empty = 1
        self.conn = None
        
        self.xdim = dimensions[0]
        self.ydim = dimensions[1] 
        self.dimensions = dimensions
                                     
        self.num_peaks_per_group = [ len(x) for x in peak_types ] 

        self.num_groups = len( self.num_peaks_per_group )
                
        # with sqlite3.connect( self.path ) as conn:

        self.connect()
        
        with open( os.path.dirname( __file__ ) + '/'
                   + schema_filename, 'rt' ) as f:

            schema = f.read()
            self.conn.executescript( schema )    

        print( 'INFO: success, now populating...' )
            
        self.populate()        
        self.write_metadata()

        self.is_empty = 0
        

    

            
    # fill the db with each fit id that we need, giving the needs update flag for
    # everything. this is meant to only be called when the table is empty 
    def populate( self ):

        for d in self.dets_used :  
            for x in range( self.xdim ):
                for y in range( self.ydim ):
                    for i in range( self.num_groups ):
                        self.insert_fit_data( d, x, y, i, None ) 



                        

    def write_metadata( self ) :
        
        time = str( datetime.datetime.now() )

        if self.is_empty :
            self.conn.cursor().execute( 'INSERT INTO metadata VALUES ( ?, ?, ?, ?, ?, ?, ? ) ',
                                        ( json.dumps( self.dets_used ), 
                                          self.xdim, self.ydim,
                                          _to_bin( self.peak_types ),
                                          _to_bin( self.constrain_det_params ),
                                          time, self.name ) )
            self.conn.commit() 

            

            
    def read_metadata( self ) :

        cursor = self.conn.cursor()
        cursor.execute( 'SELECT * FROM metadata' )
        metadata  = dict( cursor.fetchone() )

        self.dets_used = json.loads( metadata[ 'dets_used' ] ) 
        self.xdim = metadata[ 'xdim' ] 
        self.ydim = metadata[ 'ydim' ]
        self.peak_types = _from_bin( metadata[ 'peak_types' ] ) 
        self.constrain_det_params = _from_bin( metadata[ 'constrain_det_params' ] )
        self.timestamp = metadata[ 'timestamp' ]
        self.name = metadata[ 'name' ]

        self.num_groups = len( self.peak_types ) 
        self.num_peaks_per_group = [ len(x) for x in self.peak_types ]


        

        
        
    # check if this db has been created. 
    def exists( self ):
        return os.path.exists( self.path ) 

    

    
    def connect( self ):

        if self.conn is not None:
            print( 'ERROR: db is already open.' )
            return 0

        # throw an error for an attempt to connect to a nonexistent
        # database, unless we insist that it is supposed to be empty.
        
        if not self.is_empty :
            if not self.exists():
                print( 'ERROR: db has not been created.' )
                return 0 

        self.conn = sqlite3.connect( self.path ) 

        # this is set to allow you to read into an sqlite3.Row,
        # which is a dict-like object. 
        self.conn.row_factory = sqlite3.Row

        return 1


    

    
    def disconnect( self ):

        if self.conn is None:
            return 0 
        
        self.conn.close()
        self.conn = None
        
        return 1

    
    
    
    # call this before doing any reading or writing.
    def assert_open( self ):
        if self.conn is None:
            raise ValueError( 'db is not open.' )

    


    

    def insert_fit_data( self, detnum, x, y, group_num, spectrum_fit ):        

        if self.conn is None:
            raise ValueError( 'Cannot insert data, sqlite3 connection is not open. Call db.connect().' )

        if self.is_empty :
            query = 'INSERT INTO fits VALUES ( ?, ?, ?, ?, ?, ?, ?, ?, ?, ? )'
            self.conn.cursor().execute( query, ( detnum, x, y, group_num, -1,
                                                 None, None, None, None, -1 ) )
            
        else:
            if spectrum_fit is None :
                query = ( 'UPDATE fits SET success=? '
                          + 'WHERE detnum=? and x=? and y=? and group_num=?' )
                
                self.conn.cursor().execute( query, ( 0, detnum, x, y, group_num ) ) 
            
            else :
                query = ( 'UPDATE fits SET success=?, '
                          + 'params_result=?, params_result_errors=?, '
                          + 'cov=?, fit_bounds=?, redchisqr=? '
                          + 'WHERE detnum=? and x=? and y=? and group_num=?' )

                self.conn.cursor().execute( query, ( spectrum_fit.success,
                                                     _to_bin( spectrum_fit.params_result ),
                                                     _to_bin( spectrum_fit.params_result_errors),
                                                     _to_bin( spectrum_fit.cov ),
                                                     _to_bin( spectrum_fit.fit_bounds ),
                                                     spectrum_fit.redchisqr,
                                                     detnum, x, y, group_num ) )

            self.conn.commit()
        
        return None


        


    def read_fit_data( self, detnum, x, y, fit_id ):

        # print( (x,y,fit_id) )

        if self.conn is None:
            raise ValueError( 'Cannot read db, sqlite3 connection is not open. Call db.connect().' )

        query = 'SELECT * FROM fits WHERE detnum=? and x=? and y=? and fit_id=?'

        cursor = self.conn.cursor()
        cursor.execute( query, ( detnum, x, y, fit_id ) )
        result = dict( cursor.fetchone() )
        
        spec_result = _from_bin( result[ 'spectrum_fit' ] )

        print( result )

        sys.exit(0)

        # print( spec_result ) 

        return spec_result 



    



    def delete( self ):
        
        if os.path.exists( self.path ):
            
            # ans = raw_input( 'PROMPT: delete ' + filename + ', are you sure (y/n) ?  ' )
            
            # if( ans == 'y' ):

            os.remove( self.path ) 
            print( 'INFO: deleted db.' )         
            return 1

            print( 'INFO: did not delete db.' )

            return 0

        return 1
    


    
    def write_mu_values( self, path ) :

        # init the data structure storing the mu values
        
        mu_values = [ [ 0 ] * self.num_peaks_per_group[i]
                      for i in range( self.num_groups ) ] 

        for i in range( self.num_groups ) :
            for j in range( self.num_peaks_per_group[i] ) : 
                mu_values[i][j] = meas.meas.empty(( self.xdim, self.ydim ))
                            
            
        cursor = self.conn.cursor()
        cursor.execute( 'SELECT * FROM ' + _tablename )

        spec_fitters = [ spec.spectrum_fitter( self.peak_types[i],
                                               constrain_det_params = { 'a' : 1 } )
                         for i in range( self.num_groups ) ]
        
        for row in cursor:

            x = row['x']
            y = row['y'] 
            group_num = row[ 'group_num' ]
            success = row[ 'success' ]
            
            if success :
                params_result = _from_bin( row[ 'params_result' ] )
                params_result_errors = _from_bin( row[ 'params_result_errors' ] )
                
                spectrum = spec_fitters[ group_num ] 

                spectrum.set_params_from_params_array( params_result, errors = 0)
                spectrum.set_params_from_params_array( params_result_errors, errors = 1)

                # print( spectrum.peak_params[0][1] )
                # print( spectrum.peak_params[1][1] )
                # print( '' ) 
                    
                
                for j in range( self.num_peaks_per_group[ group_num ] ) :
                    mu_values[ group_num ][j][x,y] = meas.meas( spectrum.peak_params[j][1],
                                                                spectrum.peak_params_errors[j][1] ) 
                    
            else :  
                for j in range( self.num_peaks_per_group[ group_num ] ) :
                    mu_values[ group_num ][j][x,y] = meas.nan 
                        
        with open( path, 'wb' ) as f :
            dill.dump( mu_values, f )




            
    def write_peak_values( self, path ) :

        # init the data structure storing the mu values
        
        peak_values = [ [ 0 ] * self.num_peaks_per_group[i]
                      for i in range( self.num_groups ) ] 

        for i in range( self.num_groups ) :
            for j in range( self.num_peaks_per_group[i] ) : 
                peak_values[i][j] = meas.meas.empty(( self.xdim, self.ydim ))
                            
            
        cursor = self.conn.cursor()
        cursor.execute( 'SELECT * FROM ' + _tablename )

        spec_fitters = [ spec.spectrum_fitter( self.peak_types[i],
                                               constrain_det_params = { 'a' : 1 } )
                         for i in range( self.num_groups ) ]
        
        for row in cursor:

            x = row['x']
            y = row['y'] 
            group_num = row[ 'group_num' ]
            success = row[ 'success' ]

            # print( 'success: ' + str( success ) ) 
                        
            if success :
                
                cov = _from_bin( row[ 'cov' ] ) 

                params_result = _from_bin( row[ 'params_result' ] )
                params_result_errors = _from_bin( row[ 'params_result_errors' ] )
                
                spectrum = spec_fitters[ group_num ] 

                spectrum.set_params_from_params_array( params_result, errors = 0)
                spectrum.set_params_from_params_array( params_result_errors, errors = 1)
                spectrum.cov = cov 
                
                for j in range( self.num_peaks_per_group[ group_num ] ) :

                    # print( spectrum.cov)
                    # print( 'params_result: ' + str( params_result ) ) 
                    # print( spectrum._construct_params_array() )
                    # print( spectrum.peak_params[0] )
                                        
                    if spectrum.cov is not None : 
                        alpha_peak = spec.find_alpha_peak( spectrum, j, 500, 0 ) 

                        if alpha_peak is None :
                            peak_values[ group_num ][j][x,y] = meas.nan

                        else : 
                            peak_values[ group_num ][j][x,y] = alpha_peak 

                    else :
                        peak_values[ group_num ][j][x,y] = meas.nan
                        
            else :  
                for j in range( self.num_peaks_per_group[ group_num ] ) :
                    peak_values[ group_num ][j][x,y] = meas.nan 
                        
        with open( path, 'wb' ) as f :
            dill.dump( peak_values, f )


            
        

            
    def read_values( self, path ) :

        if not os.path.exists( path ) :
            print( 'ERROR: path not found' ) 
            sys.exit(0) 

        with open( path, 'rb' ) as f :
            values = dill.load( f )

        return values 

    


    
    # todo                 
    def gen_csv( self, path ) :
        pass 



    
    # get the array of mu values for each pixel and for each peak.
    
    def get_all_mu_grids( self, read_from_file = 0 ) :

        # construct appropriately sized container for all the mu grids
        mu_grids = [ [ 0 ] * self.num_peaks_per_feature[i]
                     for i in range( self.num_features ) ]

        for i in range( self.num_features ):
            mu_grids[i] = [ meas.meas.empty( self.dimensions ) ] * self.num_peaks_per_feature[i] 


        # meas.meas.empty( self.num_peaks_per_fit + self.dimensions )

        if read_from_file :

            mu_paths = [ [ [  self.mu_vals_dir + self.name + '_%d_%d_%s.bin' % (i,j,s)
                              for s in [ 'x', 'dx' ] ]
                           for j in range( self.num_peaks_per_feature[i] ) ]
                         for i in range( self.num_features) ]
                        
            if all( [ os.path.exists( path ) for path in np.array( mu_paths ).flatten() ] ) :

                for i in range( self.num_features ) :
                    for j in range( self.num_peaks_per_feature[i] ) :

                        mu = np.fromfile( mu_paths[i][j][0] ).reshape( self.dimensions )
                        mu_delta = np.fromfile( mu_paths[i][j][1] ).reshape( self.dimensions )
                        mu_grids[i][j] = meas.meas( mu, mu_delta )
                        
                return mu_grids

            else:
                print( 'INFO: requested to read mu and mu_delta from files, but they aren\'t there. constructing them now...' )


                

        # construct the array from the db if the files don't exist
        # yet, or we requested a fresh fetch.

        disconnect_conn_when_done = 0
        if self.conn is None:
            self.connect()
            disconnect_conn_when_done = 1

            
        # populate the grid 
        for i in range(3):
            for j in range(2):
                mu_grids[i][j] =  meas.meas.empty( self.dimensions )
        
        cursor = self.conn.cursor()
        cursor.execute( 'SELECT * FROM ' + _tablename )

        for row in cursor:

            x = row['x']
            y = row['y'] 
            feature = row[ 'fit_id' ]
            
            # if fit didn't converge then all the mu values for that
            # feature are assigned np.nan            
            if not row[ 'successful_fit' ]: 

                for i in range( self.num_peaks_per_feature[ feature ]  ):
                    mu_grids[ feature ][ i ][ x,y ] = meas.nan
                continue

            params = spec.get_alpha_params_dict( _from_bin( row['model'] ) ) 

            if not spec.alpha_model_valid_params_predicate( params ) :
            
                for i in range( self.num_peaks_per_feature[ feature ]  ):
                    mu_grids[ feature ][ i ][ x,y ] = meas.nan
                continue

            mu = params[ 'mu' ]

            if len( mu ) > 1 :
                if ( np.abs( mu[1].x - mu[0].x - 20 ) > 10
                     or any( mu.dx < 0.1 ) ):
                    
                    for i in range( self.num_peaks_per_feature[ feature ]  ):
                        mu_grids[ feature ][ i ][ x,y ] = meas.nan
                    continue
            
            # check if the fit used the alternate shape 
            if len( mu ) != self.num_peaks_per_feature[ feature ] :

                # todo: this doesn't work in the general case.
                # this should depend on the specified alternate shape.
                mu = meas.append( meas.nan, mu ) 


            # populate the corresponding entry of the grid.
            for i in range( self.num_peaks_per_feature[ feature ] ):
                mu_grids[ feature ][ i ][ x,y ] = mu[i]

                
        if disconnect_conn_when_done:
            self.disconnect()


        # write to a file if requested. note that we have already constructed
        # all the file names.
        if read_from_file :
            
            for i in range( self.num_features ) :
                for j in range( self.num_peaks_per_feature[i] ) :
                    
                    mu_grids[i][j].x.tofile(  mu_paths[i][j][0] )
                    mu_grids[i][j].dx.tofile( mu_paths[i][j][1] )
                    
        
        return mu_grids



    




    
    

    def get_all_peak_grids( self, read_from_file = 1 ) :
        
        # construct appropriately sized container for all the peak grids
        peak_grids = [ [ 0 ] * self.num_peaks_per_feature[i]
                     for i in range( self.num_features ) ]

        for i in range( self.num_features ):
            peak_grids[i] = [ meas.meas.empty( self.dimensions ) ] * self.num_peaks_per_feature[i] 


        # meas.meas.empty( self.num_peaks_per_fit + self.dimensions )

        if read_from_file :

            peak_paths = [ [ [  self.peak_vals_dir + self.name + '_%d_%d_%s.bin' % (i,j,s)
                              for s in [ 'x', 'dx' ] ]
                           for j in range( self.num_peaks_per_feature[i] ) ]
                         for i in range( self.num_features) ]
                        
            if all( [ os.path.exists( path ) for path in np.array( peak_paths ).flatten() ] ) :

                for i in range( self.num_features ) :
                    for j in range( self.num_peaks_per_feature[i] ) :

                        peak = np.fromfile( peak_paths[i][j][0] ).reshape( self.dimensions )
                        peak_delta = np.fromfile( peak_paths[i][j][1] ).reshape( self.dimensions )
                        peak_grids[i][j] = meas.meas( peak, peak_delta )
                        
                return peak_grids

            else:
                print( 'INFO: requested to read peak and peak_delta from files, but they aren\'t there. constructing them now...' )


                

        # construct the array from the db if the files don't exist
        # yet, or we requested a fresh fetch.

        disconnect_conn_when_done = 0
        if self.conn is None:
            self.connect()
            disconnect_conn_when_done = 1

            
        # populate the grid 
        for i in range(3):
            for j in range(2):
                peak_grids[i][j] =  meas.meas.empty( self.dimensions )
        
        cursor = self.conn.cursor()
        cursor.execute( 'SELECT * FROM ' + _tablename )

        rownum = 0
        total_iterations = ( self.dimensions[0] * self.dimensions[1]
                             * self.num_features ) 
        
        for row in cursor:

            rownum += 1
            
            estimate_time_left( rownum, total_iterations, num_updates = 100 ) 

            x = row['x']
            y = row['y'] 
            feature = row[ 'fit_id' ]

            print( (x,y,feature)) 
            
            # if fit didn't converge then all the peak values for that
            # feature are assigned np.nan
            
            if not row[ 'successful_fit' ]: 
                for i in range( self.num_peaks_per_feature[ feature ]  ):
                    peak_grids[ feature ][ i ][ x,y ] = meas.nan
                continue

            
            model = _from_bin( row['model'] )


            # throw out the data if we didn't estimate errorbars.
            if not model.errorbars :
                for i in range( self.num_peaks_per_feature[ feature ]  ):
                    peak_grids[ feature ][ i ][ x,y ] = meas.nan
                continue
                
            # do a monte carlo simulation for each peak in this feature.
                
            peaks = spec.estimate_alpha_peakpos( model, plot = 0 )

            # check if the fit used the alternate shape
            if len( peaks ) != self.num_peaks_per_feature[ feature ] :
                    
                # todo: this doesn't work in the general case.
                # this should depend on the specified alternate shape.
                peaks = meas.append( meas.nan, peaks ) 
                    
                    
            # populate the corresponding entry of the grid.
            for i in range( self.num_peaks_per_feature[ feature ] ):
                peak_grids[ feature ][ i ][ x,y ] = peaks[i]

        # reset the counter 
        estimate_time_left( 0, 0, reset = 1 ) 
                
        if disconnect_conn_when_done:
            self.disconnect()


        # write to a file if requested. note that we have already constructed
        # all the file names.
        if read_from_file :
            
            for i in range( self.num_features ) :
                for j in range( self.num_peaks_per_feature[i] ) :
                    
                    peak_grids[i][j].x.tofile(  peak_paths[i][j][0] )
                    peak_grids[i][j].dx.tofile( peak_paths[i][j][1] )
                            
        return peak_grids    
    





    



    
# return binary dump that can be
# inserted in the DB as a BLOB
def _to_bin( data ): 
    return sqlite3.Binary( dill.dumps( data, protocol = 2 ) )




# read binary as written in the _to_bin format
def _from_bin( bin ):
    return dill.loads( bin )
    

