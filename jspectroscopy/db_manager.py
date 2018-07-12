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

from jutils import meas as meas

import jutils 

import sys
import jspectroscopy as spec 

import matplotlib.cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt 
import datetime

import dill # allow lambda function pickle
import scipy.optimize
from scipy.stats import chi2



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

    
    def __init__( self, name, storage_path, dets_used = None, dimensions = None,
                  peak_types = None, constrain_det_params = None ) :
        
        self.conn = None
        self.name = name

        create_new_db = ( peak_types is not None ) and ( dimensions is not None )

        if storage_path == '' :
            storage_path = './'

        storage_path += '/' + name + '/'


        # create storage directories
        self.storage_path = storage_path
        self.heatmap_dir_path = storage_path + 'heatmaps/'
        self.strip_plots_path = storage_path + 'strip_plots/'
        
        for path in [ self.storage_path, self.heatmap_dir_path, self.strip_plots_path ] :
            if not os.path.exists( path ) :
                os.mkdir( path ) 
                
        # create the database if it doesn't exist
        self.db_path = storage_path + name + '.db'
        
        if not os.path.exists( self.db_path ) :
            print( self.db_path ) 
            if create_new_db :
                self.is_empty = 1 
                self.create( dets_used, dimensions, peak_types, constrain_det_params )
            else :
                print( 'ERROR: attempted to load db, but the path does not exist' )
                sys.exit(1)
        
        else : 
            self.is_empty = 0
            self.connect()
            self.read_metadata()

        
        for det in self.dets_used :
            tmp = self.heatmap_dir_path + '/' + str( det ) 
            if not os.path.exists( tmp ) :
                os.mkdir( tmp ) 

            
        # this is where all the extracted fit parameters are stored 
        self.param_paths = {}

        self.peak_indep_param = { 'mu' : 0, 'A' : 0, 'sigma' : 0,
                                  'eta' : 1, 'tau1' : 1, 'tau2' : 1,
                                  'sigma.calibrated' : 0,
                                  'tau1.calibrated'  : 1,
                                  'tau2.calibrated' : 1,
                                  'peaks' : 0,
                                  'peaks.calibrated' : 0 }

        self.fit_params = [ 'mu', 'tau1', 'tau2', 'A', 'eta', 'sigma' ];

        self.calibrated_fit_params = [ x + '.calibrated' for x in
                                       [ 'tau1', 'tau2', 'sigma', 'peaks' ] ]
        
        self.extra_values = [ 'success', 'hitmap', 'all_fit_params', 'peaks',
                              'calibration', 'peaks.calibrated' ]

        for param_str in self.fit_params + self.extra_values + self.calibrated_fit_params :

            # param_dir_path = storage_path + '/' + param_str

            self.param_paths[ param_str ] = ( storage_path + self.name + '.'
                                              + param_str + '.dill' ) 

        self.detnum_to_detidx_dict = dict( zip( self.dets_used, range( self.num_dets ) ) )
        self.detidx_to_detnum_dict = dict( zip( range( self.num_dets ), self.dets_used ) )

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

        if not os.path.exists( self.db_path ) :
            print( 'INFO: creating DB and schema for ' + self.db_path + '...' )
            
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
        # self.name = metadata[ 'name' ]

        self.num_groups = len( self.peak_types ) 
        self.num_peaks_per_group = [ len(x) for x in self.peak_types ]
        self.num_dets = len( self.dets_used ) 

        self.num_records = ( self.num_dets * self.xdim * self.ydim
                             * self.num_groups ) 

        
        
    # check if this db has been created. 
    def exists( self ):
        return os.path.exists( self.db_path ) 

    

    
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

        self.conn = sqlite3.connect( self.db_path ) 

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

        # print( result )

        sys.exit(1)

        # print( spec_result ) 

        return spec_result 




        



    def delete( self ):
        
        if os.path.exists( self.db_path ):
            
            # ans = raw_input( 'PROMPT: delete ' + filename + ', are you sure (y/n) ?  ' )
            
            # if( ans == 'y' ):

            os.remove( self.db_path ) 
            print( 'INFO: deleted db.' )         
            return 1

            print( 'INFO: did not delete db.' )

            return 0

        return 1
    


    
    # peak_indep_param: 0 if the container should have peak indices,
    # 1 if the container shouldn't have peak indices (i.e. just det and
    # group indices 
    
    def create_empty_data_container( self, peak_indep_param ) :

        if not peak_indep_param : 

            ret = [ [ 0 ] * self.num_peaks_per_group[i]
                    for i in range( self.num_groups ) ]
                    

            for i in range( self.num_groups ) :
                for j in range( self.num_peaks_per_group[i] ) : 
                    ret[i][j] = meas.zeros(( self.num_dets, self.xdim, self.ydim ))

        else :
            ret = [ [ 0 ] for i in range( self.num_groups ) ]

            for i in range( self.num_groups ) :
                ret[i] = meas.zeros(( self.num_dets, self.xdim, self.ydim ))
                    
        return ret
                


    
    # mode = output format
    # 0 = np ndarray
    # 1 = csv
    # 2 = both
    def write_fit_params( self, param_str, mode ) :
#    def write_mu_values( self ) :
        print( param_str ) 

        peak_indep = self.peak_indep_param[ param_str ] 
        values = self.create_empty_data_container( peak_indep )
        
        # init the data structure storing the mu values
        if param_str == 'mu' :
            param_getter = lambda spec, j : spec.get_mu(j)

        elif param_str == 'tau1' :
            param_getter = lambda spec, j : spec.get_tau1(j)

        elif param_str == 'tau2' :        
            param_getter = lambda spec, j : spec.get_tau2(j)

        elif param_str == 'sigma' :        
            param_getter = lambda spec, j : spec.get_sigma(j)
            
        elif param_str == 'A' :        
            param_getter = lambda spec, j : spec.get_A(j)
            
        elif param_str == 'eta' : 
            param_getter = lambda spec, j : spec.get_eta(j)
                        
        else :
            print( 'ERROR: param not found: ' + param_str )
            sys.exit(1)

        cursor = self.conn.cursor()
        cursor.execute( 'SELECT * FROM ' + _tablename )

        spec_fitters = [ spec.spectrum_fitter( self.peak_types[i],
                                               constrain_det_params = { 'a' : 1 } )
                         for i in range( self.num_groups ) ]

        dets_used_dict = dict( zip( self.dets_used, range( self.num_dets ) ) )
        
        for row in cursor:

            d = dets_used_dict[ row[ 'detnum' ] ]
            x = row['x']
            y = row['y'] 
            group_num = row[ 'group_num' ]
            success = row[ 'success' ]
            
            if success == 1 :
                params_result = _from_bin( row[ 'params_result' ] )
                params_result_errors = _from_bin( row[ 'params_result_errors' ] )
                
                spectrum = spec_fitters[ group_num ] 

                spectrum.set_params_from_params_array( params_result, errors = 0)
                spectrum.set_params_from_params_array( params_result_errors, errors = 1)
                    
                if not self.peak_indep_param[ param_str ] : 
                    for j in range( self.num_peaks_per_group[ group_num ] ) :
                        values[ group_num ][j][d,x,y] = param_getter( spectrum, j )

                else :
                    values[ group_num ][d,x,y] = param_getter( spectrum, 0 )
                                            
                        
            else :
                if not self.peak_indep_param[ param_str ] : 
                    for j in range( self.num_peaks_per_group[ group_num ] ) :
                        values[ group_num ][j][d,x,y] = meas.nan
                else :
                    values[ group_num ][d,x,y] = meas.nan


        if mode == 0 :
            with open( self.param_paths[ param_str ], 'wb' ) as f :
                dill.dump( values, f )

        elif mode == 1 :

            raise NotImplementedError
            
            # outpath = self.param_dir_path + param_str + '.csv'
            # with open( outpath, 'w' ) as f :
            #     values.write_csv( f )


            
    
    def write_peak_values( self, estimate_time = 0, dets = None ) :

        peak_path = self.param_paths[ 'peaks' ] 
        
        # init the data structure storing the mu values
        if estimate_time : 
            time_estimator = jutils.time_estimator( self.num_records, 20 ) 
        
        peak_values = self.create_empty_data_container(0)

        if dets is None :
            dets = self.dets_used 
                    
        cursor = self.conn.cursor()
        cursor.execute( 'SELECT * FROM ' + _tablename )

        spec_fitters = [ spec.spectrum_fitter( self.peak_types[i],
                                               constrain_det_params = { 'a' : 1 } )
                         for i in range( self.num_groups ) ]

        for row in cursor:

            time_estimator.update() 

            x = row['x']
            y = row['y'] 
            detnum = row[ 'detnum' ]
            d = self.detnum_to_detidx_dict[ detnum ]
            group_num = row[ 'group_num' ]
            success = row[ 'success' ]

            # print( d ) 
            
            if success > 0 and detnum in dets : 
                
                cov = _from_bin( row[ 'cov' ] ) 

                params_result = _from_bin( row[ 'params_result' ] )
                params_result_errors = _from_bin( row[ 'params_result_errors' ] )
                
                spectrum = spec_fitters[ group_num ] 

                spectrum.set_params_from_params_array( params_result, errors = 0)
                spectrum.set_params_from_params_array( params_result_errors, errors = 1)
                spectrum.cov = cov 
                
                for j in range( self.num_peaks_per_group[ group_num ] ) :
                                        
                    if spectrum.cov is not None : 
                        alpha_peak = spec.find_alpha_peak( spectrum, j, 500, 0 ) 

                        if alpha_peak is None :
                            peak_values[ group_num ][j][d,x,y] = meas.nan

                        else : 
                            peak_values[ group_num ][j][d,x,y] = alpha_peak 

                    else :
                        peak_values[ group_num ][j][d,x,y] = meas.nan
                        
            else :  
                for j in range( self.num_peaks_per_group[ group_num ] ) :
                    peak_values[ group_num ][j][d,x,y] = meas.nan 
                        
        with open( peak_path, 'wb' ) as f :
            dill.dump( peak_values, f )


            
    
    def write_all_fit_params( self, mode = 0 ) :
        for param_str in self.fit_params : 
            self.write_fit_params( param_str, mode )
            
        

            
    def read_params( self, param_str ) :

        path = self.param_paths[ param_str ] 
        
        if not os.path.exists( path ) :
            print( 'ERROR in db_manager.read_values(): path not found: %s'
                   % path ) 
            sys.exit(1) 

        with open( path, 'rb' ) as f :
            values = dill.load( f )

        return values 


    def read_all_fit_params( self ) :

        fit_params = {}
        
        for param_str in self.fit_params :
            fit_params[ param_str ] = self.read_params( param_str )

        return fit_params 


    
    def package_all_fit_params( self ) :
        all_fit_params = self.read_all_fit_params()
        with open( self.param_paths[ 'all_fit_params' ], 'wb' ) as f :
            dill.dump( all_fit_params, f )


            
    def read_packaged_fit_params( self ) :

        with open( self.param_paths[ 'all_fit_params' ], 'rb' ) as f :
            all_fit_params = dill.load( f )

        return all_fit_params 


    
    def calibrate_pixels( self, peaks_to_use, actual_energies ) :

        mu_values = self.read_params( 'mu' )

        flattened_actual_energies = np.array( [ a for b in range( len( actual_energies ) )
                                                for a in actual_energies[b] ] )
        
        total_num_peaks = len( flattened_actual_energies ) 
        
        # timer = jutils.time_estimator( self.num_dets * self.xdim * self.ydim, 10 ) 

        calibrations = meas.zeros( ( 2, self.num_dets, self.xdim, self.ydim ) )
        calibrations[:] = meas.nan
        
        for det in range( self.num_dets ) :

            mu_flattened = meas.empty( ( total_num_peaks, self.xdim, self.ydim ) )

            k = 0 
            for i in range( self.num_groups ) :
                for j in peaks_to_use[i] :
                    # print( det ) 
                    # print( mu_values.shape ) 
                    mu_flattened[k] = mu_values[i][j][det]
                    k += 1

            for x in range( self.xdim ) :
                for y in range( self.ydim ) :

                    # timer.update()


                    # skip data with any nan-entries
                    if np.sum( np.isnan( mu_flattened[:,x,y].x ) ) > 0 :
                        continue

                    params_guess = [ 1.0, 0.0 ] 
                    
                    ret = scipy.optimize.leastsq( self._linear_resid,
                                                  params_guess,
                                                  args = ( flattened_actual_energies,
                                                           mu_flattened[:,x,y].x,
                                                           mu_flattened[:,x,y].dx ),
                                                  full_output = 1 )

                    params_result, cov, info, msg, status = ret
        
                    success = ( status >= 1 and status <= 4
                                and ( cov is not None ) )

#                    print( success )

                    if status :

                        chisqr = np.sum( info['fvec']**2 )
                        nfree = len( flattened_actual_energies ) - len( params_result ) 
                        redchisqr = chisqr / nfree
                        cov = cov
            
                        params_result_errors = np.sqrt( np.diag( cov ) * redchisqr )
                        # calibrations[ det, x, y ] = params_result[:] 

                        pvalue = 1 - chi2.cdf( chisqr, nfree )
                        
                        # print() 
                        # print( mu_flattened[:,x,y].x ) 
                        # print( params_result ) 
                        # print( params_result_errors )
                        # print( redchisqr )

                        # print( pvalue )
                        if pvalue > 0.05 :
                            calibrations[ :, det, x, y ] = meas.meas( params_result ,
                                                                   params_result_errors )

        with open( self.param_paths[ 'calibration' ], 'wb' ) as f :
            dill.dump( calibrations, f )

                            # print( 

                        

                            
    def write_calibrated_params( self ) :

        calibration = self.read_params( 'calibration' ) 

        for param_str in [ 'tau1', 'tau2', 'sigma' ] :

            data = self.read_params( param_str )
            
            if not self.peak_indep_param[ param_str ] :
                calibrated = self.create_empty_data_container( 0 ) 
                for i in range( self.num_groups ) :
                    for j in range( self.num_peaks_per_group[i] ) :
                        calibrated[i][j] = data[i][j] / calibration[0]

            else :
                calibrated = self.create_empty_data_container(1)         
                for i in range( self.num_groups ) :
                    calibrated[i] = data[i] / calibration[0]
            
            outname = param_str + '.calibrated'            
            
            with open( self.param_paths[ outname ], 'wb' ) as f :
                dill.dump( calibrated, f )

        # do peaks separately because they need the b offset
        data = self.read_params( 'peaks' )
        calibrated = self.create_empty_data_container( 0 ) 
        for i in range( self.num_groups ) :
            for j in range( self.num_peaks_per_group[i] ) :
                calibrated[i][j] = ( data[i][j] - calibration[1].x ) / calibration[0].x

        with open( self.param_paths[ 'peaks.calibrated' ], 'wb' ) as f :
            dill.dump( calibrated, f )
            

            
    
    def plot_all_params( self, source_names = None, secant_matrices = None ) :

        max_peaks = max( self.num_peaks_per_group ) 

        cmap = matplotlib.cm.jet
        
        # plot heatmaps to get an idea of how the parameters look
        for param_str in [ 'tau1', 'tau2', 'sigma', 'A', 'eta', 'mu',
                           'peaks', 'peaks.calibrated',
                           'tau1.calibrated', 'tau2.calibrated', 'sigma.calibrated'] :

            print( param_str ) 

            data = self.read_params( param_str )
                
            for d in range( self.num_dets ) :

                det = self.detidx_to_detnum_dict[ d ]

                if not self.peak_indep_param[ param_str ] : 

                    f, heatmap_axarr = plt.subplots( self.num_groups, max_peaks,
                                             figsize = ( 12, 8 ) )
                    
                    for i in range( self.num_groups ) :
                        for j in range( max_peaks ) : 

                            if j < self.num_peaks_per_group[i] : 
                                im = heatmap_axarr[i,j].imshow( data[i][j][d].x, cmap = cmap )
                                divider = make_axes_locatable( heatmap_axarr[i,j] )
                                cax = divider.append_axes("right", size="5%", pad=0.05)
                                f.colorbar(im, cax=cax)

                                if source_names is not None :
                                    heatmap_axarr[i,j].set_title( source_names[i] + ': Peak %d'
                                                          % j )
                                
                            else :
                                heatmap_axarr[i,j].axis( 'off' )

                else :
                    f, heatmap_axarr = plt.subplots( 1, self.num_groups, figsize = (10,12) )
                    for i in range( self.num_groups ) :
                        
                        if source_names is not None :
                            heatmap_axarr[i].set_title( source_names[i] )

                        im = heatmap_axarr[i].imshow( data[i][d].x, cmap = cmap )
                        divider = make_axes_locatable( heatmap_axarr[i] )
                        cax = divider.append_axes("right", size="5%", pad=0.05)
                        f.colorbar(im, cax=cax)
                                
                f.suptitle( 'Det %d: %s' % ( det, param_str ) )
                f.suptitle( 'Det %d: %s' % ( det, param_str ) )
                
                
                plt.savefig( self.heatmap_dir_path + '%d/%s.png'
                             % ( det, param_str ) )

                plt.close( 'all' ) 
                
        if secant_matrices is None :
            return

        self.plot_strips( source_names, secant_matrices ) 




    # plot individual fstrips vs sec theta 
    def plot_strips( self, source_names, secant_matrices ) :

        max_peaks = max( self.num_peaks_per_group ) 
        
        for param_str in [ 'tau1', 'tau2', 'sigma', 'A', 'eta', 'mu',
                           'peaks', 'peaks.calibrated',
                           'tau1.calibrated', 'tau2.calibrated', 'sigma.calibrated'] :

            print( param_str ) 

            data = self.read_params( param_str )

            for x in range( self.xdim ) : 

                for d in range( self.num_dets ) :

                    no_data = 1 

                    detnum = self.detidx_to_detnum_dict[ d ]

                    if not self.peak_indep_param[ param_str ] : 

                        f, axarr = plt.subplots( self.num_groups, max_peaks,
                                                         figsize = ( 12,8 ) )

                        f.subplots_adjust( wspace = 0.5, hspace = 0.5 )

                        for i in range( self.num_groups ) :
                            for j in range( max_peaks ) :

                                if j < self.num_peaks_per_group[i] : 


                                    nan_count = np.count_nonzero(
                                        np.isnan( secant_matrices[i][d][x] )
                                        | np.isnan( data[i][j][d,x].x ) ) 
                                    
                                    if nan_count < self.xdim :
                                        no_data = 0

                                    axarr[i,j].errorbar( secant_matrices[i][d][x],
                                                         data[i][j][d,x].x,
                                                         data[i][j][d,x].dx,
                                                         ls = 'none' ) 

                                    if source_names is not None :
                                        axarr[i,j].set_title( source_names[i] + ': Peak %d'
                                                              % j )

                                else :
                                    axarr[i,j].axis( 'off' )

                    else :
                                                
                        f, axarr = plt.subplots( 1, self.num_groups, figsize = (12,8) )
                        for i in range( self.num_groups ) :

                            # print( secant_matrices[i][d][x] )
                            # print( data[i][d,x].x ) 
                            
                            if source_names is not None :
                                axarr[i].set_title( source_names[i] )

                            nan_count = np.count_nonzero(
                                np.isnan( secant_matrices[i][d][x] )
                                | np.isnan( data[i][d,x].x ) ) 
                            
                            if nan_count < self.xdim :
                                no_data = 0 
                                 
                            axarr[i].errorbar( secant_matrices[i][d][x],
                                               data[i][d,x].x,
                                               data[i][d,x].dx,
                                               ls = 'none' ) 

                    f.suptitle( 'Det %d: %s' % ( detnum, param_str ) )

                    outdir = self.strip_plots_path + str(detnum) + '/' + param_str + '/'

                    # if not os.path.exists( outdir ) :
                    os.makedirs( outdir, exist_ok = 1 )

                    if not no_data : 
                        plt.savefig( outdir + '%d.%s.%d.png' % ( detnum, param_str, x ) )

                    plt.close( 'all' ) 
        

                             
    # todo                 
    def gen_csv( self, path ) :
        pass 



    def _linear_fitfunc( self, params, x ) :
        return params[0] * x + params[1]

    def _linear_resid( self, params, x, y, yerr ) :
        return ( y - self._linear_fitfunc( params, x ) ) / yerr




    
    def compute_fit_success_rate( self ) :

        query = 'SELECT * FROM fits'

        cursor = self.conn.cursor()

        cursor.execute( query ) 

        total_count = 0
        success_count = 0 
        
        for row in cursor.fetchall() :

            row = dict( row )

            success = row[ 'success' ]

            if success == 1 :
                success_count += success

            total_count += 1


        # print( success_count )
        # print( total_count )
        
        return success_count / total_count 
            


    def plot_success_pixels( self ) :
        query = 'SELECT * FROM fits'

        cursor = self.conn.cursor()

        cursor.execute( query ) 

        total_count = 0
        success_count = 0

        all_values = np.zeros( ( self.num_groups, self.num_dets, 32, 32 ) ) 
        
        for row in cursor.fetchall() :

            row = dict( row )

            success = row[ 'success' ]
            g = row['group_num']
            d = row['detnum'] - 1 
            x = row['x']
            y = row['y']

            # print( g,d,x,y)
            if success == 1 :
                all_values[g,d,x,y] = 1
            
        f, axarr = plt.subplots( self.num_groups, self.num_dets )

        for g in range( self.num_groups ) :
            for d in range( self.num_dets ) :
                axarr[g,d].imshow( all_values[g,d] )
                
        plt.show()
        

    

    



    
# return binary dump that can be
# inserted in the DB as a BLOB
def _to_bin( data ): 
    return sqlite3.Binary( dill.dumps( data, protocol = 2 ) )




# read binary as written in the _to_bin format
def _from_bin( bin ):
    return dill.loads( bin )
    

