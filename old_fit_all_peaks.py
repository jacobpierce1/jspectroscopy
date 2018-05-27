
# the purpose of this script is to bulid off of fit_5_peaks to run
# through all the data. this poses a challenge because it is difficult
# to verify that all the fits actually worked. the strategy is to
# display two images of 16 plots (32 total) for each of the 32 pixel
# rows. alarming chisq values ( > 2 ) are reported, along with
# instances in which a fit fails to converge. all the pf's are written
# to a file for later processing.

# the first step is to identify the peaks and use this to make guesses
# for p0, the parameter guess for the least squares fit.

# date: 9.06.17
# author: jacob pierce


# default config. do not adjust here, make a dedicated function that sets 
# config parameters in an abnormal situtation (e.g. debugging). any mistake in
# the db can be very costly in time, either through repairing db which will 
# probably required a specialized function or through rerunnnig the program.



SHOW_HISTO = 0
PRINT_PEAK_POSITIONS = 0
PRINT_PF = 0
PRINT_PFERR = 0
PRINT_FIT_STATUS = 0
PRINT_FILE_NAMES = 0
BREAK_AFTER_1_PLOT = 0
SAVE_DATA = 1
UPDATE_DB = 1
PRINT_FIT_ATTEMPT = 0




# my files 
import libjacob.jmath as jmath
import libjacob.jpyplot as jplt
import libjacob.jutils as jutils


import jspectroscopy as spec

# from peakdetect import peakdetect 
import deadlayer_helpers.sql_db_manager as dbmgr





# from lmfit import Model, Parameters 




## includes 
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"

import time 
import sys
import sqlite3
import os
from enum import Enum 
from scipy.stats import chi2




class spectrum_fitter_inputs( object ) :

    def __init__( self, dimensions, data_fetcher, 
                  group_ranges, peak_locations,
                  num_peaks_to_detect, primary_peak_detector,
                  peak_sizes_guesses, det_params_guesses, peak_mu_offset,
                  peak_position_tolerance = None, db_path = None,
                  name = None,
                  fit_acceptor = None ) :

        self.dimensions = dimensions 

        # for 2D case: a function where the input is (i,j)
        # and the output is a histogram of counts vs channel

        # if not callable( data_fetcher ) :
        #     print( '''data_fetcher must be a function that returns a histogram for 
        #     pixel (i,j), which are the 2 args of this function. currently data_fetcher
        #     is not callable''' )

        self.data_fetcher = data_fetcher

        self.num_groups = len( group_ranges )

        if len(peak_locations) != self.num_groups : 

            print( '''error: inconsistent size of group_ranges, peak_structures, or peak_locations.
            they should all have length equal to group_ranges''' )
            sys.exit(0)

        self.peaks_per_group = [ len(peak_locations[i]) for i in range( self.num_groups ) ]
            
        self.group_ranges = group_ranges
        # self.peak_structures = peak_structures
        self.peak_locations = peak_locations
        # self.primary_peak_ids = primary_peak_ids

        self.num_peaks_to_detect = num_peaks_to_detect
        
        if callable( primary_peak_detector ) :
            self.primary_peak_detector = primary_peak_detector 
        else :
            sys.exit(0)


        self.peak_sizes_guesses = peak_sizes_guesses
        self.det_params_guesses = det_params_guesses
        self.peak_mu_offset = peak_mu_offset
        
        self.name = name
        
        # construct a DB in which to store the data if a name for the db
        # is supplied. 
        if db_path is not None :
            pass
        else:
            self.db = None



            

# NUM_PEAKS_PER_FEATURE = [ 2, 2, 2 ] 
# NUM_PEAKS = np.sum( NUM_PEAKS_PER_FEATURE )
# NUM_FEATURES = len( NUM_PEAKS_PER_FEATURE )




def make_model_and_apply_fit( npeaks, params_guess, x, y, dy ):
    pass


# these are the functions that attempt to modify p0_attempt and fit_bounds_attempt.
def default( p0, p0_attempt, fit_bounds, fit_bounds_attempt, npeaks ):
    pass 
    
def increase_left_fit_bound( p0, p0_attempt, fit_bounds, fit_bounds_attempt, npeaks ):
    fit_bounds_attempt[0] += 10
    
def increase_left_fit_bound_more( p0, p0_attempt, fit_bounds, fit_bounds_attempt, npeaks ):
    fit_bounds_attempt[0] += 20
    
def smaller_A( p0, p0_attempt, fit_bounds, fit_bounds_attempt, npeaks ):
    # print( npeaks ) 
    for i in range( npeaks ):
        p0_attempt['A' + str(i)] = 0.4 * p0[ 'A' + str(i) ]

def larger_A( p0, p0_attempt, fit_bounds, fit_bounds_attempt, npeaks ):
    # fit_bounds_attempt[1] += 5
    for i in range( npeaks ):
        p0_attempt[ 'A' + str(i) ] = 3 * p0[ 'A' + str(i) ]

def next_fit( p0, p0_attempt, fit_bounds, fit_bounds_attempt, npeaks ):
    pass



class fit_attempt_status( Enum ):
    SUCCESS = 1
    FAILURE = 2
    NO_UPDATE = 3



# this function creates another guess for p0 in the case of a failed fit based on
# observations about what is gonig wrong. the particular attempts we use depend on the peak id.
def next_fit_attempt( p0, p0_attempt, fit_bounds, fit_id,
                      fit_bounds_attempt, npeaks, attempt ):

    
    # these will be modified and used later as the initial guess.
    for key, val in p0.items():
        p0_attempt[ key ] = val     

    fit_bounds_attempt[:] = fit_bounds[:]    


        
    # array of functions that will be used to modify p0_attempt and fit_bounds_attempt
    if fit_id == 1:
        modifier_functions = [ default, increase_left_fit_bound, larger_A, next_fit ]

    else:
        modifier_functions = [ default, increase_left_fit_bound, 
                                    increase_left_fit_bound_more, smaller_A ]

        
    # this covers the case in which we have already used all the functions as recorded in the 
    # db, so we don't bother to try again. the only way to get out of this is reset the function
    # number in db or add another function (latter almost certainly what you need).
    # note that we don't update the db
    if( attempt >= len(modifier_functions ) ):
        return fit_attempt_status.NO_UPDATE

    
    # this checks if we had not previously attempted all functions, but none of them succeeded
    # on this run. in that case put it in the DB. 
    if attempt == len(modifier_functions):
        return fit_attempt_status.FAIL

    
    # otherwise we call the current attempt function
    # to generate a new fit parameters attempt.
    modifier_functions[ attempt ] ( p0, p0_attempt, fit_bounds, fit_bounds_attempt, npeaks )
    
    if PRINT_FIT_ATTEMPT:
        print( 'INFO: testing fit attempt ' + str(attempt) )
        print( 'p0_attempt = ' + str(p0_attempt) )
        print( 'fit_bounds_attempt = ' + str(fit_bounds_attempt) )
        
    return 1
    







# if db connection supplied, then we assume that the entry is supposed
# to be put in, it makes more sense for the check on whether the data
# is in the db to be performed before calling this function. when
# called: apply peak fit, make sure it worked, try again if not, put
# it in db if successful, and then plot. return 1 if successful fit
# obtained, otherwise 0. last_attempt is the last attempted + 1,
# i.e. it is the first fit that will be attempted.

def apply_peak_fit( ax, input_params, indices,
                    group_number,
                    xvals, yvals,
                    fit_bounds,
                    p0, npeaks, last_attempt = -1, 
                    params_bounds = None, reduc_chisq_all=None,
                    db=None, peak_detect=None ):
    
    
    fitfunc = spec.sum_n_fitfuncs( spec.fitfunc_n_alpha_peaks_free_det_params, npeaks )

    
    # these will be populated based on the current attempt number and
    # the given p0 and fit_bounds, which are used as the first
    # guesses. they are passed as parameters so they can be used after
    # without returning them.

    p0_attempt = {}
    fit_bounds_attempt = []

    
    # start attempt is the one after the most recently tested attempt.
    current_attempt = last_attempt + 1

    # keep trying new fit parameters till they are exhausted and all fits
    # fail, or if we have a convergent fit.
    
    while( 1 ):

        status = next_fit_attempt( p0, p0_attempt,
                                   fit_bounds, group_number, fit_bounds_attempt, 
                                   npeaks, current_attempt )

        # our fit attempts all failed in this case 
        if status == fit_attempt_status.FAILURE:

            
            if PRINT_FIT_STATUS:
                print( "INFO: fit failed, adding to DB " )
            
            if UPDATE_DB and db is not None:
                db.insert_fit_data( x, y, group_number,
                                    last_attempt=attempt )

        # in this case, the fits all failed but it was already in DB 
        if status == fit_attempt_status.NO_UPDATE:
            
            if PRINT_FIT_STATUS:
                print( "INFO: peak failed to converge, already detected in db.." )

            return 0


        # if those 2 cases are not reached, then we have a new fit attempt.
        # try doing least squares. 
        model = jmath.jleast_squares( xvals, yvals, np.sqrt( yvals ),
                                      p0_attempt, fitfunc,
                                      fit_bounds = fit_bounds_attempt,
                                      reduc_chisq_max = 2,
                                      params_bounds = params_bounds,
                                      successful_fit_predicate = successful_alpha_fit_predicate,
                                      pvalue = 0.05,
                                      print_results = 0 )
        
        # if done, do final checks on the fit. if successful, break out of
        # loop and proceed to add to DB / create plot

        if model is not None:
            if PRINT_FIT_STATUS:
                print( "INFO: success, breaking " )

            # take the data from the model and apply a peakdetect. make sure that
            # the peaks are the same. model.best_fit is the model evaluated on the
            # x data provided, which is xvals.
        
            if peak_detect is not None:

                # peaks predicted by model: 
                model_peaks = jmath.get_n_peak_positions( npeaks, model.best_fit )

                # peaks observed, corrected for the new start channel:
                original_peaks = np.array( peak_detect ) - fit_bounds_attempt[0]

                # print( 'model_peaks: ' + str( model_peaks ) )
                # print( 'peak_detect: ' + str( peak_detect ) )
                # print( 'original_peaks: ' + str( original_peaks ) )

                
                # try again if wrong number of peaks
                if len( model_peaks ) != npeaks :
                    current_attempt += 1
                    continue

                # now make sure the peaks are the same
                restart = 0 
                for a in range( npeaks ) :
                    if abs( original_peaks[a] - model_peaks[a] ) >= 2 :
                        restart = 1
                        break

                if restart :
                    current_attempt += 1
                    continue

            # if this is reached, then a fit converged which passed all the tests
            # now break and add to DB / create plot.
            break

            
        # otherwise try new fit params.
        else:
            if PRINT_FIT_STATUS:
                print( "WARNING: Unsuccessful fit, trying new parameters...\n\n" )
            current_attempt += 1
            continue


    # now we have broken out of the loop after a successful fit. unpack the results.
    reduc_chisq = model.redchi

    print('')
    print( model.params ) 

    # model.plot_fit( ax = ax, datafmt = '-r', numpoints = 100 * model.ndata )

    # dof = model.nfree
    # pf = model.params
    # pferr = model.         


    
    
    # write all relevant data to the db
    if UPDATE_DB and db is not None:
        db.insert_fit_data( x, y, group_number,
                            successful_fit = 1,
                            npeaks = npeaks, 
                            last_attempt = current_attempt,
                            params_guess = p0_attempt,
                            fit_bounds = fit_bounds_attempt, 
                            peak_guesses = peak_detect,
                            model = model ) 
        

    jplt.add_fit_to_plot( ax, xvals, fit_bounds_attempt, jmath.model_func( model ) )


    if reduc_chisq_all is not None:
        if PRINT_FIT_ATTEMPT :
            print( 'p value: ' + str( 1 - chi2.cdf( model.chisqr, model.nfree ) ) )
            print( model.params )
            
        reduc_chisq_all.append( reduc_chisq )


        
    return 1
    
    






# only allow a successful fit if we have sub-2 channel precision
# on the mu values.

def successful_alpha_fit_predicate( model ):

    params = model.params
    
    i = 0
    keys = params.keys()

    while( 1 ) :

        mu = 'mu' + str(i)
        
        if mu in keys: 
            if params[ mu ].stderr > 4:
                return 0
            i += 1 
            
        else:
            break
            
    return 1 









# main function, goes through all the data.
def fit_all_spectra( dbs ):

    # dimesions of detector
    totalx = 32
    totaly = 32

    # dimensions of output images: number of plots in each dimension
    dimx = 4
    dimy = 4

    a = -1

    for db in dbs :

        # count number of db's we have processed. 
        a += 1 
        print( db.name )
        

        # construct db if not already existing 
        if not db.exists():
            db.create()

        db.connect()

       
        # make dir for output images if it doesn't exist 
        current_outdir = '../../images/current_fit_images/' + db.name + '/'

        if not os.path.exists( current_outdir ):
            os.mkdir( current_outdir )


        logfile = current_outdir + 'log.txt'
        
        
        # these 2 loops loop over grids of dimx x dimy images. 
        for x in range( totalx ):
            
            # these loops create the 4x4 grids. 2 per strip row.
            for k in range( totaly // (dimx * dimy ) ):
                
                plt.clf()
                plt.close()
                f, axarr = plt.subplots( dimx, dimy, figsize=(20,10) )    
                
                for i in range(dimx):
                    for j in range(dimy):
                        
                        y = ( i * dimx + j ) + ( k * dimx * dimy )
                        
                        print( str( [x,y] ) )    
                        
                        # this estimates the time remaining for the program to terminate
                        jutils.estimate_time_left( x * totaly + y + ( a * totalx * totaly ), 
                                                   len( dbmgr.all_dbs ) * totalx * totaly,
                                                   num_updates=100 )
                        
                        # current pixel coords
                        # coords = ( x, y )
                        
                        # file containing the data
                        current_file = db.name + "_%d_%d.bin" % ( x, y )
                        
                        if PRINT_FILE_NAMES:
                            print( "INFO: processing file: " + current_file )
                            
                        current_file = ( "../../data/extracted_ttree_data/"
                                         + db.name + '/' + current_file )
                        
                        process_file( axarr[i,j], current_file, x, y,
                                      db = db, logfile = logfile ) # sql_conn = sql_conn, logfile = logfile)
                        
                        if( BREAK_AFTER_1_PLOT ):
                            plt.show()
                            return 1
                        
                jplt.saveplot_low_quality( current_outdir, db.name + '_x=%d.%d' % (x, k) )
                            
            
        db.disconnect()
                        # plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
                            # plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
                            
                            



                            

# dimensions, data_fetcher, 
#                   group_ranges, peak_locations,
#                   peak_position_tolerance = None, db_path = None ) :

                  
    
# this function is to be used for debugging the fits. once you have found a fit that 
# does not converge, add a new fit and 
def auto_fit_spectrum( histo, group_ranges, peak_locations,
                       num_peaks_to_detect, primary_peak_detector,
                       peak_sizes_guesses, det_params_guesses, peak_mu_offset,
                       peak_position_tolerance = None, db_path = None,
                       fit_acceptor = None, ax = None ) :
    
    # data_fetcher = lambda x: histo

    inputs = spectrum_fitter_inputs( None, None,
                                     group_ranges, peak_locations,
                                     num_peaks_to_detect, primary_peak_detector,
                                     peak_sizes_guesses, det_params_guesses, peak_mu_offset,
                                     peak_position_tolerance, db_path,
                                     fit_acceptor = None )

    # inputs = spectrum_fitter_inputs( None, data_fetcher, group_ranges ) 

    plt.figure(figsize=(10,12))
    
    ax = plt.axes() 
    
    process_spectrum( inputs, None, ax, histo ) 

    plt.show() 


    
    # if db not in dbmgr.all_dbs :
    #     print( 'ERROR: db_id given is not in the list of db ids. ' )
    #     return 0
    
    # set_globals_for_debugging( test_db )
    
    # current_file = db.name +  "_%d_%d.bin" % (x,y)
    
    # if PRINT_FILE_NAMES:
    #     print( "INFO: processing file: " + current_file )
    
    # current_file = ( "../../data/extracted_ttree_data/"
    #                  + db.name + '/'  + current_file )
        
    # ax = plt.axes()
    

    # if test_db:

    #     # construct db if not already existing 
    #     if not db.exists():
    #         db.create()

    #     db.connect()
    #     process_file( ax, current_file, x, y, db = db )
    #     db.disconnect()

    # else:
    #     process_file( ax, current_file, x, y, nice_format = 1 ) 

        
    # plt.show()

    # return 1
    
    
       




# print message to the console and write it to a log file.
def log_message( logfile, msg ):
    
    if log_message.first_msg or not os.path.exists( logfile ):
        log_message.first_msg = 0
        mode = 'w'
    else:
        mode = 'a'
        
    with open( logfile, mode ) as log:
        log.write( msg + '\n' )
    
log_message.first_msg = 1





# set debug globals, called by make_one_plot()
def set_globals_for_debugging( test_db = 0 ):

    if not test_db:
        global UPDATE_DB
        UPDATE_DB = 0
    
    global PRINT_PEAK_POSITIONS
    PRINT_PEAK_POSITIONS = 1
    
    global PRINT_FIT_ATTEMPT
    PRINT_FIT_ATTEMPT = 1
    
    global PRINT_FIT_STATUS
    PRINT_FIT_STATUS = 1

    global SHOW_HISTO
    SHOW_HISTO = 0
     
       
# make_one_plot( dbmgr.angled, 24, 30, test_db = 0 )
# fit_all_peaks( dbmgr.all_dbs )


