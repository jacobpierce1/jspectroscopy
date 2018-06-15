# this module contains a class for constructing fitfuncts
# with conversions between their input array params
# and dicts for interpretation / data analysis later on.


# TODO:
# * option to automatically construct a sum of fitfuncs
# * * suboption to hold some params constant when doing this

import matplotlib 
matplotlib.use('Agg')


import jspectroscopy.spec_utils as spec_utils




import numpy as np
np.warnings.filterwarnings('ignore')

import scipy.special as special
# from lmfit import Model

# import libjacob.jmeas as meas
from jutils import meas 

import matplotlib.pyplot as plt

import scipy.optimize
from scipy.stats import chi2

import os 

from jutils import time_estimator

# from libjacob.jpyplot import saveplot_low_quality

import copy


from enum import Enum



# alpha fit params :  (A, mu, sigma) * N, eta, tau1, tau2



class peak_types( Enum ) :
    a = 0
    b = 1 
    g = 2



_fitters_num_params = { 'a' : (3,3),
                        'b' : None,
                        'g' : None,
                        'cb' : None } 

     
# __alpha_num_peak_params = 2
# __alpha_num_det_params = 4







def xcut( x, y, newx_bounds, xsorted = 0 ):

    # if x is sorted, then find newx_bounds[0] from the left
    # and newx_bounds[1] from the right. should speed up the cut,
    # but has not been tested.

    if xsorted:
        
        left = np.searchsorted( x, newx_bounds[0] )
        right = np.searchsorted( x, newx_bounds[1], side = 'right' )
        return y[ left : right ]

    # otherwise find applicable indices by brute force.
    # should be avoided.

    return np.asarray( y[ (x >= newx_bounds[0]) & (x <= newx_bounds[1]) ] )












def _resid( f, params, x, y, dy ) :
    ret = ( f( params, x ) - y ) / dy

    # print( dy ) 
    # print( ret )
    
    return ret 






def alpha_fitter( peak_params, det_params, x ) :
    return alpha_fit( *peak_params, *det_params, x ) 


def beta_fitter( peak_params, det_params, x ) :
    raise NotImplementedError()


def gamma_fitter( peak_params, det_params, x ) :
    raise NotImplementedError()



# check if params are valid

def alpha_fit_params_checker( peak_params, det_params ) :

    # check mu and A 
    if np.any( peak_params < 0 ) :
        return 0

    # check sigma, eta, tau1, tau2
    if np.any( det_params < 0 ) :
        return 0

    # check on eta param
    if det_params[0] > 1 :
        return 0

    return 1 

    
def beta_fit_params_checker( peak_params, det_params, x ) :
    raise NotImplementedError()


def gamma_fit_params_checker( peak_params, det_params, x ) :
    raise NotImplementedError()



_available_fits = [ 'a', 'b', 'g' ]

_num_available_fits = len( _available_fits ) 


_all_fitters = { 'a' : alpha_fitter,
                 'b' : beta_fitter,
                 'g' : gamma_fitter } 


_all_fit_params_checkers = { 'a' : alpha_fit_params_checker,
                             'b' : beta_fit_params_checker,
                             'g' : gamma_fit_params_checker } 


class spectrum_fitter( object ) :


    def __init__( self, peak_types, constrain_det_params = None, group_num = None ) :

        if constrain_det_params is None :
            constrain_det_params = { 'a' : 0, 'b' : 0, 'g' : 0 } 

        self.constrain_det_params = copy.deepcopy( constrain_det_params )
            
        self.peak_types = copy.deepcopy( peak_types ) 

        self.num_peaks = len( peak_types )

        self.group_num = group_num
        
        # self.peak_type_indices = {}

        # for t in peak_types  :
        #     peak_type_indices[t] = [ i for i, x in enumerate( peak_types )
        #                              if x == t ]
            
        
        self.fitters = [ _all_fitters[t] for t in peak_types ] 

        self.params_array_peak_indices = np.empty( ( self.num_peaks, 2 ), dtype='int' ) 
        self.params_array_det_indices = np.empty( ( self.num_peaks, 2 ), dtype='int' )

        # set the indices of the relevant params for the flattened array

        peak_idx = 0 
        det_idx = 0

        constrained_det_param_indices = { 'a' : -1, 'b' : -1, 'g' : -1 }

        
        self.peak_params = [0] * self.num_peaks
        self.det_params = [0] * self.num_peaks
        self.peak_params_errors = [0] * self.num_peaks
        self.det_params_errors = [0] * self.num_peaks

        
        for i in np.arange( self.num_peaks ) :

            peak_type = self.peak_types[i]

            num_peak_params, num_det_params = _fitters_num_params[ peak_type ]

            # set these defaults 
            self.peak_params[i] = np.zeros( num_peak_params )
            self.det_params[i] = np.zeros( num_det_params )
            self.peak_params_errors[i] = np.zeros( num_peak_params )
            self.det_params_errors[i] = np.zeros( num_det_params ) 
            
            self.params_array_peak_indices[i] = np.array( [ peak_idx,
                                                            peak_idx + num_peak_params ] )
            peak_idx += num_peak_params

            if constrain_det_params[ peak_types[i] ] :

                constrained_det_param_idx = constrained_det_param_indices[ peak_type ]
                
                if constrained_det_param_idx == -1 :

                    self.params_array_det_indices[i] = np.array( [ det_idx,
                                                                   det_idx + num_det_params ] )

                    constrained_det_param_indices[ peak_type ] = det_idx
                    det_idx += num_det_params

                else : 
                    self.params_array_det_indices[i] = np.array( [ constrained_det_param_idx,
                                                                   constrained_det_param_idx
                                                                   + num_det_params ] )

            else : 
                self.params_array_det_indices[i] = ( det_idx, det_idx + num_det_params )
                det_idx += num_det_params
            
    
        self.params_array_det_indices += peak_idx
        self.params_array_size = peak_idx + det_idx
             

        # defaults
        self.success = -1 
        self.params_result = None
        self.params_result_errors = None
        self.redchisqr = -1
        self.cov = None



        

    def fit( self, x, y, dy, 
             peak_params_guess, det_params_guess, xbounds = None,
             fit_acceptor = None, params_shuffler = None,
             ax = None, plot_bounds = None, logscale = 1,
             print_output = 0 ) :

        
        self.peak_params_guess = copy.deepcopy( peak_params_guess )
        self.det_params_guess = copy.deepcopy( det_params_guess )
        self.xbounds = copy.copy( xbounds )
        self.plot_bounds = copy.copy( plot_bounds ) 

        if ax is not None:
            ax.plot( x, y, ls = 'steps', zorder = 1, c = 'k', linewidth = 0.1 ) 

            if( logscale ):
                ax.set_yscale('log')
                ax.set_ylim( bottom = 1 ) 
            
            if plot_bounds is not None :
                ax.set_xlim( * plot_bounds ) 
            
        self.fit_acceptor = fit_acceptor
        self.params_shuffler = params_shuffler 
    
        #self.ax = ax
        
        self.fit_attempted = 1

        
        if xbounds is not None :
            y = xcut( x, y, xbounds )
            dy = xcut( x, dy, xbounds )
            x = xcut( x, x, xbounds ) 

            
        self.fit_bounds = xbounds
        # self.x = x
        
        for peak_type, det_params in det_params_guess.items() : 

            # convert everything to an array
            det_params_guess_arr = np.asarray( det_params )
            det_params_guess[ peak_type ] = det_params_guess_arr
            
            # if self.constrain_det_params[ peak_type ] : 
            #     det_params_guess[ peak_type ] = det_params_guess_arr

            
        for i in np.arange( self.num_peaks ) :

            peak_type = self.peak_types[i]
            self.peak_params[i] = np.asarray( peak_params_guess[i] )
            self.det_params[i] = np.copy( det_params_guess[ peak_type ] )
            self.peak_params_errors[i] = np.zeros_like( self.peak_params[i] )  
            self.det_params_errors[i] = np.zeros_like( self.det_params[i] )  
            

        
            # if self.constrain_det_params[ peak_type ] : 
            #     self.det_params[i] = det_params_guess[ peak_type ]  
            # else :
            #     self.det_params[i] = np.copy( det_params_guess[ peak_type ] ) 
                
        params_array = self._construct_params_array()
       
        # print( params_array ) 

        # print( 'det_params: ' + str( self.det_params ) )
        # print( 'peak_params: ' + str( self.peak_params ) )
        # print( 'params_array: ' + str( params_array ) ) 

        # # call scipy.curve_fit with appropriate input
        # params_final, cov = scipy.optimize.curve_fit( self.__fit_eval, x, y,
        #                                               p0 = params_array,
        #                                               sigma = dy )

        objective = lambda _params, _x, _y, _dy : _resid( self.__fit_eval, _params,
                                                          _x, _y, _dy ) 

        
        ret = scipy.optimize.leastsq( objective, params_array, args = (x, y, dy),
                                      full_output = 1 )

        params_result, cov, info, msg, status = ret
        
        success = ( status >= 1 and status <= 4
                    and ( cov is not None ) )

        self.success = success

        # we think that the fit converged. now compute the chisqr, pvalue and
        # check anything else that could be fishy about the fit in an
        # external function

        if success :
            
            if print_output :
                print( '' ) 
                print( 'fit converged!' )
                print( self.xbounds ) 
                    

            if not self.check_valid_fit_params() :
                self.success = 0
                if print_output :
                    print( 'invalid fit params' )
                    print( [ '%.2f' % z for z in params_array ] )
                    print( [ '%.2f' % z for z in params_result ] )
                    
                return self 
            
            self.chisqr = np.sum( info['fvec']**2 )
            self.nfree = len( x ) - len( params_array ) 
            self.redchisqr = self.chisqr / self.nfree
            self.cov = cov

            params_result_errors = np.sqrt( np.diag( cov ) * self.redchisqr )

            self.set_params_from_params_array( params_result_errors, 1 )
            
            self.pvalue = 1 - chi2.cdf( self.chisqr, self.nfree )

            if print_output : 
                print( [ '%.2f' % z for z in params_array ] )
                print( [ '%.2f' % z for z in params_result ] )
                print( [ '%.2f' % z for z in params_result_errors ] )
                
                print( 'redchisqr: %.3f' % self.redchisqr )
                print( 'pvalue: %.3f' % self.pvalue )

            
            if fit_acceptor is not None :
                if not fit_acceptor( x, y, dy, self ) :
                    if print_output :
                        print( 'fit acceptor failed!' ) 
                    self.success = 0
                    return self

            self.params_result = params_result
            self.params_result_errors = params_result_errors



            # print(cov) 
                
            if self.success and ax is not None :
                ax.plot( x, self.eval( x ), c = 'r', zorder = 2 ) 


        else :
            if print_output :
                print( '' ) 
                print( 'fit failed to converge' ) 
                
        return self
        

                
    def eval( self, x ) :

        ret = 0
        for i in np.arange( self.num_peaks ) :
            ret += self.fitters[i]( self.peak_params[i], self.det_params[i], x ) 

        return ret

    

    
    def __fit_eval( self, params_array, x ) :

        self.set_params_from_params_array( params_array, 0 )

        return self.eval(x)

    
            
    def _construct_params_array( self ) :

        p = np.empty( self.params_array_size )

        for i in np.arange( self.num_peaks ) :

            # print( self.params_array_peak_indices[i] )
            # print( self.peak_params[i] ) 

            p[ slice( * self.params_array_peak_indices[i] ) ] = self.peak_params[i]
            p[ slice( * ( self.params_array_det_indices[i] ) ) ] = self.det_params[i]

        return p 


    

    def set_params_from_params_array( self, params_array, errors = 0 ) :

        if not errors : 
            peak_array = self.peak_params
            det_array = self.det_params

        else :
            peak_array = self.peak_params_errors
            det_array = self.det_params_errors
        
        for i in np.arange( self.num_peaks ) :

            # print( ''  )
            # print( 'num_peaks: ' + str( self.num_peaks ) ) 
            # print( params_array ) 
            # print( slice( * self.params_array_peak_indices[i] ) )
            # print( slice( * self.params_array_det_indices[i] ) )
            # print( peak_array[i] )
            # print( det_array[i] )
 
            peak_array[i][:] = params_array[ slice( * self.params_array_peak_indices[i] ) ]
            det_array[i][:] = params_array[ slice( * self.params_array_det_indices[i] ) ]
            
        
            

    def check_valid_fit_params( self ) :
        
        for i in range( self.num_peaks ) :
            if not _all_fit_params_checkers[ self.peak_types[i] ]( self.peak_params[i],
                                                                   self.det_params[i] ) :
                return 0
        
        return 1

    

    def plot( self, ax, c = 'r', **kw_args ) :
        ax.plot( self.x, self.eval( self.x ), c=c, **kw_args )



        

    # return a spectrum_fitter with one peak, equal to the peaknum
    # peak in self.
    def subpeak( self, peaknum ) :

        peak_type = self.peak_types[ peaknum ]
        
        single_peak = spectrum_fitter( peak_type )

        single_peak.num_peaks = 1 
        
        single_peak.peak_params = [ np.copy( self.peak_params[ peaknum ] ) ] 
        single_peak.peak_params_errors = [ np.copy( self.peak_params_errors[ peaknum ] ) ]
        single_peak.det_params = [ np.copy( self.det_params[ peaknum ] ) ]
        single_peak.det_params_errors = [ np.copy( self.det_params_errors[ peaknum ] ) ]
        
        num_peak_params, num_det_params = _fitters_num_params[ peak_type ]

        single_peak.params_array_peak_indices = [ [ 0, num_peak_params ] ]
        single_peak.params_array_det_indices = [ [ num_peak_params,
                                                   num_peak_params + num_det_params ] ]

        # single_peak.set_params_from_params_array( self._construct_params_array() ) 

        # single_peak._construct_params_array()

        single_peak.success = self.success

        if self.cov is not None :
            
            cov_indices = np.concatenate( [ np.arange( * self.params_array_peak_indices[ peaknum ],
                                                       dtype = int ),
                                            np.arange( * self.params_array_det_indices[ peaknum ],
                                                       dtype = int ) ] )
            
            single_peak.cov = np.copy( self.cov[ cov_indices ][ :, cov_indices ] )

        else :
            single_peak.cov = None
            
        return single_peak 




    
        
# reference: equation 10 in Bortels 1987  
# this function is meant to be applied to scalar x, not list

def alpha_fit( A, mu, sigma, eta, tau1, tau2, x ):

    # print( 'A: ' + str( A ) )
    # print( 'mu: ' + str( mu ) ) 
    # print( 'sigma: ' + str( sigma ) )
    # print( 'eta: ' + str( eta ) ) 
    # print( 'tau1: ' + str( tau1 ) )
    # print( 'tau2: ' + str( tau2 ) ) 

    
    # prevent overflow by computing logs and then exponentiating, at the expense of some
    # floating pt error. logtmpz is the log of the 2 analagous terms in the integrand.
    tauz = [tau1, tau2]

    logtmpz = [ (x-mu)/tau + sigma**2.0/(2*tau**2.0)
                + np.log( special.erfc( (x-mu)/sigma + sigma/tau) / np.sqrt(2) )
                for tau in tauz ]

    return ( A / 2.0 ) * ( (1-eta)/tau1 * np.exp(logtmpz[0])
                           + (eta/tau2) * np.exp(logtmpz[1])   ) 











def fit_spectrum( peak_types, x, y, dy, xbounds,
                  peak_params_guess, det_params_guess,
                  constrain_det_params = None,
                  fit_acceptor = None, params_shuffler = None,
                  ax = None, plot_bounds = None, logscale = 1,
                  group_num = None, print_output = 0 ) :

    spec_fitter = spectrum_fitter( peak_types, constrain_det_params, group_num )

    spec_fitter.fit( x, y, dy, peak_params_guess,
                     det_params_guess, xbounds,
                     fit_acceptor, params_shuffler,
                     ax, plot_bounds, logscale = logscale,
                     print_output = print_output ) 
    
    return spec_fitter










def auto_fit_spectrum( x, y, dy,
                       group_ranges, peak_locations,
                       num_peaks_to_detect, primary_peak_detector,
                       peak_sizes_guesses, peak_width_guesses,
                       det_params_guesses, peak_mu_offset,
                       fit_acceptor = None,
                       params_shuffler = None,
                       ax = None,
                       rel_plot_bounds = None,
                       logscale = 1, print_output = 0 ) :
    
    num_groups = len( group_ranges )

    if len(peak_locations) != num_groups : 

        print( '''ERROR: inconsistent size of group_ranges, peak_structures, or peak_locations.
        they should all have length equal to group_ranges''' )
        sys.exit(0)

    peaks_per_group = [ len(peak_locations[i]) for i in range( num_groups ) ]
    
    # find main peaks
    our_peaks = np.asarray( spec_utils.get_n_peak_positions( num_peaks_to_detect, y ) )

    # print( our_peaks ) 
    
    primary_peaks = primary_peak_detector( our_peaks, y )

    # self.primary_peaks = primary_peaks 
    
    if primary_peaks is None :
        return None 
    
    # print( primary_peaks )
    
    # print( primary_peaks ) 
    
    num_peaks_found = len( our_peaks )


    # do a check on peak values: energy differenc for the pairs should
    # be constant. no check necessary on the largest peak since it
    # always dominates, the only potential problem is really the
    # smallest peak in the lowest energy pair. this occurs after
    # plotting so that we can return with the plot alread made.
    # basically, we are saying that if we cannot determine where the
    # peak positions are to 0th order, we are not going to bother
    # fitting them since it will definitely not work if we misidentify
    # a peak position.
    
        
        
    # determine which fits to perform, 1 = fit must be attempted. by default if no 
    # conn is supplied we process all the fits.

    fit_attempts = [ -1 ] * num_groups  # first fit to try, which is the one we left off on.

    fits_to_perform = [ 1 ] * num_groups 

    reduc_chisq_all = [] # store all reduced chisq values.
    

    fit_bounds = [ np.array( group_ranges[a] ) + primary_peaks[a]
                   for a in range( num_groups )  ]

    if rel_plot_bounds is not None :
        plot_bounds = [ fit_bounds[0][0] + rel_plot_bounds[0],
                        fit_bounds[-1][1] + rel_plot_bounds[1] ] 

    else :
        plot_bounds = None
        
    # list of detected peaks, to be added to DB (not yet implemented) 
    # peak_detect = [ our_peaks[0:2], our_peaks[2:4], our_peaks[4:] ]

    ret = [0] * num_groups 
    
    # loop through the fits that were not in the db and add them if successful.
    for i in range( num_groups ):

        if fits_to_perform[i]:

            # print( fit_bounds[i] )
            # print( primary_peaks[i] ) 

            mu_array_guess = ( primary_peaks[i]
                               + np.array( peak_locations[i] )
                               + peak_mu_offset )
            
            peak_params_guess = [ [ peak_sizes_guesses[i][d], mu_array_guess[d],
                                    peak_width_guesses[i][d] ]
                                  for d in range( peaks_per_group[i] ) ]


            # print( peak_params_guesses )
            # print( det_params_guess )

                        
            spec_result = fit_spectrum( [ 'a' ] * peaks_per_group[i],
                                        x, y, dy, fit_bounds[i],
                                        peak_params_guess, det_params_guesses[i],
                                        constrain_det_params = { 'a' : 1 },
                                        params_shuffler = params_shuffler,
                                        fit_acceptor = fit_acceptor,
                                        ax = ax,
                                        plot_bounds = plot_bounds,
                                        logscale = logscale,
                                        group_num = i,
                                        print_output = print_output )

            ret[i] = spec_result 

    return ret






        
        



def auto_fit_many_spectra( spec_db, data_retriever,
                           image_path, image_dimensions, 
                           group_ranges, peak_locations,
                           num_peaks_to_detect, primary_peak_detector,
                           peak_sizes_guesses, peak_width_guesses, det_params_guesses,
                           peak_mu_offset,
                           fit_acceptor = None,
                           params_shuffler = None,
                           rel_plot_bounds = None,
                           logscale = 1,
                           time_estimator = None,
                           print_output = 0,
                           dets_used = None ) :
    
    # dimensions of output images: number of plots in each dimension

    if dets_used is None :
        dets_used = [ -1 ] 
    
    xdim = spec_db.xdim
    ydim = spec_db.ydim
    
    im_xdim = image_dimensions[0]
    im_ydim = image_dimensions[1]

    num_groups = len( group_ranges ) 
    
    if not os.path.exists( image_path ):
        os.mkdir( image_path )
        
    # these 2 loops loop over grids of dimx x dimy images.
    for d in dets_used : 
        for x in range( xdim ):

            # these loops create the 4x4 grids. 2 per strip row.
            for k in range( ydim // (im_xdim * im_ydim) ) :

                f, axarr = plt.subplots( im_xdim, im_ydim, figsize=(20,10) )    

                for i in range( im_xdim ):
                    for j in range( im_ydim ):

                        y = ( i * im_xdim + j ) + ( k * im_xdim * im_ydim )

                        # print( str( [x,y] ) )    

                        # this estimates the time remaining for the program to terminate
                        if time_estimator : 
                            time_estimator.update() 

                        tmp = data_retriever( d, x, y )

                        if tmp is not None :
                            xdata, ydata, dydata = tmp

                        else :
                            continue

                        spec_fits = auto_fit_spectrum( xdata, ydata, dydata,
                                                       group_ranges, peak_locations,
                                                       num_peaks_to_detect, primary_peak_detector,
                                                       peak_sizes_guesses, peak_width_guesses,
                                                       det_params_guesses,
                                                       peak_mu_offset,
                                                       fit_acceptor, params_shuffler,
                                                       axarr[i,j], rel_plot_bounds, logscale,
                                                       print_output )

                        for l in range( num_groups ) :
                            if spec_fits is not None : 
                                spec_db.insert_fit_data( d, x, y, l, spec_fits[l] ) 

                            else :
                                spec_db.insert_fit_data( d, x, y, l, None ) 

                plt.savefig( image_path + ( '%d_%d_%d' % (d,x,k) ) + '.png', format='png')
                # plt.clf()
                plt.close( f )


    # spec_db.disconnect()

    return 1 

    # plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
                            # plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
                            


    # # option to pass None as the DB 
    # if db is not None:

    #     for i in range( num_groups ):
            
    #         # extract
    #         db_data = db.read_fit_data( x, y, i )

    #         successful_fit = db_data[ 'successful_fit' ]

    #         # signal that this fit will have to be re-attempted.
    #         fits_to_perform[i] = not successful_fit

    #         # last attempt that we left off on 
    #         fit_attempts[i] = db_data[ 'last_attempt' ]
            
    #         if successful_fit:
    #             # fitfunc = spec.sum_n_fitfuncs( spec.fitfunc_n_alpha_peaks, NUM_PEAKS_PER_FEATURE[i] )
                
    #             model = db_data[ 'model' ] 
    #             jplt.add_fit_to_plot( ax, x, db_data[ 'fit_bounds' ],
    #                                   jmath.model_func( model ),
    #                                   logscale = 1 )

    #             # print( model.params ) 
                
    #             reduc_chisq_all.append( model.redchi )
                

    #     # no further processing required if all the fits converged.
    #     if not any( fits_to_perform ):
    #         _add_text( ax, x, y, reduc_chisq_all )
    #         return 1
                
         

    pass 
