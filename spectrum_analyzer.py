import jspectroscopy as spec
import numpy as np 
import matplotlib.pyplot as plt
import scipy.optimize
import libjacob.jmeas as meas


# input: alpha spectrum with one peak or more peaks
# output: tuple containing location of the peak and
# confidence interval

def find_alpha_peak( alpha_spectrum, peaknum = 0,
                     num_iterations = 1000, plot = 0, max_failed_samples = 1000 ) :
    
    # print( 'alpha_spectrum params:' + str( alpha_spectrum._construct_params_array() ) )

    single_peak_spectrum = alpha_spectrum.subpeak( peaknum )
    
    mean_params = np.copy( single_peak_spectrum._construct_params_array() )

    peakpos_arr = np.empty( num_iterations ) 

    failed_samples = 0 
    
    for i in np.arange( num_iterations ) :

        status = 0

        while not status : 
        
            sampled_params = np.random.multivariate_normal( mean_params,
                                                            single_peak_spectrum.cov )

            sampled_params[0] = 1
            
            single_peak_spectrum.set_params_from_params_array( sampled_params ) 
            
            if single_peak_spectrum.check_valid_fit_params() : 
                break

            elif failed_samples >= max_failed_samples :
                return None 

            failed_samples += 1 
            
        # now we have valid params. compute the min of this function.
        
        # print( 'mean params: ' + str( mean_params ) )
        # print( 'sampled_params: ' + str( sampled_params )  ) 

        current_inverted_f = lambda x_: 0 - single_peak_spectrum.eval( x_ ) 
        
        result = scipy.optimize.fmin( current_inverted_f, sampled_params[1] - 4, xtol = 0.01, disp=0 )
        
        peakpos_arr[i] = result

    m = np.mean( peakpos_arr ) 
    s = np.std( peakpos_arr  )

    print( m, s ) 

    
    if plot :
        ax = plt.axes() 
        ax.set_title( 'Estimated Peak Channel Distribution', fontsize = 20 )
        ax.set_xlabel( 'Peak Channel', fontsize = 18 )
        ax.set_ylabel( 'Counts', fontsize = 18 )
        ax = plt.axes()
        ax.hist( peakpos_arr, bins = 20 ) 
        plt.show() 

            
    return meas.meas( m, s ) 












    # peakpos_results = meas.meas.empty( npeaks ) 

    # for peaknum in range( npeaks ) :

    #     params_array = meas.meas.from_array( np.array( [ A[peaknum],
    #                                                      params['mu'][peaknum],
    #                                                      params['sigma'],
    #                                                      params['eta'],
    #                                                      params['tau1'],
    #                                                      params['tau2'] ] ) )
        
    #     # if the starting params_array is invalid, then don't
    #     # do the simulation on it.
        
    #     if not alpha_model_valid_params_array_predicate( params_array ) :
    #         peakpos_results[peaknum] = meas.nan
    #         continue


    #     # otherwise proceed with the simulation.
        
    #     peakpos_arr = np.empty( num_iterations, dtype=np.float64 )


    #     # print( params_array ) 
    
    #     for i in range(num_iterations):
            
    #         # keep picking p until the amplitude is positive
    #         while 1:

    #             # randomize all params but don't bother with A,
    #             # since the peak location doesn't depend on A.
                
    #             current_params = np.empty( len( params_array ) )
    #             current_params[0] = params_array[0].x
    #             current_params[1:] = np.random.normal( params_array.x[1:],
    #                                                     params_array.dx[1:] )
                
    #             if alpha_model_valid_params_array_predicate( current_params ) :
    #                 break
                
    #         current_inverted_f = lambda x_: 0 - alpha_fit_array( current_params, x_ )
                
    #         result = scipy.optimize.fmin( current_inverted_f, current_params[1], disp=0 )
            
    #         peakpos_arr[i] = result
            


        
    #     sim_result = meas.meas( peakpos_arr.mean(), peakpos_arr.std() )

    #     # print( sim_result )

    #     peakpos_results[peaknum] = sim_result

        
