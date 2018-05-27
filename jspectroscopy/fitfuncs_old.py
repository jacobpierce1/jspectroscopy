
# from functools import partial 


# global 'vars' that can be used for rapidly
# applying fitfuncs.

# alpha_fitfunc = fitfunc( 
# beta_fitfunc =

# gamma_fitfunc is synonymous with gaussian.
# gamma_fitfunc = 



    
    
# # class for all the fitfuncs. for maximal efficiency
# # and compabibility with scipy functions, they all
# # take array parameters as input. a function can also
# # be provided to convert the array to a dict
# # and if this is specified, another func to convert
# # such a dict to an array for future use. fitfunc
# # must as its first argument 

# class fitfunc( object ):

#     def __init__( self, fitfunc,
#                   params_array_to_dict_func = None,
#                   params_dict_to_array_func = None ):


#         if not callable( fitfunc ):
#             raise ValueError( 'fitfunc must be a function' )


#         # verify that f is a function that accepts 2 params.
#         if f.__code__.co_argcount != 2:
#             raise ValueError( 'fitfunc must take 2 parameters: ' +
#                               '( params_array, input_array )' )
        
        
#         self._fitfunc = fitfunc


#         # define the parameter conversion functions:
#         if params_array_to_dict_func is not None:

#             if not callable( params_array_to_dict_func ):
#                 raise ValueError( 'params_array_to_dict_func must be a function' )
        
#             self._array_to_dict_func = params_array_to_dict_func

            
#         if params_dict_to_array_func is not None:

#             if not callable( params_dict_to_array_func ):
#                 raise ValueError( 'params_dict_to_array_func must be a function' )

#             self._dict_to_array_func = params_dict_to_array_func


#     def dict_to_array( self, params_array ):
#         return self._dict_to_array_func( params_array )


#     def array_to_dict( self, params_array ):
#         return self._array_to_dict_func( params_array )
    

#     def apply( self, pf, input_array ):
#         return self._fitfunc( pf, input_array ) 




    
    
# # in this case you supply all args as fit params, even the det ones which should be fixed. 
# def fitfunc_alpha_free_det_params( p, x ): 

#     return alpha_fit( p[0],p[1],p[2],p[3],p[4],p[5], x ) 

# #    return np.array( map( lambda z: alpha_fit( *p, z ), x ) )
#     #return np.array( map( lambda z: alpha_fit( p[0],p[1],p[2],p[3],p[4],p[5],z ), x ) )
#     # return alpha_fit( p[0],p[1],p[2],p[3],p[4],p[5] , x )  # expand array as function args 
    



# return an inline of fitfunc_n_alpha_peaks given the number of peaks to fit.
# the first 4 params are the detector params, and the last 2*n are (A1, u1, ... )
# previously fitfunc_n_alpha_peaks_abstract
def sum_n_fitfuncs( fitfunc, n ):
    return lambda x, **params: fitfunc( n, x, **params )
    # return partial( fitfunc, n )




# # fit n alpha peaks given vector p and array x 
# # p format: sigma, eta, tau1, tau2, A1, mu1, ..., A_n, mu_n
# def fitfunc_n_alpha_peaks( n, p, x ):

#     ret = 0

#     for i in range( n ):
#         ret += alpha_fit( p[4+2*i],p[4+2*i+1],p[0],p[1],p[2],p[3], x )

#     return ret





def fitfunc_n_alpha_peaks_free_det_params( n, x, **params ) :

    ret = 0
    for i in range( n ):
        ret += alpha_fit( params[ 'A' + str(i) ], params[ 'mu' + str(i) ],
                          params[ 'sigma' + str(i) ], params[ 'eta' + str(i) ], 
                          params['tau1' + str(i)], params['tau2' + str(i)], x )

    return ret





# construct dict for input to fitfunc_n_alpha_peaks.

def construct_n_alpha_peaks_params_free_det_params(
        sigma, eta, tau1, tau2, A, mu ):

    n = len( A )

    params = {}
    param_fit_bounds = {}

    for i in range( n ):
        
        params[ 'mu' + str(i) ] = mu[i]
        params[ 'A' + str(i) ] = A[i]
        params[ 'sigma' + str(i) ] = sigma[i]
        params[ 'tau1' + str(i) ] = tau1[i]
        params[ 'tau2' + str(i) ] = tau2[i]
        params[ 'eta' + str(i) ] = eta[i]
        
        param_fit_bounds[ 'mu' + str(i) ] = [ 0, None ]
        param_fit_bounds[ 'A' + str(i) ] = [ 0, None ] 
        param_fit_bounds[ 'sigma' + str(i) ] = [ 0, None ] 
        param_fit_bounds[ 'tau1' + str(i) ] = [ 0, None ] 
        param_fit_bounds[ 'tau2' + str(i) ] = [ 0, None ] 
        param_fit_bounds[ 'eta' + str(i) ] = [ 0, None ] 

    return ( params, param_fit_bounds ) 





# # port the parameters of the model to a meas.meas
# def get_alpha_params_dict_free_det_params( model ): 

#     params = model.params

#     npeaks = ( len( params ) - 4 ) // 2  

#     ret = { 'sigma' : 0, 'tau1' : 0, 'tau2' : 0, 'eta' : 0 } 

#     for key in ret:
#         ret[ key ] = meas.meas( params[key].value, params[key].stderr )

#     for key_base in [ 'mu', 'A' ]:

#         vals = np.empty( npeaks )
#         deltas = np.empty( npeaks ) 
            
#         for i in range( npeaks ):

#             param = params[ key_base + str(i) ] 
#             vals[i] = param.value
#             deltas[i] = param.stderr

#         ret[ key_base ] = meas.meas( vals, deltas ) 
            

#         #    print( ret ) 
#     return ret








# return an lmfit model constructed from n alpha
# fit functions. keep sigma, eta, tau1, tau2 the same
# in all of them.

def fitfunc_n_alpha_peaks( n, x, **params ):
    
    ret = 0
    for i in range( n ):
        ret += alpha_fit( params[ 'A' + str(i) ], params[ 'mu' + str(i) ],
                          params[ 'sigma' ], params[ 'eta' ], 
                          params['tau1'], params['tau2'], x )

    return ret





# construct dict for input to fitfunc_n_alpha_peaks.

def construct_n_alpha_peaks_params(
        sigma, eta, tau1, tau2, A_array, mu_array ):

    n = len( mu_array )

    param_fit_bounds = { 'sigma' : [ 0, None ],
                         'tau1' : [ 0, None ],
                         'tau2' : [ 0, None ],
                         'eta' : [ 0, 1 ] }
    
    params =  { 'sigma' : sigma, 'eta' : eta,
                'tau1' : tau1, 'tau2' : tau2 }

    for i in range( n ):

        # keys for the 2 dicts 
        mu = 'mu' + str(i)
        A = 'A' + str(i) 
        
        params[ mu ] = mu_array[i]
        params[ A ] = A_array[i]

        param_fit_bounds[ mu ] = [ 0, None ]
        param_fit_bounds[ A ] = [ 0, None ] 

    return ( params, param_fit_bounds ) 





# port the parameters of the model to a meas.meas
def get_alpha_params_dict( model ): 

    params = model.params

    npeaks = ( len( params ) - 4 ) // 2  

    ret = { 'sigma' : 0, 'tau1' : 0, 'tau2' : 0, 'eta' : 0 } 

    for key in ret:
        ret[ key ] = meas.meas( params[key].value, params[key].stderr )

    for key_base in [ 'mu', 'A' ]:

        vals = np.empty( npeaks )
        deltas = np.empty( npeaks ) 
            
        for i in range( npeaks ):

            param = params[ key_base + str(i) ] 
            vals[i] = param.value
            deltas[i] = param.stderr

        ret[ key_base ] = meas.meas( vals, deltas ) 
            

        #    print( ret ) 
    return ret






# estimate the peak position for each peak if they were isolated.

def estimate_alpha_peakpos( model, num_iterations = 1000, plot=0 ) :
    
    params = get_alpha_params_dict( model )

    A = params['A']
    
    npeaks = len( A )
    
    peakpos_results = meas.meas.empty( npeaks ) 

    for peaknum in range( npeaks ) :

        params_array = meas.meas.from_array( np.array( [ A[peaknum],
                                                         params['mu'][peaknum],
                                                         params['sigma'],
                                                         params['eta'],
                                                         params['tau1'],
                                                         params['tau2'] ] ) )
        
        # if the starting params_array is invalid, then don't
        # do the simulation on it.
        
        if not alpha_model_valid_params_array_predicate( params_array ) :
            peakpos_results[peaknum] = meas.nan
            continue


        # otherwise proceed with the simulation.
        
        peakpos_arr = np.empty( num_iterations, dtype=np.float64 )


        # print( params_array ) 
    
        for i in range(num_iterations):
            
            # keep picking p until the amplitude is positive
            while 1:

                # randomize all params but don't bother with A,
                # since the peak location doesn't depend on A.
                
                current_params = np.empty( len( params_array ) )
                current_params[0] = params_array[0].x
                current_params[1:] = np.random.normal( params_array.x[1:],
                                                        params_array.dx[1:] )
                
                if alpha_model_valid_params_array_predicate( current_params ) :
                    break
                
            current_inverted_f = lambda x_: 0 - alpha_fit_array( current_params, x_ )
                
            result = scipy.optimize.fmin( current_inverted_f, current_params[1], disp=0 )
            
            peakpos_arr[i] = result
            


        
        sim_result = meas.meas( peakpos_arr.mean(), peakpos_arr.std() )

        # print( sim_result )

        peakpos_results[peaknum] = sim_result

        
        if plot :
            ax = plt.axes() 
            ax.set_title( 'Estimated Peak Channel Distribution', fontsize = 20 )
            ax.set_xlabel( 'Peak Channel', fontsize = 18 )
            ax.set_ylabel( 'Counts', fontsize = 18 )
            ax = plt.axes()
            ax.hist( peakpos_arr, bins = 20 ) 
            plt.show() 

            
    return peakpos_results








# return 1 if the alpha_model has valid parameters and 0 otherwise.

def alpha_model_valid_params_predicate( alpha_params_dict ) :


    for mu in alpha_params_dict['mu'] :
        if mu.x < 0:
            return 0

    for A in alpha_params_dict['A'] :
        if A.x < 0 or A.dx > A.x :
            return 0 
    
    for i in range(1,3) :
        tau = alpha_params_dict['tau' + str(i) ]
        if tau.x < 0 or tau.dx > tau.x :
            return 0

    eta = alpha_params_dict['eta']
    if eta.x < 0 or eta.x > 1 or eta.dx > eta.x :
        return 0

    sigma = alpha_params_dict['sigma']
    if sigma.x < 0 or sigma.dx > sigma.x : 
        return 0

    return 1 







# use for increased efficiency

def alpha_model_valid_params_array_predicate( alpha_params_array ) :

    if meas.ismeas( alpha_params_array ) :

        A, mu, sigma, eta, tau1, tau2 = [ x for x in alpha_params_array ] 
        
        if A.x < 0 :
            return 0
        
        if mu.x < 0 :
            return 0

        if sigma.x < 0 or sigma.dx > sigma.x :
            return 0
        
        if eta.x < 0 or eta.x > 1 or eta.dx > eta.x :
            return 0 
        
        if tau1.x < 0 or tau1.dx > tau1.x :
            return 0
        
        if tau2.x < 0 or tau2.dx > tau2.x :
            return 0

    else:
    
        A, mu, sigma, eta, tau1, tau2 = alpha_params_array
        
        if A < 0 :
            return 0
        
        if mu < 0 :
            return 0

        if sigma < 0 :
            return 0
        
        if eta < 0 or eta > 1 :
            return 0 
        
        if tau1 < 0 :
            return 0
        
        if tau2 < 0 :
            return 0
        
    return 1 

    






# TODO: DO THIS

# # input: a function that takes array of parameters and a scalar
# # variable x, same as input of optimize.least_sq; pf and pferr,
# # obtained from jacob_least_squares; peakpos_guess, estimate of the
# # peak positions; number of iterations to perform.

# # behavior: assume that pferr are standard deviations of a normal
# # distribution of which pf values are the means; do a monte carlo
# # simulation in which an array is formed with values from those normal
# # distributions. then find the maximum of the function and add it to a
# # list.  return: peakpos (average), peakpos_delta (std of mean),
# # peakval (function at peakpos), peakval_delta (estimated using 2nd
# # order taylor expansion of f; first order normally works, but in this
# # case we know that f'(x) = 0 at the max so it will give 0.

# def estimate_peakpos( spec_model, num_iterations=1000 ):

#     peakpos_arr = np.empty( num_iterations, dtype=np.float64 )
#     # print( 'p: ' + str(p) ) 
    
#     for i in range(num_iterations):

#         # keep picking p until the amplitude is positive
#         while 1:
#             current_p = np.random.normal( p, p_delta )

#             # break out of the while if certain entries are
#             # not physical. TODO: abstract this.
#             if current_p[4] > 0 and current_p[1] < 1:
#                 break

#         # now we construct a function from the random p
#         current_inverted_f = lambda x_: 0 - f( current_p, x_ )  
#         result = scipy.optimize.fmin( current_inverted_f, peakpos_guess, disp=0 )

#         peakpos_arr[i] = result
        
#     return peakpos_arr




# # fit n alpha peaks given vector p and array x 
# # p format: sigma, tau, A1, mu1, ..., A_n, mu_n
# def fitfunc_n_alpha_peaks_eta1( n, p, x ):

#     ret = 0
    
#     for i in range( n ):
#         ret += alpha_fit_eta1( p[2+2*i],p[2+2*i+1],p[0],p[1], x )

#     return ret


#     # return np.sum( [ map( lambda y: alpha_fit_eta1(p[2+2*i],p[2+2*i+1],p[0],p[1], y), x )  for i in range(n) ], axis=0 )  



    
# # same as alpha_fit but with eta fixed at one. this means we are taking only one of the terms. 
# # experimented with this after noting that the fits seem indep. of eta, which can cause trouble.
# def alpha_fit_eta1( A, mu, sigma, tau, x ):
#     logtmp = (x-mu)/tau + sigma**2.0/(2*tau**2.0) + np.log( special.erfc( (x-mu)/sigma + sigma/tau) / np.sqrt(2) ) 
#     return A/(2.0*tau) * np.exp( logtmp )





# wrapper allows for compact abstractions of the different fitfuncs.

def alpha_fit_array( params_array, x ) :
    return alpha_fit( * params_array, x ) 
