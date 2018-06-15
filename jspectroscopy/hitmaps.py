

# similar to spectrum_fitter: peak detection used, but
# no other information required. find peaks, check that peaks are valid,
# then integrate over the given range



import numpy as np 
import jutils 
import matplotlib.pyplot as plt
import matplotlib.colors
import dill 
import os
import sys
import jspectroscopy.spec_utils as utils 






def compute_hitmaps( dimensions,
                     data_retriever,
                     group_ranges, 
                     num_peaks_to_detect, primary_peak_detector,         
                     # peak_sizes_guesses, peak_width_guesses, det_params_guesses,
                     # peak_mu_offset,
                     # fit_acceptor = None,
                     # params_shuffler = None,
                     # rel_plot_bounds = None,
                     # logscale = 1,
                     # time_estimator = None,
                     # print_output = 0,
                     # dets_used = None,
                     filter_data = 0,
                     save_path = None,
                     reset = 0,
                     image_path = None,
                     plot = 0,
                     debug_coords = None,
                     rel_plot_bounds = None ) :
    
    # load saved data if it exists 
    if save_path is not None and not reset and debug_coords is None :
        if os.path.exists( save_path ) :
            with open( save_path, 'rb' ) as f :
                hitmaps = dill.load( f )

            # if filter_data :
                
                
            # print( hitmaps ) 
                
            if plot :
                plot_hitmaps( hitmaps ) 
            
            return hitmaps 

    num_maps = len( group_ranges )

    hitmaps = np.zeros( ( num_maps, * dimensions ) )

    if debug_coords :

        xaxis, histo, dy = data_retriever( * debug_coords )
        
        # find peaks and check if valid detection
        peaks = np.asarray( utils.get_n_peak_positions( num_peaks_to_detect, histo ) ) 
        primary_peaks = primary_peak_detector( peaks, histo )
        
        # handle invalid peak detect 
        if primary_peaks is None :
            print( 'Unable to find primary peaks.' )

        ax = plt.axes()

        ax.set_xlim( primary_peaks[0] + rel_plot_bounds[0],
                  primary_peaks[-1] + rel_plot_bounds[1] )
        
        ax.semilogy( xaxis, histo, c = 'k' )

        for i in range( num_maps ) :
            ax.axvline( x = primary_peaks[i], c = 'r' )
            tmp = primary_peaks[i] + np.array( group_ranges[i] )
            print( tmp ) 
            ax.plot( tmp, [10]*2, c = 'r' )

        plt.show() 
        
        return None 

    
    for x in range( dimensions[0] ) :
        for y in range( dimensions[1] ) :

            xaxis, histo, dy = data_retriever( x, y )
         
            # find peaks and check if valid detection
            peaks = np.asarray( utils.get_n_peak_positions( num_peaks_to_detect, histo ) ) 
            primary_peaks = primary_peak_detector( peaks, histo )

            # handle invalid peak detect 
            if primary_peaks is None :
                hitmaps[:, x, y] = np.nan
                continue

            # else count within group 
            else :
                for i in range( num_maps ) :                    
                    hitmaps[i,x,y] = np.sum( histo[ primary_peaks[i] + group_ranges[i][0] :
                                                     primary_peaks[i] + group_ranges[i][1] ] )
            
        
    # write the data
    if save_path is not None :
        with open( save_path, 'wb' ) as f :
            dill.dump( hitmaps, f )

    if plot :
        plot_hitmaps( hitmaps ) 
        
    return hitmaps 







def plot_hitmaps( hitmaps ) :

    num_maps = len( hitmaps )

    f, axarr = plt.subplots( 1, num_maps ) 
    
    for i in range( num_maps ) :

        axarr[i].imshow( hitmaps[i], norm = matplotlib.colors.Normalize() )

    plt.show() 






def plot_peakdetect_and_group_range( histo, peaks ) :
    pass 
