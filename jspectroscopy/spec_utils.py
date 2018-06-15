import peakdetect.peakdetect as peakdetect
import numpy as np


def get_n_peak_positions( n, data ):

    output_peak_positions = [0] * n
    
    # peakdetect returns 2 tuples: positions and counts of peaks
    peak_positions = peakdetect.peakdetect( data, lookahead=10 )[0]

    # it is possible that not all n peaks are found. in that case
    # we will still populate our_peaks with the ones that were
    # found as it may still be useful. 
    num_peaks_found = min( n, len( peak_positions ) )

        
    # now find the 5 largest and sort by x position.
    # indices is the indices of the peaks as found in data
    indices = np.argpartition( [ z[1] for z in peak_positions ],
                               -num_peaks_found )[ -num_peaks_found : ]
    
    output_peak_positions[:] = [ np.asscalar( peak_positions[z][0] )
                                 for z in sorted( indices ) ]

    
    # return number of peaks detected.
    return output_peak_positions

