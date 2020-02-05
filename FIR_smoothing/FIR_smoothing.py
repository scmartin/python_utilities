import numpy as np
from scipy.optimize import least_squares


def FIRed_value(epsilon, A, window):
    window_size = window.shape[1]
    if window_size % 2 == 0:
        raise ValueError('Window size must be an odd number')
    return A*np.sum(np.exp(-((np.arange(window_size)-(window_size-1)/2)/epsilon)**2)*window,axis=1)

def residual(epsilon, A, N, data):
    epsilon = epsilon[0]
    residual_val = (FIRed_value(epsilon, A, np.roll(data[:,:2*N+1],1,axis=0)) 
                    + FIRed_value(epsilon, A, np.roll(data[:,:2*N+1],-1,axis=0)) 
                    - 2*FIRed_value(epsilon, A, data[:,:2*N+1]))
    return residual_val

def parameter_finder(fine_scale, N=25, e_0=5.0):
    """ Find the parameters for the FIR smoothing using gaussian smoothing 
    of every point of the form A*exp(-(k / \epsilon)^2) where k is the index
    of each point in the smoothing window relative to the central point (i.e.
    -N >= k >= N). This routine returns epsilon and A, where A is chosen such
    that \sum_{k} A*exp(-(k / \epsilon)^2) = 1.
    
    
    Arguments:
    fine_scale  -- raw histogram data to be smoothed
    N           -- number of data points on either side of central point
                   to include in smoothing (default 25)
    e_0         -- initial guess for epsilon (default 5.0)
    
    Returns:
    smoothed    -- (numpy array) smoothed version of the histogram
    parameters  -- (list) [epsilon, A] optimal value of epsilon
                   and A which gives proper normalization as described above
    
    For reference, see Laird and Davidchack, 1998."""

    rolled = np.array([np.roll(fine_scale,-i) for i in range(fine_scale.shape[0])])
    result = least_squares(residual, e_0, 
                       args=(1.0, N, rolled), method='lm')
    epsilon = result.x
    A = 1.0/np.sum(np.exp(-(np.arange(-N,N+1)/epsilon)**2))
    
    smoothed = np.roll(FIRed_value(epsilon, A, rolled[:,:2*N+1]), N+1)

    return smoothed, [epsilon, A]
    

if __name__ == '__main__':
    pass
