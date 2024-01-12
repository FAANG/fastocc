from typing import Iterator
import numpy as np


def nonzero_intervals(x: np.ndarray) -> Iterator[tuple]:
    """
    Calculate continuous regions in a numpy array, returning
    a generator of (start, stop) tuples.

    e.g. [0, 0, 1, 1, 1, 0, 1, 1, 0] => (2,5); (6,8)

    NB If the interval starts or ends with 1, this is handled, as a result the indexing
    is basically BED-like, ie (0, 1) represents an interval containing just the first
    item. This should be very rare for our exact use case, it means that the peak caller
    has missed a nucleosome free region.

    e.g. [1, 0, 0, 1, 1] => (0, 1); (3, 5)
    """
    d = np.diff(x)
    idx = d.nonzero()[0] + 1

    if x[0]:
        idx = np.r_[0, idx]

    if x[-1]:
        idx = np.r_[idx, x.size]

    for start, stop in np.reshape(idx, (-1, 2)):
        yield start, stop


# should be a rolling window function and a separate calculation
def rolling_window(arr: np.ndarray, window_size: int):
    """
    Sum counts over overlapping 110bp (or 11 bin) windows. The resulting array is 5 bins
    shorter on either side than the input.
    """
    # create matrix representing 11-bin windows into counts
    index = np.arange(window_size)[None, :] + np.arange(arr.shape[0] - window_size + 1)[:, None]

    # sum over windows
    return np.sum(arr[index], 1)


def smooth(x: np.ndarray, gaussian: np.ndarray, gaussian_norm):
    """
    Smoothes input 1D numpy array using hardcoded gaussian convolution. Output is
    normalised using a hardcoded constant derived from convolving the kernel with a
    vector of 1s. NaNs in x will result in errors. Returns an array with 5 fewer values
    on each side so that only 'valid' results are returned.
    """
    
    smoothed = np.convolve(gaussian, x, mode = 'valid')

    return smoothed / gaussian_norm


def gaussian(M, std):
    """
    Adapted from scipy.signal.gaussian
    """
    n = np.arange(0, M) - (M - 1.0) / 2.0
    sig2 = 2 * std * std
    w = np.exp(-n ** 2 / sig2)

    return w
