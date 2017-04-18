from scipy.signal import butter, lfilter, decimate, hilbert
from numpy import square, sqrt, mean, abs
from math import floor


"""
Python implementation of the DEMON algorithm for estimating the
envelope of an amplitude modulated noise. Two versions are given,
the square law and hilbert transform DEMON algorithms. Depending 
on the use case either algorithm may be more desirable. Both expect
an input array of the original data, and filter parameters. Both
return an output array containing the estimated envelope of the 
input data.

author: Alex Pollara

Algorithm is described in:
Pollara, A., Sutin, A., & Salloum, H. (2016). 
Improvement of the Detection of Envelope Modulation on Noise (DEMON) 
and its application to small boats. In OCEANS 2016 MTS-IEEE Monterey
IEEE.
    
"""


def square_law(x, cutoff=1000.0, high=30000, low=20000, fs=200000):
    """
    :param x: numpy.ndarray
    :param cutoff: float
    :param high: float
    :param low: float
    :param fs: float
    :return: numpy.ndarray
    """

    # check that parameters meet bandwidth requirements
    if (high+low)/2 <= 2*(high-low):
        raise Exception("Error, band width exceeds pass band center frequency")
    else:

        # Bandpass filter parameters
        nyq = .5 * fs  # band limit of signal Hz

        # Passband limits as a fraction of signal band limit
        high /= nyq
        low /= nyq
        order = 3

        # Butterworth bandpass filter coefficients
        # noinspection PyTupleAssignmentBalance
        b, a = butter(order, [low, high], btype='band')

        # filter signal
        x = lfilter(b, a, x)

        # square signal
        x = square(x)

        # calculate decimation rate
        n = int(floor(fs / (cutoff * 2)))
        print n
        # decimate signal by n using a low pass filter
        x = decimate(x, n, ftype='fir')

        # square root of signal
        x = sqrt(x)

        # subtract mean
        x = x - mean(x)

        return x


def hilbert_detector(x, cutoff=1000.0, high=30000, low=20000, fs=200000):
    """
    :param x: numpy.ndarray
    :param cutoff: float
    :param high: float
    :param low: float
    :param fs: float
    :return: numpy.ndarray
    """
    # Bandpass filter parameters
    nyq = .5 * fs  # band limit of signal Hz

    # Passband limits as a fraction of signal band limit
    high /= nyq
    low /= nyq
    order = 3

    # Butterworth bandpass filter coefficients
    # noinspection PyTupleAssignmentBalance
    b, a = butter(order, [low, high], btype='band')

    # filter signal
    x = lfilter(b, a, x)

    # hilbert transform of signal
    x = hilbert(x)

    # absolute value of signal
    x = abs(x)

    # calculate decimation rate
    n = int(floor(fs / (cutoff * 2)))

    # decimate signal by n using a low pass filter
    x = decimate(x, n, ftype='fir')

    # square root of signal
    x = sqrt(x)

    # subtract mean
    x = x - mean(x)

    return x
