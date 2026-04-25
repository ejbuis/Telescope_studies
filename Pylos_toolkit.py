from math import *
import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt
from scipy import signal
from scipy.optimize import minimize_scalar
import sys
# from my_styles import *
from scipy.io.wavfile import read

from numpy.random import seed
from numpy.random import shuffle

###########################################
# TOOLKIT FOR PROCESSING PYLOS NOISE DATA
###########################################


########################### load data  ###########################
def load_Pylos_data(filename, trace_start_s = 2325., trace_length_s = 180):
    """
    Load Pylos data from a WAV file.

    Parameters:
    filename : str
        The name of the WAV file to load.
    trace_start_s : float
        The start time of the trace in seconds.
    trace_length_s : float
        The length of the trace in seconds.

    Returns:
    t : ndarray
        The time vector of the loaded trace.
    y : ndarray
        The amplitude values of the loaded trace.
    fs : float
        The sampling frequency of the loaded trace.
    """

    time_trace    = read(filename)
    sampling_rate = time_trace[0]
    time_series   = time_trace[1]

    print("Sampling rate:         ", sampling_rate, "Hz")
    print("Length of time series: ", round(len(time_series)/sampling_rate/60), "minutes")

    trace_end   = trace_start_s + trace_length_s  # [s]

    i_trace_start = int(trace_start_s * sampling_rate)
    i_trace_end   = int(trace_end * sampling_rate)

    t = np.arange((i_trace_end - i_trace_start),dtype=float)
    t /= sampling_rate

    y = np.array(time_series[i_trace_start : i_trace_end],dtype=float)

    # Processing
    # y *= scaling # Remove scaling, as we will rescale later
    y -= np.mean(y) # Remove DC offset
    print("Length of output time series: ", round(len(y)/sampling_rate), "seconds")

    return t, y, sampling_rate, filename


########################### PSD values ###########################
def PSD_using_scipy(ampl, fs, nperseg=1024):
    """Calculate the Power Spectral Density (PSD) using scipy's Welch method.

    Parameters:
    ampl : array_like
        The amplitude values of the signal.
    fs : float
        The sampling frequency of the signal.
        
    Returns:
    f : ndarray
        Frequencies at which the PSD is computed.
    Pxx : ndarray
        Power spectral density of the signal.
    """
    return signal.welch(ampl, fs, nperseg=nperseg, scaling= 'density',window='hann')


def Knudsen_curve(seastate = 0, frequency=1000):
    """Calculate the Knudsen curve for a given seastate and frequency.
    The seastate should be an integer from 0 to 6, where 0 is calm and 6 is stormy.
    The frequency should be in Hz.
    Returns the Knudsen curve value in dB re 1 uPa^2/Hz.
    Valid for f > 1000 Hz
    """
    assert seastate < 7
    ss = [44.5, 50, 55., 61.5, 64.5, 66.5, 70]
    
    return -17.* np.log10(np.maximum(abs(frequency), 1000)/1000) + ss[seastate]


def Thermal_noise(frequency=1000):
    """
    Calculate the thermal noise level for a given frequency or array of frequencies.
    The frequency should be in Hz (scalar or ndarray).
    Returns the thermal noise level in dB re 1 uPa/√Hz.
    Valid for f >> 1 Hz, so clip to 100 Hz
    """
    return -15 + 20 * np.log10(np.maximum(abs(frequency), 1) / 1000)

def add_PSDs(*PSDs):
    if len(PSDs) == 1:
        return PSDs[0]
    else:
        return add_PSDs(10 * np.log10(10**(PSDs[0]/10) + 10**(add_PSDs(*PSDs[1:])/10)))


############################ Correction factors ###########################
def correction_factor_Square(SS):
    """Calculate the correction factor for a given seastate using the Square method.
    The seastate should be an integer from 0 to 6, where 0 is calm and 6 is stormy.
    Returns the correction factor in dB re 1 uPa^2/Hz.
    """
    assert SS < 7
    scale_factors = np.array([2.45e+02, 1.30e+02, 7.31e+01, 3.46e+01, 2.45e+01, 1.95e+01, 1.30e+01])
    return scale_factors[SS]


def correction_factor_Welch(SS):
    """Calculate the correction factor for a given seastate using the Welch method.
    The seastate should be an integer from 0 to 6, where 0 is calm and 6 is stormy.
    Returns the correction factor in dB re 1 uPa^2/Hz.
    """
    assert SS < 7
    scale_factors = np.array([4.62e+02, 2.45e+02, 1.38e+02, 6.53e+01, 4.62e+01, 3.67e+01, 2.45e+01])
    
    return scale_factors[SS]


############################ Generate subtraces from time series ###########################
def random_subtrace(times_series, length, sampling_rate=144e3):
    """
    Randomly select a subtrace of a given length from the time series.
    """
    start = np.random.randint(0, len(times_series) - length)
    if start + length > len(times_series):
        raise ValueError("The length of the subtrace exceeds the available data.")
    
    t = np.linspace(0, length / sampling_rate, length, endpoint=False)

    ts = times_series[start:start + length]
    ts = np.array(ts, dtype=float)
    
    ts -= np.mean(ts)

    return t, ts

def random_subtrace_ss(seastate, times_series, length, sampling_rate=144e3,method='Welch'):
    """
    Randomly select a subtrace of a given length from the time series
    Scaled to match the seastate level.
    The timeseries should be zero-mean

    The factors are computed from:      "201650198.180131085932.wav" [2325, 2335]s

    Outputs:
        t : ndarray [s]
            The time vector for the subtrace.
        p : ndarray [Pa]
            The scaled pressure subtrace.
    """

    if method=='Square':
        # scale_factors = np.array([1.18e-03, 6.26e-04, 3.52e-04, 1.64e-04, 1.17e-04, 9.34e-05, 6.40e-05]) # From 0-10 seconds (muPa)
        # scale_factors = np.array([2.45e-04, 1.31e-04, 7.20e-05, 3.64e-05, 2.62e-05, 1.95e-05, 1.24e-05]) # From 2325, 2335 seconds (muPa)
        scale_factor = correction_factor_Square(seastate)
    elif method=='Welch':
        scale_factor = correction_factor_Welch(seastate)
    else:
        raise ValueError("Method must be 'square' or 'welch'.")


    start = np.random.randint(0, len(times_series) - length)
    if start + length > len(times_series):
        raise ValueError("The length of the subtrace exceeds the available data.")
    
    assert seastate < 7

    t = np.linspace(0, length / sampling_rate, length, endpoint=False)

    p = np.array(times_series[start:start + length],dtype=float)
    p /= scale_factor  # Scale to SS

    return t, p


########################### Generate fit parameters for Knudsen curve ###########################
def time_domain_scaling_factor_from_SS(time_series, SS, sampling_rate=144e3, fmin=1000, fmax=25000):
    """
    Calculate the scaling factor for the time series based on the seastate.
    The time series should be zero-mean.
    The function uses the Knudsen curve to find the optimal scaling factor
    that minimizes the difference between the power spectral density of the
    time series and the Knudsen curve for the given seastate.

    Parameters:
    time_series : array_like
        The time series data to be processed.
    SS : int
        The seastate level (0 to 6).
    sampling_rate : float, optional
        The sampling rate of the time series (default is 144e3 Hz).
    fmin : float, optional
        The minimum frequency for fitting (default is 1000 Hz).
    fmax : float, optional
        The maximum frequency for fitting (default is 25000 Hz).

    Returns:
    alpha : float
        The optimal scaling factor for the time series.
    freqs : ndarray
        The frequencies corresponding to the power spectral density.
    """

    Nt = len(time_series)

    freqs       = np.fft.fftfreq(Nt, 1/sampling_rate)
    ss          = Knudsen_curve(seastate=SS, frequency=freqs)
    fit_mask    = (freqs >= fmin) & (freqs <= fmax)

    def loss_fn(alpha):
        ycorr_scaled = time_series / alpha
        Y_scaled     = np.fft.fft(ycorr_scaled)
        m_dat_scaled = np.abs(Y_scaled)**2
        p_dat_scaled = 10 * np.log10( m_dat_scaled / Nt / sampling_rate * 1e12)
  
        # Mean squared error between spectrum and Knudsen curve within fitting range
        return np.mean((p_dat_scaled[fit_mask] - ss[fit_mask]) ** 2)

    res     = minimize_scalar(loss_fn, bounds=(1e-6, 1e6), method='bounded')
    alpha   = res.x

    if res.success:
        return alpha, freqs
    else:
        raise RuntimeError("Optimization failed to find a suitable scaling factor.")
    
def time_domain_scaling_factor_from_SS_hanning(time_series, SS, sampling_rate=144e3, fmin=1000, fmax=25000, nperseg=1024):
    """
    Calculate the scaling factor for the time series based on the seastate.
    The time series should be zero-mean.
    The function uses the Knudsen curve to find the optimal scaling factor
    that minimizes the difference between the power spectral density of the
    time series and the Knudsen curve for the given seastate.

    Instead of a square window, the function uses a Hanning window to reduce spectral leakage.

    Parameters:
    time_series : array_like
        The time series data to be processed.
    SS : int
        The seastate level (0 to 6).
    sampling_rate : float, optional
        The sampling rate of the time series (default is 144e3 Hz).
    fmin : float, optional
        The minimum frequency for fitting (default is 1000 Hz).
    fmax : float, optional
        The maximum frequency for fitting (default is 25000 Hz).

    Returns:
    alpha : float
        The optimal scaling factor for the time series.
    freqs : ndarray
        The frequencies corresponding to the power spectral density.
    """

    Nt          = len(time_series)
    freqs,_     = signal.welch(time_series, sampling_rate, nperseg=nperseg, scaling= 'density')
    ss          = Knudsen_curve(seastate=SS, frequency=freqs) # Note this is in dB re 1 muPa^2/Hz
    fit_mask    = (freqs >= fmin) & (freqs <= fmax)

    def loss_fn(alpha):
        ycorr_scaled    = time_series / alpha
        _,Y_PSD         = signal.welch(ycorr_scaled, sampling_rate, nperseg=nperseg, scaling= 'density')
        p_dat_scaled    = 10 * np.log10(Y_PSD * 1e12)  

        # Mean squared error between spectrum and Knudsen curve within fitting range
        return np.mean((p_dat_scaled[fit_mask] - ss[fit_mask]) ** 2)

    res     = minimize_scalar(loss_fn, bounds=(1e-6, 1e6), method='bounded')
    alpha   = res.x

    if res.success:
        return alpha, freqs
    else:
        raise RuntimeError("Optimization failed to find a suitable scaling factor.")
    

########################### Scramble phase ###########################
def random_scrambled_trace_ss(time_series, SS, sample_frequency=144e3, method='Welch'):
    """
    Generate a phase-scrambled version of the input time series.
    """
    if method == 'Welch':
        alpha = correction_factor_Welch(SS)
    elif method == 'Square':
        alpha = correction_factor_Square(SS)
    else:
        raise ValueError("Unknown method specified.")

    Nt = len(time_series)
    time_series -= np.mean(time_series)  # Remove DC offset
    t0 = np.arange(Nt) / sample_frequency

    m_yFFT = np.abs(np.fft.fft(time_series))

    # Generate random phases
    phase = 2*pi* np.random.rand(Nt) - pi
    phase[0] = 0  # Set the phase of the DC component to zero to avoid an offset when taking the real part of the IFFT

    Z = m_yFFT * (np.cos(phase) + 1j*np.sin(phase))
    z = np.fft.ifft(Z)
    # mock_noise = (np.real(z) - np.mean(np.real(z)))
    mock_noise = np.real(z) / alpha

    return t0, mock_noise


def random_scrambled_trace_ss_Gauss(PSD, SS, df, sample_frequency=144e3, fmin=1000, fmax=72000, method='Welch'):
    """
    Generate a phase & amplitude-scrambled trace from PSD.
    PSD in dB re 1 Pa^2/Hz
    """
    if method == 'Welch':
        alpha = correction_factor_Welch(SS)
    elif method == 'Square':
        alpha = correction_factor_Square(SS)
    else:
        raise ValueError("Unknown method specified.")

    Nfs = len(PSD)
    Nt  = (Nfs - 1) * 2

    Ampl = np.sqrt(PSD * df / (2 * alpha * alpha)) * Nt # as the Welch method gives a 1/N factor, which we will have in the ifft as well
    kmin = int(fmin / df)
    kmax = int(fmax / df)

    t0 = np.arange(Nt) / sample_frequency

    # Generate random numbers
    FFTreal = np.random.normal(0, 1, Nfs) * Ampl
    FFTimag = np.random.normal(0, 1, Nfs) * Ampl

    Z = FFTreal + 1j * FFTimag
    Z[:kmin] = 0
    Z[kmax:] = 0
    z = np.fft.irfft(Z)

    # mock_noise = (np.real(z) - np.mean(np.real(z)))
    mock_noise = np.real(z)#  / alpha

    return t0, mock_noise


def random_scrambled_trace_ss_zerophase(time_series, SS, sample_frequency=144e3, method='Welch'):
    """
    Generate a 0-phase version of the input time series.
    """
    if method == 'Welch':
        alpha = correction_factor_Welch(SS)
    elif method == 'Square':
        alpha = correction_factor_Square(SS)
    else:
        raise ValueError("Unknown method specified.")

    Nt = len(time_series)
    time_series -= np.mean(time_series)  # Remove DC offset
    t0 = np.arange(Nt) / sample_frequency

    m_yFFT = np.abs(np.fft.fft(time_series))

    # Generate random phases
    phase = np.zeros(Nt)

    Z = m_yFFT * (np.cos(phase) + 1j*np.sin(phase))
    z = np.fft.ifft(Z)
    # mock_noise = (np.real(z) - np.mean(np.real(z)))
    mock_noise = np.real(z) / alpha

    return t0, mock_noise