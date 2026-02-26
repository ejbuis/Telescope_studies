import numpy as np
from fibre_hydrophone import *

## ============================================================================
## Python implementation of the matched filter algorithm as designed for
## integration in JPP, where the computation of random phases for the noise
## realization is integrated in the frequency domain Matched filter computation
##
## ============================================================================

def matched_filter_JPP(p_test, padded_template, psd, sampling_rate, flow, fhigh, alpha=None):
  '''
  Python version of the JPP matched filtering implementation where the computation of
  random phases for a noise realization is included in the frequency-domain matched filter
  computation.

  Inputs:

  '''
  # Check on input size
  assert len(p_test) == len(padded_template), "data and template must match in length"
  assert (len(psd) - 1) * 2 == len(p_test), "PSD must be size len(p_test) / 2 + 1"
  L_test        = len(p_test)

  # FFTs
  sig_FFT       = np.fft.fft(p_test)
  templ_FFT     = np.fft.fft(padded_template)

  # Filter coefficients
  if alpha:
  
    indices = np.arange(L_test)

    Homega  = kspace_coefficient(alpha, indices, L_test, sampling_rate)
  
  else:
    Homega = np.ones(L_test)

  # Frequency grid -> as done in PyCBC
  dff   = sampling_rate / L_test
  kmin  = int(flow / dff)
  kmax  = int(fhigh / dff)

  # limit kmin and kmax to valid range
  if kmin <= 0:
    kmin = 1
  if kmax > int((L_test + 1)/2.):
    kmax = int((L_test + 1)/2.)

  mask            = np.zeros(L_test, dtype=bool)
  mask[kmin:kmax] = True

  psd_full                = np.zeros(L_test)
  psd_full[:len(psd)]     = psd
  psd_full[-len(psd)+1:]  = psd[1:][::-1] # mirror for negative frequencies
  
  # Core matched filter
  tmp           = sig_FFT * np.conj(templ_FFT) / psd_full
  tmp[~mask]    = 0
  Homega[~mask] = 0

  q             = 2 * np.fft.ifft(tmp * Homega) # factor 2 to account for taking only positive spectrum

  # Normalization (must apply same mask!)
  hnorm         = 4 * dff * np.sum((np.abs(templ_FFT[mask])**2) / psd_full[mask])
  norm          = 2 / np.sqrt(hnorm)

  # Final SNR time series
  snr           = q * norm

  return tmp, snr
