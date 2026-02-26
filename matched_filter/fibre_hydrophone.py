import numpy as np
from scipy import signal

## ============================================================================
## Python implementation of the costant-gain resonator model fit to the fibre
## hydrophone directionality response
## hydrophone data from:
##    https://www.sciencedirect.com/science/article/pii/S0927650525000325
##
## resonator model from:
##    https://ccrma.stanford.edu/~jos/filters/Constant_Peak_Gain_Resonator.html
##
## ============================================================================

# ============================================================================
# General properties
# ============================================================================  
HYDROFIT_SAMPLINGRATE = 144e3
HYDROFIT_ANGLES       = np.array([0, 30, 60, 90, 120, 150, 180])
HYDROFIT_GPEAK        = np.array([9.9381, 14.659, 26.621, 28.543, 21.574, 12.330, 8.7958])
HYDROFIT_GCONST       = np.array([0.8447, 2.626, 5.755, 5.486, 3.643, 1.859, 0.6858])
HYDROFIT_FPEAK        = np.array([16200, 16200, 16200, 16200, 16200, 16200, 16200])
HYDROFIT_FWIDTH       = np.array([750, 750, 750, 750, 750, 750, 750])

def hydrophonetf(angle, trace, sampling_rate = HYDROFIT_SAMPLINGRATE):
  '''
  Hydrophone constant peak-gain resonator transfer function for a given
  signal incidence angle to be applied to a full time trace

  Inputs:
    angle         [deg] - signal zenith angle of incidence
    trace         [Pa]  - 1D array of pressure values sampled at sampling_rate
    sampling_rate [Hz]  - sampling rate of the time trace

  Outputs:
    p_filt        [Pa]  - trace with the effect of the hydrophone response applied 
  
  Note that this function performs linear interpolation on the experimentally
  obtained values of the peak gain and constant gain.
  '''

  HYDROFIT_GPEAK_interp  = np.interp(angle, HYDROFIT_ANGLES, HYDROFIT_GPEAK)
  HYDROFIT_GCONST_interp = np.interp(angle, HYDROFIT_ANGLES, HYDROFIT_GCONST)
  HYDROFIT_FPEAK_interp  = np.interp(angle, HYDROFIT_ANGLES, HYDROFIT_FPEAK)
  HYDROFIT_FWIDTH_interp = np.interp(angle, HYDROFIT_ANGLES, HYDROFIT_FWIDTH)

  g_lin  = 10.0 ** (HYDROFIT_GCONST_interp / 20.0)
  g_peak = 10.0 ** (HYDROFIT_GPEAK_interp / 20.0)
  g_r    = g_peak - g_lin

  R      = np.exp(-np.pi * HYDROFIT_FWIDTH_interp / sampling_rate)
  theta  = 2 * np.pi * HYDROFIT_FPEAK_interp / sampling_rate
  g_corr = g_r * (1 - R * R) / 2

  a = np.array([1, -2 * R * np.cos(theta), R * R])
  b = np.array([g_corr + a[0] * g_lin, a[1] * g_lin, -g_corr + a[2] * g_lin])


  _, zi = signal.lfilter(b, a, trace, zi = [0,0])
  print(zi)
  p_filt,_ = signal.lfilter(b, a, trace, zi=zi)

  return p_filt


def kspace_coefficient(angle, k, NT, sampling_rate = HYDROFIT_SAMPLINGRATE):
  '''
  Hydrophone constant peak-gain resonator transfer function absolute value
  of k-space  coefficient for a given entry k of NT

  Inputs:
    angle         [deg] - signal zenith angle of incidence
    k             [-]   - k-space index
    NT            [-]   - total length of the time trace
    sampling_rate [Hz]  - sampling rate of the time trace

  Outputs:
    |H|           [-]   - absolute value of the hydrophone response
  
  Note that this function performs linear interpolation on the experimentally
  obtained values of the peak gain and constant gain.
  '''
  HYDROFIT_GPEAK_interp  = np.interp(angle, HYDROFIT_ANGLES, HYDROFIT_GPEAK)
  HYDROFIT_GCONST_interp = np.interp(angle, HYDROFIT_ANGLES, HYDROFIT_GCONST)
  HYDROFIT_FPEAK_interp  = np.interp(angle, HYDROFIT_ANGLES, HYDROFIT_FPEAK)
  HYDROFIT_FWIDTH_interp = np.interp(angle, HYDROFIT_ANGLES, HYDROFIT_FWIDTH)

  g_lin  = 10.0 ** (HYDROFIT_GCONST_interp / 20.0)
  g_peak = 10.0 ** (HYDROFIT_GPEAK_interp / 20.0)
  g_r    = g_peak - g_lin

  R      = np.exp(-np.pi * HYDROFIT_FWIDTH_interp / sampling_rate)
  theta  = 2 * np.pi * HYDROFIT_FPEAK_interp / sampling_rate
  g_corr = g_r * (1 - R * R) / 2

  a = np.array([1, -2 * R * np.cos(theta), R * R])
  b = np.array([g_corr + a[0] * g_lin, a[1] * g_lin, -g_corr + a[2] * g_lin])

  # Compute k space coefficient
  omega = 2 * np.pi * k / NT

  ejw = np.exp(-1j * omega)
  num = b[0] + b[1] * ejw + b[2] * (ejw * ejw)
  den = a[0] + a[1] * ejw + a[2] * (ejw * ejw)

  return np.abs(num / den)
    