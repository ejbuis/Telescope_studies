import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scipy import interpolate, optimize
from mpmath import gammainc
import sys
sys.path.append('..')
import report_plotstyle

## =====================================================================
## Python package that contains all the main scripts for computing the 
## noise directionality density F(theta) and the various approximations 
## of an effective noise angle.
##
## The definition of F(theta) and the normalization condition have been
## taken from:
## https://pubs.aip.org/asa/jasa/article/134/2/950/950939
## =====================================================================

# ======================================================================
# Directionality fits to the hydrophone properties of:
# https://www.sciencedirect.com/science/article/pii/S0927650525000325
# ======================================================================

HYDROANGLESDEG          = np.array([0,30,60,90,120,150,180])
HYDRORESONANCEdB        = np.array([9.9381, 14.659, 26.621, 28.543, 21.574, 12.330, 8.7958])
HYDROCONSTdB            = np.array([0.8447, 2.626, 5.755, 5.486, 3.643, 1.859, 0.6858])
HYDRORESONANCEdBPCHIP   = PchipInterpolator(HYDROANGLESDEG, HYDRORESONANCEdB)
HYDROCONSTdBPCHIP       = PchipInterpolator(HYDROANGLESDEG, HYDROCONSTdB)

# ======================================================================
# Compute the various integrands of the normalization condition
# ======================================================================

def Ftheta(theta_deg, z_m, f_Hz):
    """
    F(theta) for no attenuation, defined on theta_deg in [0, 90]

        F(theta) = 4 cos(theta)
    
    """
    if theta_deg > 90:
        return 0.
    else:
        theta_rad = theta_deg * np.pi / 180
        return 4 * np.cos(theta_rad)

def alpha_approx(f_Hz):
    """
    Approximate linear fit to the attenuation curve used by: 
        Michael J. Buckingham. Theory of the directionality and spatial coherence of wind-driven
        ambient noise in a deep ocean with attenuation

    """
    alpha_min   = 1.5e-7
    alpha_max   = 1e-3
    fmin        = 1e-1
    fmax        = 1e2
    slope       = (alpha_max - alpha_min) / (fmax - fmin)
    
    f_kHz = f_Hz / 1e3
    if (f_kHz < fmin or f_kHz > fmax):
        assert("frequency outside range")
    
    return alpha_min + slope * f_kHz

def Ftheta_attenuation(theta_deg, z_m, f_Hz):
    """
    F(theta) considering an attenuation of alpha(z, f), defined on theta_deg in [0,90]

        F(theta) = [2 alpha z^2 Gamma(-2, 2 alpha z)]^(-1) cos(theta) exp(-2 alpha z / cos(theta) )

    From Buckingham
    
    """
    if theta_deg > 90:
        return 0.
    else:
        theta_rad   = theta_deg * np.pi / 180
        alpha       = alpha_approx(f_Hz)
        beta        = 2 * alpha * z_m
        gammaval    = float(gammainc(-2, beta))
        B           = 1 / (2 * (alpha**2) * (z_m**2) * gammaval)
        return B * np.cos(theta_rad) * np.exp(-beta / np.cos(theta_rad))

# =====================================================
# Helper functions
# =====================================================

def reflect(theta_deg, Ftheta):
    """
    Helper function to reflect all rays perfectly in a mirroring sea floor

    """

    theta_ext   = (theta_deg + theta_deg[-1])[1:]
    theta_refl  = np.append(theta_deg[:], theta_ext)
    Ftheta_ext  = Ftheta[::-1]
    Ftheta_refl = np.append(Ftheta[:]/2, Ftheta_ext[1:]/2)
    return theta_refl, Ftheta_refl

def hydrophone_resonance(theta_deg):
    """
    Helper function that returns the height (in linear units) of the hydrophone resonance peak at a given angle
    Angle defined on [0, 180]

    """
    G_res_dB        = np.interp(theta_deg, HYDROANGLESDEG, HYDRORESONANCEdB)
    resonance_lin   = 10 ** (G_res_dB / 20)

    return resonance_lin


def hydrophone_response(theta_deg, f_Hz):
    """
    Helper function that returns the height (in linear units) of the hydrophone response for any combination of theta and f
    Angle defined on    [0, 180     ] deg
    f defined on        [0, 72000   ] Hz

    """
    sampling_rate   = 144e3
    f_peak          = 16200
    f_width         = 750

    G_res_dB        = HYDRORESONANCEdBPCHIP(theta_deg) 
    G_HYDROCONSTdB  = HYDROCONSTdBPCHIP(theta_deg)

    g_lin  = 10.0 ** (G_HYDROCONSTdB / 20.0)
    g_peak = 10.0 ** (G_res_dB / 20.0)
    g_r    = g_peak - g_lin

    R      = np.exp(-np.pi * f_width / sampling_rate)
    theta  = 2 * np.pi * f_peak / sampling_rate
    g_corr = g_r * (1 - R * R) / 2

    a = np.zeros((len(theta_deg), 3))
    b = np.zeros((len(theta_deg), 3))

    for i in range(len(theta_deg)):
        a[i,0] = 1 
        a[i,1] = -2 * R * np.cos(theta)
        a[i,2] = R * R

        b[i,0] = g_corr[i] + a[i,0] * g_lin[i]
        b[i,1] = a[i,1] * g_lin[i]
        b[i,2] = -g_corr[i] + a[i,2] * g_lin[i]

    # Compute k space coefficient
    omega = 2 * np.pi * f_Hz / sampling_rate
    ejw = np.exp(-1j * omega)
    num = b[:,0] + b[:,1] * ejw + b[:,2] * (ejw * ejw)
    den = a[:,0] + a[:,1] * ejw + a[:,2] * (ejw * ejw)

    return np.abs(num / den)

# =====================================================
# Evaluation metric functions
# =====================================================

def compute_half_width(theta_deg, Ftheta):
    """
    Compute The value of theta_deg for which

        F(theta_HW) = F(0) / 2

    return

        theta_HW
        F(0)

    """
    peak = np.max(Ftheta)
    HW_val = peak / 2

    if np.argmax(Ftheta) != 0:
        print("NOTICE: the peak value not found at i =0 \n")

    theta_HW = np.interp(HW_val, Ftheta[::-1], theta_deg[::-1])

    return theta_HW, HW_val


def compute_half_width_hydrophone(theta_deg, Ftheta, f_Hz):
    """
    Compute The value of theta_deg for which

        F(theta_HW)H(theta_HW, f) = F(0)H(0, f) / 2

    return

        theta_HW
        F(0)
        theta_max

    """
    FthetaHtheta = Ftheta * hydrophone_response(theta_deg, f_Hz)

    peak        = np.max(FthetaHtheta)
    HW_val      = peak / 2

    theta_max   = theta_deg[np.argmax(FthetaHtheta)]

    theta_HW    = np.interp(HW_val, FthetaHtheta[::-1], theta_deg[::-1])

    return theta_HW, HW_val, theta_max


def compute_3dB(theta_deg, Ftheta):
    """
    Compute The value of theta_deg for which

        F(theta_3dB) = F(0) / sqrt(2)

    return

        theta_3dB
        F(0)

    """
    peak = np.max(Ftheta)
    thres_val = peak / np.sqrt(2)

    if np.argmax(Ftheta) != 0:
        print("NOTICE: the peak value not found at i =0 \n")

    theta_3dB = np.interp(thres_val, Ftheta[::-1], theta_deg[::-1])
    
    return theta_3dB, thres_val


def compute_F1(theta_deg, Ftheta):
    """
    Compute The value of theta_deg for which

        F(theta_F1) = 1

    return

        theta_F1
        1

    """    
    HW_val = 1

    if np.argmax(Ftheta) != 0:
        print("NOTICE: the peak value not found at i =0 \n")

    theta_F1 = np.interp(HW_val, Ftheta[::-1], theta_deg[::-1])

    return theta_F1, HW_val


def compute_halfint(theta_deg, Ftheta):
    """
    Compute The value of theta_deg for which

        C = int_0^90 ( F(theta) sin(theta) )

        int_0^theta_halfint ( F(theta) sin(theta) ) = C / 2

    return

        theta_halfint
        F(theta_halfint)
        
    """
    theta_rad = np.radians(theta_deg)
    integrand = Ftheta * np.sin(theta_rad)

    # Cumulative integral using the trapezoidal rule
    cumulative = cumulative_trapezoid(integrand, theta_rad, initial=0)
    tot_int    = cumulative[-1]
    
    idx = np.searchsorted(cumulative, tot_int/2)

    # Linear interpolation for better accuracy
    if idx == 0 or idx >= len(theta_deg):
        theta_halfint   = theta_deg[idx]
        Ftheta_halfint  = Ftheta[idx]
    else:
        x0, x1          = theta_deg[idx-1], theta_deg[idx]
        y0, y1          = cumulative[idx-1], cumulative[idx]
        theta_halfint   = x0 + (0.5 - y0) * (x1 - x0) / (y1 - y0)

        # Interpolate Ftheta at that angle
        Ftheta_halfint  = np.interp(theta_halfint, theta_deg, Ftheta)

    return theta_halfint, Ftheta_halfint

def compute_halfint_hydrophone_res(theta_deg, Ftheta):
    """
    Compute The value of theta_deg for which

        C_res = int_0^90 ( F(theta) H(theta, f_res) sin(theta) )

        int_0^theta_halfint ( F(theta) H(theta, f_res) sin(theta) ) = C_res / 2

    return

        theta_halfint
        F(theta_halfint)
        
    """
    theta_rad = np.radians(theta_deg)
    integrand = Ftheta * np.sin(theta_rad) * hydrophone_resonance(theta_deg)

    # Cumulative integral using the trapezoidal rule
    cumulative = cumulative_trapezoid(integrand, theta_rad, initial=0)
    tot_int    = cumulative[-1]
    
    idx = np.searchsorted(cumulative, tot_int/2)

    # Linear interpolation for better accuracy
    if idx == 0 or idx >= len(theta_deg):
        theta_halfint = theta_deg[idx]
        Ftheta_halfint = Ftheta[idx]
    else:
        x0, x1 = theta_deg[idx-1], theta_deg[idx]
        y0, y1 = cumulative[idx-1], cumulative[idx]
        theta_halfint = x0 + (0.5 - y0) * (x1 - x0) / (y1 - y0)

        # Interpolate Ftheta at that angle
        Ftheta_halfint = np.interp(theta_halfint, theta_deg, Ftheta)

    return theta_halfint, Ftheta_halfint


def compute_const_hydrophone_res(theta_deg, Ftheta):
    """
    Compute The value of theta_deg for which

        C_res = int_0^90 ( F(theta) H(theta, f_res) sin(theta) )

        F(theta_const) H(theta_const, f_res) = C_res / 2

    return

        theta_const
        C_res
        
    """
    theta_rad       = np.radians(theta_deg)
    FthetaHtheta    = np.sin(theta_rad) * hydrophone_resonance(theta_deg)
    integrand       = Ftheta * FthetaHtheta

    # Cumulative integral using the trapezoidal rule
    cumulative = cumulative_trapezoid(integrand, theta_rad, initial=0)
    tot_int    = cumulative[-1]

    # The total power is given by integral of F(theta)H(theta)sin(theta)
    # The same total power is obtained for a constan value of theta_const, given by
    f_const = tot_int / 2

    # Interpolant for root finding (difference = 0 where F(theta)*H(theta) == f_const)
    interp = interpolate.interp1d(theta_rad, FthetaHtheta - f_const, kind='cubic')

    roots = []
    for i in range(len(theta_rad) - 1):
        if np.sign(FthetaHtheta[i] - f_const) != np.sign(FthetaHtheta[i+1] - f_const):
            root = optimize.brentq(interp, theta_rad[i], theta_rad[i+1])
            roots.append(root)

    # Convert to degrees for readability
    theta_const = np.degrees(roots)

    return theta_const, tot_int


def compute_const_hydrophone_full(theta_deg, Ftheta, f_Hz):
    """
    Compute The value of theta_deg for which

        C_res = int_0^90 ( F(theta) H(theta, f) sin(theta) )

        F(theta_halfint_f) H(theta_halfint_f, f) = C_res / 2

    return

        theta_const_f
        C_res
        
    """

    theta_rad       = np.radians(theta_deg)
    FthetaHtheta    = np.sin(theta_rad) * hydrophone_response(theta_deg, f_Hz)
    integrand       = Ftheta * FthetaHtheta

    # Cumulative integral using the trapezoidal rule
    cumulative = cumulative_trapezoid(integrand, theta_rad, initial=0)
    tot_int    = cumulative[-1]

    # The total power is given by integral of F(theta)H(theta)sin(theta)
    # The same total power is obtained for a constan value of theta_const, given by
    f_const = tot_int / 2 

    # Interpolant for root finding (difference = 0 where F(theta)*H(theta) == f_const)
    interp = interpolate.interp1d(theta_rad, FthetaHtheta - f_const, kind='cubic')

    roots = []
    for i in range(len(theta_rad) - 1):
        if np.sign(FthetaHtheta[i] - f_const) != np.sign(FthetaHtheta[i+1] - f_const):
            root = optimize.brentq(interp, theta_rad[i], theta_rad[i+1])
            roots.append(root)

    # Convert to degrees for readability
    theta_const_f = np.degrees(roots)

    return theta_const_f, tot_int


# =====================================================
# Plotting functions
# =====================================================

def plot_Ftheta_vs_f_theta(z, f_Hz, theta_deg):
    """
    2D plot: x=theta, y=f_Hz, for fixed depth z.
    Adds isoline showing the half-width for each frequency.
    """
    F, THETA = np.meshgrid(f_Hz, theta_deg, indexing='ij')

    F_val = np.vectorize(lambda f, t: Ftheta(t, z, f))(F, THETA)
    F_att_val = np.vectorize(lambda f, t: Ftheta_attenuation(t, z, f))(F, THETA)   

    # --- Compute local HWs (for each frequency)
    theta_HW_F = []
    theta_HW_F_att = []
    for i in range(len(f_Hz)):
        theta_HW_F.append(np.interp(np.max(F_val[i,:])/2, F_val[i,::-1], theta_deg[::-1]))
        theta_HW_F_att.append(np.interp(np.max(F_att_val[i,:])/2, F_att_val[i,::-1], theta_deg[::-1]))
    theta_HW_F = np.array(theta_HW_F)
    theta_HW_F_att = np.array(theta_HW_F_att)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

    # --- Ftheta
    im1 = axes[0].pcolormesh(theta_deg, f_Hz, F_val, shading='auto')
    axes[0].plot(theta_HW_F, f_Hz, 'r-', linewidth=2, label='Half width')  # NEW
    axes[0].set_title(f"Ftheta (z = {z:.0f} m)")
    axes[0].set_xlabel("θ [deg]")
    axes[0].set_ylabel("f [Hz]")
    fig.colorbar(im1, ax=axes[0], label="Ftheta")
    axes[0].legend()

    # --- Ftheta_attenuation
    im2 = axes[1].pcolormesh(theta_deg, f_Hz, F_att_val, shading='auto')
    axes[1].plot(theta_HW_F_att, f_Hz, 'r-', linewidth=2, label='Half width')  # NEW
    axes[1].set_title(f"Ftheta_attenuation (z = {z:.0f} m)")
    axes[1].set_xlabel("θ [deg]")
    axes[1].set_ylabel("f [Hz]")
    fig.colorbar(im2, ax=axes[1], label="Ftheta_attenuation")
    axes[1].legend()

    plt.show()


def plot_Ftheta_vs_theta_z(f_fixed, z_m, theta_deg):
    """
    2D plot: x=theta, y=z_m, for fixed frequency f_fixed.
    Adds isoline showing the half-width for each depth.
    """
    Z, THETA = np.meshgrid(z_m, theta_deg, indexing='ij')

    F_val = np.vectorize(lambda z, t: Ftheta(t, z, f_fixed))(Z, THETA)
    F_att_val = np.vectorize(lambda z, t: Ftheta_attenuation(t, z, f_fixed))(Z, THETA)

    # --- Compute local HWs (for each depth)
    theta_HW_F = []
    theta_HW_F_att = []
    for i in range(len(z_m)):
        theta_HW_F.append(np.interp(np.max(F_val[i,:])/2, F_val[i,::-1], theta_deg[::-1]))
        theta_HW_F_att.append(np.interp(np.max(F_att_val[i,:])/2, F_att_val[i,::-1], theta_deg[::-1]))
    theta_HW_F = np.array(theta_HW_F)
    theta_HW_F_att = np.array(theta_HW_F_att)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

    im1 = axes[0].pcolormesh(theta_deg, z_m, F_val, shading='auto')
    axes[0].plot(theta_HW_F, z_m, 'r-', linewidth=2, label='Half width')  # NEW
    axes[0].set_title(f"Ftheta (f = {f_fixed:.0f} Hz)")
    axes[0].set_xlabel("θ [deg]")
    axes[0].set_ylabel("z [m]")
    fig.colorbar(im1, ax=axes[0], label="Ftheta")
    axes[0].legend()

    im2 = axes[1].pcolormesh(theta_deg, z_m, F_att_val, shading='auto')
    axes[1].plot(theta_HW_F_att, z_m, 'r-', linewidth=2, label='Half width')  # NEW
    axes[1].set_title(f"Ftheta_attenuation (f = {f_fixed:.0f} Hz)")
    axes[1].set_xlabel("θ [deg]")
    axes[1].set_ylabel("z [m]")
    fig.colorbar(im2, ax=axes[1], label="Ftheta_attenuation")
    axes[1].legend()

    plt.show()


def plot_Ftheta_vs_theta_1D(f_fixed, z_fixed, theta_deg, refl=0):
    """
    1D plot: Ftheta and Ftheta_attenuation vs theta.
    """
    F_vals      = np.array([Ftheta(theta, z_fixed, f_fixed) for theta in theta_deg])
    F_att_vals  = np.array([Ftheta_attenuation(theta, z_fixed, f_fixed) for theta in theta_deg])
    
    theta_HW_F, HW_F            = compute_half_width(theta_deg, F_vals)
    theta_HW_F_att, HW_F_att    = compute_half_width(theta_deg, F_att_vals)

    theta_rad = np.deg2rad(theta_deg)
    print(np.trapezoid(F_vals*np.sin(theta_rad), theta_rad), np.trapezoid(F_att_vals*np.sin(theta_rad), theta_rad))

    if refl:
        [theta_deg_plot, F_vals] = reflect(theta_deg, F_vals)
        [_, F_att_vals]           = reflect(theta_deg, F_att_vals)
    else:
        theta_deg_plot = theta_deg    

    plt.figure(figsize=(8, 5))
    plt.plot(theta_deg_plot, F_vals, 'r', label="Ftheta", linewidth=2)
    plt.plot(theta_HW_F, HW_F, 'xr', markersize=8, markeredgewidth=2)
    plt.plot(theta_deg_plot, F_att_vals, 'b--', label="Ftheta_attenuation", linewidth=2)
    plt.plot(theta_HW_F_att, HW_F_att, 'xb', markersize=8, markeredgewidth=2)

    plt.xlabel("θ [deg]")
    plt.ylabel("Amplitude")
    plt.title(f"Fθ and Fθ_attenuation vs θ\n(f = {f_fixed:.0f} Hz, z = {z_fixed:.0f} m)")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


def plot_extremes_pylos(refl=0):
    """
    1D plot: extreme cases with HW crosses.
    """

    theta_deg = np.linspace(0, 90, 100)
    
    F_vals_undeep   = np.array([Ftheta_attenuation(theta, 1000, 1000) for theta in theta_deg])
    F_vals_deep     = np.array([Ftheta_attenuation(theta, 4000, 25000) for theta in theta_deg])

    theta_HW_undeep, HW_undeep = compute_half_width(theta_deg, F_vals_undeep)
    theta_HW_deep, HW_deep     = compute_half_width(theta_deg, F_vals_deep)

    if refl:
        [theta_deg_plot, F_vals_undeep] = reflect(theta_deg, F_vals_undeep)
        [_, F_vals_deep]                = reflect(theta_deg, F_vals_deep)
    else:
        theta_deg_plot = theta_deg    

    plt.figure(figsize=(6.4, 4.8))
    plt.plot(theta_deg_plot, F_vals_undeep, label="z = 1 km, f = 1 kHz", linewidth=2)
    plt.plot(theta_deg_plot, F_vals_deep, label="z = 4 km, f = 25 kHz", linewidth=2)

    # Undeep
    plt.plot([theta_HW_undeep, theta_HW_undeep], [0, HW_undeep],
            '--', color='grey')
    plt.plot([0, theta_HW_undeep], [HW_undeep, HW_undeep],
            '--', color='grey')

    # Deep
    plt.plot([theta_HW_deep, theta_HW_deep], [0, HW_deep],
            '--', color='grey')
    plt.plot([0, theta_HW_deep], [HW_deep, HW_deep],
            '--', color='grey')
    
    plt.plot(theta_HW_undeep, HW_undeep, 'x', color='k', markersize=8, markeredgewidth=2)
    plt.plot(theta_HW_deep, HW_deep, 'x', color='k', markersize=8, markeredgewidth=2)

    plt.xlabel(r"$\theta$ [°]")
    plt.ylabel(r"F($\theta$) [-]")
    plt.xlim(0, 90)
    plt.ylim(0,10)
    # plt.yscale('log') 
    # plt.ylim(1e-5, 20)
    # plt.title(f"Fθ vs θ at typical Pylos specs.")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()
