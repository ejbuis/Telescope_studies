from directionality_core import *

# =============================================================
# Main run-file to generate a plot of F_theta for the Pylos
# operating conditions and/or a set of provided f and z values
# =============================================================

def main(argv):
    refl_bool = 0
    theta_deg = np.linspace(0, 90, 100)
    z_m       = np.linspace(1000, 4000, 100)
    f_Hz      = np.geomspace(1e3, 25e3, 100)
    
    plot_extremes_pylos(refl_bool)

    if len(argv) >= 3:
        f_fixed = float(argv[1])
        z_fixed = float(argv[2]) 
        print(f"Running with f = {f_fixed:.2e} Hz and z = {z_fixed:.2f} m")
        plot_Ftheta_vs_theta_1D(f_fixed=f_fixed, z_fixed=z_fixed, theta_deg=theta_deg, refl=refl_bool)

    if len(argv) == 4:
        plot_Ftheta_vs_f_theta(z=z_fixed, f_Hz=f_Hz, theta_deg=theta_deg)
        plot_Ftheta_vs_theta_z(f_fixed=f_fixed, z_m=z_m, theta_deg=theta_deg)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
