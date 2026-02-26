The Modified sources folder is used to analyse the effect of the source 
geometry, parameterisation and binning on the acoustic pulse.
The Showerparam_extensions foler includes plot and scripts to generate 
spherical, disk and cylinder sources

The main scripts in this folder are:

1) compound_fit_shower.m
    This script fits a combination of the amplitude decay of a spherical
    source with the decay of a cylindrical source to a typical 1e9 GeV shower

2) pulse_multimple_nmc.m
    This scripts evaluates the impact of the number of mc points to the pulse
    shape at an individual pulse level

3) shower_vs_sphere_and_cylinder.m
    A more elaborate version of compound_fit_shower.m, also showing the exact
    profile obtained through the cylindrical and spherical sources from the
    Showerparam_extensiosn folder

4) zmax_analysis.m
    Obtain the energy deposition maximum for the CORSIKA hadronic shower
    parameterizations as included by the ACoRNE contribution
