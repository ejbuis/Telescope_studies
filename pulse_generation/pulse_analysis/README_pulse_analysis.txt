The scripts in this folder use the kernelfrmod function from the main pulse_generation
folder to generate ther rztp maps both in .mat and binary format to be used
inside the JPP full telescope simulations. The main scrips are:

1) rztp_tablefunc.m
    based on a neutrino energy Eo in GeV and a grid specification (to be
    obtained from auxiliary/generate_coordinates.m) it creates a rztp map
    in binary format (and in .mat if desired)

2) resample_rztmap.m
    based on a .mat rztp inputfile, a original sampling frequency and a 
    target sampling frequency this function generates a resampled map and 
    stores it in a .bin format

3) template_bankgen.m
    This script creates a template bank to be used as an input for the
    matched filtering used in JPP

4) generate_resampled_map.m
    This is an example script for how to use the resampling function