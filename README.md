# Telescope_studies
This repository contains the software developed by Huib Baetsen to aid in the research of the performance of an acoustic neutrino telescope. Together with the full telescope simulation in JPP, these files cover the topics of:
- Acoustic pulse generation
- Acoustic noise study
- Matched filtering
- Full telescope data inspection
- Fiber hydrophone parameterization
As well as auxiliary topics such as interaction lengths.

This README will give an overview of each of the main folder, along with their status and the programming language(s) used.

## attenuation (python)
Single jupyter notebook to compare different (complex) acoustic attenuation models

## directionality (python)
Python code to enhance the directionality of ambient acoustic noise as theoretically modelled by Buckingham (https://pubs.aip.org/asa/jasa/article/134/2/950/950939)

## interaction length (python)
Single jupyter notebook to evaluate interaction lengths at high energies

## matched filter (python)
Python implementation of the JPP matched filter and performance comparisons with the PyCBC package. This analysis also includes the effect of the fibre hydrophone response

## noise processing (python)
Analysis and parameterization of the acoustic noise measured at the Pylos site

## pulse_generation (matlab)
The full ACoRNE-based matlab toolchain used to generate the pressure maps for the full acoustic telescope simulations. This is a quite elaborate folder. It contains:
- **Auxiliary** scripts, used throughout the folder
- **frequency_analysis**: scripts used to evaluate the effect of the frequency parameter in the generation of the pulses, as well as the impact of downsampling and the hydrophone response
- **pulse_analysis**: all the scripts used to generate the rztp-maps
- **source_analysis**: evaluation of different source geometries (sphere, cylinder), compared to the neutrino-induced shower

