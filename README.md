RTSOS — Radiative Transfer model based on Successive Orders of Scattering
Copyright © 2025 Pengwang Zhai, University of Maryland Baltimore County
This software is released under the Creative Commons Attribution–NonCommercial 4.0 International (CC BY-NC 4.0) license.
You are free to use, modify, and share it for research and educational purposes, provided that you give appropriate credit.
Commercial use is not permitted without prior written consent from the author.
📧 For commercial licensing inquiries, contact Pengwang Zhai at [pwzhai@umbc.edu].
🔗 Full license text: Creative Commons BY-NC 4.0

RTSOS can solve multiple scattering radiative transfer equation from UV, visible, to infrared. It can handle atmosphere-land or atmosphere-ocean coupled systems. The atmosphere can be a mixture of molecules, aerosols, and cloud droplets. The land bottom can be Lambertian, snow surface, Ross-Li, and a number of other surfaces. The ocean waters are modeled by a mixture of pure ocean water, phytoplankton, and colored dissolved, organic matter, and other hydrsols. The sensors can be placed at arbitrary levels in the Earth system. The output of the sensor can include the full polarized Stokes parameters (I, Q, U, V). For more information, see References at the end of this document as well as the documentations located at 
/Documentation/Describtion_of_input_files




PACE simulator RSR files:

PACE_OCI_L1B_LUT_RSR_baseline_1.1.1.nc
HARP2_RSR_Composite_v01b.csv
spexone_isrf.nc

which can be downloaded from:

https://oceancolor.gsfc.nasa.gov/data/pace/characterization/

Section 1: COMPILE
The steps to use it is: compile the package, configure the input files, and run the code with the executable file.
To compile, open a terminal and go to:
$cd RTSOS/SOS_Callable/compile
then try to compile the monochromatic code:
$make
If this is good, then edit the PACE Simulator make file to set the hdf5 library paths:
$emacs makefile_PACE_Simulator_DoubleK

HDF5DIR = /usr/local/Cellar/hdf5@1.10/1.10.7
HDF5LIB=-I$(HDF5DIR)/include
H5FC=$(HDF5DIR)/bin/h5fc
LIBSHDF= $(HDF5LIB) -L$(HDF5DIR)/lib/ -lhdf5 -lhdf5_fortran

If these are correctly set, you can compile the pace simulator:
$make -f makefile_PACE_Simulator_DoubleK
This will compile the pace simulator.

Section 2: RUN the MODELS

Section 2.1 PACE Simulator
Next we may try out the pace simulator:
$cp rtsos_PACE_Simulator_DoubleK.exe ../test/PACE_Simulator/
$cd ../test/PACE_Simulator/

Check out “notes_to_run_pace_simulator” for the use of
pace simulator, which is located in
RTSOS/SOS_Callable/test/PACE_Simulator/
In the same folder you may also find a folder “python_script”
which contains three python scripts including “pace_simulator_inputfile_prepare.py”, 
which can be used to generate a batch of pace simulator input files.

Current version can conveniently generate synthetic datasets for OCI, HARP, and SPEX
for flexible datasets for OCI, HARP, and SPEX for flexible atmospheric and ocean condition. 
The input parameters are specified in input files.
Right now the aerosol models are: IAEROSOL=1-20, where
IAEROSOL=1-10 is the Shettle&Fenn model:
E. P. Shettle and R. W. Fenn, “Models for the aerosols of the lower atmosphere and the effects of humidity variations on their optical properties,” AFGL-TR 790214, U. S. Air Force Laboratory, Hanscom Air Force Base, Mass. (1979).

and IAEROSOL=11-20 is Zia’s new aerosol model. 
IRH=1-8 where RH(IRH)=(/0.30,0.50,0.70,0.75,0.80,0.85,0.90,0.95/).
Ziauddin Ahmad, Bryan A. Franz, Charles R. McClain, Ewa J. Kwiatkowska, Jeremy Werdell, Eric P. Shettle, and Brent N. Holben, "New aerosol models for the retrieval of aerosol optical thickness and normalized water-leaving radiances from the SeaWiFS and MODIS sensors over coastal regions and open oceans," Appl. Opt. 49, 5545-5560 (2010)



References:
Zhai, P., Hu, Y., Trepte, C. R., Lucker, P. L. (2009). A vector radiative transfer model for coupled atmosphere and ocean systems based on successive order of scattering method. Optics express, 17(4), 2057-2079.
Zhai, P., Hu, Y., Chowdhary, J., Trepte, C. R., Lucker, P. L., Josset, D. B. (2010). A vector radiative transfer model for coupled atmosphere and ocean systems with a rough interface. Journal of Quantitative Spectroscopy and Radiative Transfer, 111(7), 1025-1040.
Zhai, P., Hu, Y., Josset, D. B., Trepte, C. R., Lucker, P. L., Lin, B. (2013). Advanced angular interpolation in the vector radiative transfer for coupled atmosphere and ocean systems. Journal of Quantitative Spectroscopy and Radiative Transfer, 115, 19-27.
Zhai, P., Hu, Y., Winker, D. M., Franz, B., Boss, E. (2015). Contribution of Raman scattering to polarized radiation field in ocean waters. Optics Express, 23(18), 23582-23596.
Zhai, P., Hu, Y., Winker, D. M., Franz, B., WERDELL, J., BOSS, E. (2017). Inelastic vector radiative transfer solution in ocean waters. Optics Express, 25, A223-A239. 
Zhai, P., Boss, E., Franz, B., Werdell, J., Hu, Y. (2018). Radiative transfer modeling of phytoplankton fluorescence quenching  processes. MDPI Remote Sensing, 10(8), 1039.
Zhai, P., Hu, Y. (2022). An improved pseudo spherical shell algorithm for vector radiative transfer. Journal of Quantitative Spectroscopy and Radiative Transfer, 282, 108132. https://www.sciencedirect.com/science/article/pii/S0022407322000693.
Zhai, P., Gao, M., Franz, B. A., Werdell, P. J., Ibrahim, A., Hu, Y., Chowdhary, J. (2022). A Radiative Transfer Simulator for PACE: Theory and Applications. Frontiers in Remote Sensing, 3. https://www.frontiersin.org/article/10.3389/frsen.2022.840188.
