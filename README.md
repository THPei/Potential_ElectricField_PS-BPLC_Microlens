# Potential_ElectricField_PS-BPLC_Microlens
The software for the publication about The Distributions of The Electric Potential and Fields in The Polymer-Stabilized Blue Phase Liquid Crystals Between Two Parallel Electrodes
We present software for calculating the analytical distributions of potential and field in microlens made of polymer-stabilized blue-phase liquid crystals. The microlens is located between two electrodes with a circular opening on the upper electrode. The program is encoded using the derived analytical solution formulas of the potential and electric field. It is faster, more energy-saving, and has less memory configuration than the finite-difference equation method solving by using the Poisson Equation with boundary conditions. The software outputs relevant Data files and a MATLAB code that can draw the distributions of potential, electric field, and directors, allowing the calculations visualized.

In this code, please compile it in a C/C++ environment. For example. you can use the open C/C++ platform Dev C++ to compile it. Then according to the parameter input processes, please input all parameters for the programm running. The instructions for the parameter input processes are listed in the following:

Please enter the ordinary refractive index of Blue-phase Liquid Crystal when there is no external potential: no=

Please enter the extra-ordinary refractive index of Blue-phase Liquid Crystal when there is no external potential: ne= 

Please enter the dielectric constant when Blue-phase Liquid Crystal is driven by AC voltage: 

Please enter the potential difference applied to Blue-phase Liquid Crystal: 

Please enter the cell gap of Blue-phase Liquid Crystal (in units of micrometer): 

Please enter the Kerr constant of Blue-Phase Liquid Crystal (in units of nm/VV):

Please enter the saturation electric field of Blue-Phase Liquid Crystal (in units of V/um): 

Please enter the Threshold voltage that causes Blue-Phase Liquid Crystal to perform Kerr Effect (>=0): 

Please enter the sign to make the Kerr effect of Blue-Phase Liquid Crystal produce Birefringence (1 or -1): 

Please enter the number of additional mirror potentials required to satisfy the upper and lower plate boundary conditions: 

Please enter the opening radius of the Blue-phase Liquid Crystal lens (in units of micrometer): 

Please enter the wavelength of the incident light (in units of nanometer): 

Please enter the grid spacing in the x direction (radial direction) (in units of micrometer): 

Please enter the grid spacing in the z direction (radius direction) (in units of micrometer): 

After those input processes, the program will run automatically. After completing the running, the program will output a Matlab code to plot some figures and Voltage.dat, E_radius.dat, Ez.dat, E_radius_Dr.dat, E_radius_Dz.dat, Ez_Dr.dat, and Ez_Dz.dat for use. 
