# Driftfusion
An open source drift diffusion code based in MATLAB for simulating solar cells

## Info
Authors: Philip Calado, Piers RF Barnes, Ilario Gelmetti, Mohammed Azzouzi, Benjamin Hilton

Imperial College London, 2019

If you use Driftfusion please let us know by emailing:
p.calado13@imperial.ac.uk

Please log bugs through GitHub.

## QuickStart Guide

To get started do the following:
 
1.	If you are new to GitHub it is highly recommended that you download GitHub desktop at: https://desktop.github.com/.
Alternatively, you can download Driftfusion standalone as a folder but you won’t easily be able to synch to the latest version.

2.	Fork the Driftfusion GitHub repository (project). Instructions can be found here: https://help.github.com/en/articles/fork-a-repo 

3.	Navigate to your repository/project folder in MATLAB.

4.	Type `initialise_df` into the command prompt. This adds the subfolders and associated functions to the file path.

5.	The parameters class `pc` defines all the parameters for your device. You can edit the parameters directly in `pc` (found in the `/Core folder` and create  a parameters object by typing:
 
`par = pc;`

The  names of properties and their units are given in the commented code in the `pc` class.
	
6. 	A better way of defining your device is to use the .csv files contained in the `/Input_files` folder. For example, to create a parameters object with default parameters to simulate an inverted PTPD/perovskite/PCBM stack type:

`par = pc('Input_files/ptpd_mapi_pcbm.csv');`

It is recommended to duplicate the existing .csv files to define your own device and use a programme like Excel or Open Office to edit them. Ensure that the material name given in the 'stack' column matches with one of the materials in the `Libraries/Index_of_Refraction_library.xls` if you wish to use the Beer-Lambert optical model.

7.	Obtain equilibrium solutions for your device by typing:
 
`soleq = equilibrate(par)`
 
Inside soleq are 2 solutions: `soleq.el` and `soleq.ion`. For perovskites you should generally use `soleq.ion`. `soleq.el` is the same device at equilibrium without mobile ions.
 
8.	Try running a JV using:
 
`JVsol = doJV(soleq.ion, 1e-2, 100, 1, 1, 0, 1.4, 3)`
 
The various input arguments are discussed in the comments of each code. At any time you can plot the JV using `dfplot.JV(JVsol,3)`.
Try using doJV as a model for writing functions of your own.

The example script shows you how to produce two devices with different transport layers and scan current-voltage scans at 50 mVs-1 for each and plot the current-voltage curve, an energy level diagram and the currents using `dfplot`.

Good luck!

## Replicating published paper's data

### Impedance Spectroscopy on homojunction model, applying an oscillating voltage profile

Submitted [on arXiv on may 2018](https://arxiv.org/abs/1805.06446), published on [Energy & Environmental Science on march 2019](https://pubs.rsc.org/en/content/articlelanding/2019/ee/c8ee02362j).

Check out the `2018-EIS` branch of this repository or [download directly just this branch](https://github.com/barnesgroupICL/Driftfusion/archive/2018-EIS.zip).

Then follow the instructions in the branch readme or on this [wiki page](https://github.com/barnesgroupICL/Driftfusion/wiki/2018-Impedance-Spectroscopy-on-homojunction-model).

