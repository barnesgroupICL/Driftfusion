# Driftfusion
An open source drift diffusion code based in MATLAB for simulating solar cells

## Info
Authors: Phil Calado, Piers RF Barnes, Ilario Gelmetti, Mohammed Azzouzi, Ben Hilman
Imperial College London, 2019

If you use Driftfusion please let us know by emailing:
p.calado13@imperial.ac.uk

Please log bug fixes through GitHub.

## QuickStart Guide

To get started do the following:
 
1.	If your new to GitHub it is highly recommended that you download GitHub desktop at: https://desktop.github.com/
You can simply download Driftfusion standalone as a folder but you won’t easily be able to synch to the latest version.

2.	Fork the Driftfusion GitHub repository (project). Instructions can be found here: https://help.github.com/en/articles/fork-a-repo 

3.	Navigate to your repository/project folder in MATLAB.

4.	The parameters class PC defines all the parameters for your device. Create a parameters object by typing:
 
`par = pc;`
	
The default parameters simulate an inverted PEDOT/perovskite/PCBM stack. The  names of properties and their units are given in the commented code in the PC class.

5.	Get equilibrium solutions – type:
 
`soleq = equilibrate(par)`
 
Inside soleq are 4 solutions. For perovskites you should generally use soleq.i_sr (i = ion and sr = surface recombination).
 
6.	Try running a JV using:
 
`JV = doJV(soleq.i_sr, 1e-2, 100, 1, 1, 0, 1.4, 2)`
 
The various input arguments are discussed in the comments of each code. At any time you can plot the JV using plotJV(JVsol).
Try using doJV as a model for writing functions of your own.

Good luck!
