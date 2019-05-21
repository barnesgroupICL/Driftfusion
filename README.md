# Driftfusion - Code for impedance spectroscopy simulation on mobile ions containing homojunction devices

An open source drift diffusion code based in MATLAB for simulating solar cells

Please see the [wiki](https://github.com/barnesgroupICL/Driftfusion/wiki) and the headers of the .m files for complete documentation.

## Info
Authors: Phil Calado, Piers RF Barnes, Ilario Gelmetti, Mohammed Azzouzi, Ben Hilman
Imperial College London, 2019

If you use Driftfusion please let us know by emailing:
p.calado13@imperial.ac.uk

Please log bug fixes through GitHub.

## QuickStart Guide

To get started do the following:
 
1.	If your new to GitHub it is highly recommended that you download GitHub desktop at: https://desktop.github.com/
You can simply download Driftfusion standalone as a folder but you won’t easily be able to sync to the latest version.

2.	Fork the Driftfusion GitHub repository (project). Instructions can be found here: https://help.github.com/en/articles/fork-a-repo 

3.	Navigate to your repository/project folder in MATLAB.

4.	The parameters class PC defines all the parameters for your device. Create a parameters object by typing:

``` matlab
par = pinParams;
```
	
The default parameters simulate an homojunction p-type/intrinsic perovskite/n-type stack. The names of properties and their units are given in the commented code in the pinParams.m file.

5.	Get equilibrium solutions – copy and paste:

``` matlab
[sol_eq, sol_eq_SR, sol_i_eq, sol_i_eq_SR, ssol_eq, ssol_eq_SR, ssol_i_eq, ssol_i_eq_SR, sol_1S, sol_1S_SR, sol_i_1S, sol_i_1S_SR, ssol_1S, ssol_1S_SR, ssol_i_1S, ssol_i_1S_SR] = equilibrate(par)
```

For perovskites you should generally use:

* `sol_i_eq_SR` for short circuit in dark;
* `sol_i_1S_SR` for short circuit at 1 sun illumination;
* `ssol_i_eq_SR` for a symmetrical solution at open circuit in dark;
* `ssol_i_1S_SR` for a symmetrical solution at open circuit at 1 sun illumination.

Where ssol = symmetrical solution, i = mobile ions, 1S = 1 sun illumination, SR = surface recombination.
 
The various input arguments are discussed in the comments of each code.

At any time you can plot the solution at the last time point using:

``` matlab
pinana(sol)
```

Good luck!

## Replicating published paper's data

### Impedance Spectroscopy on homojunction model, applying an oscillating voltage profile

Submitted [on arXiv on may 2018](https://arxiv.org/abs/1805.06446), published [on Energy & Environmental Science on march 2019](https://pubs.rsc.org/en/content/articlelanding/2019/ee/c8ee02362j).

Check out the `2018-EIS` branch of this repository or [download directly just this branch](https://github.com/barnesgroupICL/Driftfusion/archive/2018-EIS.zip).

``` bash
git checkout 2018-EIS
```

Open MATLAB and set this as the working directory.

In MATLAB, generate stabilised solutions with [equilibrate](https://github.com/barnesgroupICL/Driftfusion/blob/2018-EIS/equilibrate.m):

``` matlab
>> [sol_eq, sol_eq_SR, sol_i_eq, sol_i_eq_SR, ssol_eq, ssol_eq_SR, ssol_i_eq, ssol_i_eq_SR, sol_1S, sol_1S_SR, sol_i_1S, sol_i_1S_SR, ssol_1S, ssol_1S_SR, ssol_i_1S, ssol_i_1S_SR] = equilibrate(pinParams)
```

or, quicker with [equilibrate_minimal](https://github.com/barnesgroupICL/Driftfusion/blob/2018-EIS/equilibrate_minimal.m):

``` matlab
>> [sol_eq, sol_i_eq, sol_i_eq_SR, ssol_i_eq, ssol_i_eq_SR, sol_i_1S_SR, ssol_i_1S_SR] = equilibrate_minimal(pinParams)
```

For simulating a different set of parameters, specify the relative pinParams file as an argument to equilibrate:

``` matlab
>> [~, ~, sol_i_eq_SR, ~, ssol_i_eq_SR, sol_i_1S_SR, ssol_i_1S_SR] = equilibrate_minimal(pinParams_10kxSRH_001xmajority)
```

For plotting the free charges and ionic defect density profile and the energy levels you can use the `pinana` routine:

``` matlab
>> pinana(ssol_i_1S_SR)
```

Simulate impedance spectroscopy with oscillating voltage on a single solution with [ISwave_EA_single_exec](https://github.com/barnesgroupICL/Driftfusion/blob/2018-EIS/ISwave_EA_single_exec.m):

``` matlab
>> asymssol_i_1S_SR = asymmetricize(ssol_i_1S_SR)
>> asymssol_i_1S_SR = stabilize(asymssol_i_1S_SR)
>> ssol_i_1S_SR_is_100Hz_50mV = ISwave_EA_single_exec(asymssol_i_1S_SR, 5e-2, 1e2, 20, 40, true, false, 1e-8)
```

Analyse the obtained solution with [ISwave_single_analysis](https://github.com/barnesgroupICL/Driftfusion/blob/2018-EIS/ISwave_single_analysis.m):

``` matlab
>> ISwave_single_analysis(ssol_i_1S_SR_is_100Hz_50mV, false, true)
```

Generate a cell with solutions at various illuminations at open circuit conditions with [genIntStructs](https://github.com/barnesgroupICL/Driftfusion/blob/2018-EIS/genIntStructs.m):

``` matlab
>> [structs_oc, VOCs, ~] = genIntStructs(ssol_i_eq_1S, 1, 1e-3, 7, true)
```

or, more accurate but slower with [genIntStructsRealVoc](https://github.com/barnesgroupICL/Driftfusion/blob/2018-EIS/genIntStructsRealVoc.m):

``` matlab
>> [structs_oc_realvoc, VOCs_real] = genIntStructsRealVoc(ssol_i_eq_SR, 1, 1e-3, 7, true)
```

Generate a cell with solutions in dark at various applied voltages with [genVappStructs](https://github.com/barnesgroupICL/Driftfusion/blob/2018-EIS/genVappStructs.m):

``` matlab
>> structs_vapp = genVappStructs(sol_i_eq_SR, [0, 0.2, 0.4, 0.6, 0.8])
```

Simulate the impedance spectroscopy with oscillating voltage on a cell containing various structures (at different illuminations or applied voltages) with [ISwave_full_exec](https://github.com/barnesgroupICL/Driftfusion/blob/2018-EIS/ISwave_full_exec.m):

``` matlab
>> ISwave_oc = ISwave_full_exec(structs_oc, 1e9, 1e-2, 56, 2e-3, false, true, true)
```

For plotting the Bode plots of capacitance and impedance versus frequency of a previous EIS simulation:

``` matlab
>> IS_full_analysis_impedance(ISwave_oc)
```

For plotting the Bode plot of phase versus frequency of a previous EIS simulation:

``` matlab
>> ISwave_full_analysis_phase(ISwave_oc)
```

For plotting the Nyquist plot of a previous EIS simulation:

``` matlab
>> ISwave_full_analysis_nyquist(ISwave_oc)
```

To save the workspace to a file:

``` matlab
>> save('EIS_20190515')
```

To load the workspace from a file:

``` matlab
>> load('EIS_20190515')
```

Complete documentation for each function is in the header of the relative file, which can be found [here](https://github.com/barnesgroupICL/Driftfusion/tree/2018-EIS).

This page can be accessed also in the Driftfusion [wiki](https://github.com/barnesgroupICL/Driftfusion/wiki/2018-Impedance-Spectroscopy-on-homojunction-model).

## Unit testing

``` matlab
>> runtests('examples_unit_test')
```
