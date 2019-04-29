% Copyright 2010 George F. Burkhard, Eric T. Hoke, Stanford University

%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

% This program calculates the field profile, exciton generation profile
% and generated current from the wavelength dependent complex indicies of
% refraction in devices using the transfer matrix method. It assumes the light source
% light source is in an n=1 environment (air) and that the first layer is
% a thick superstrate, so that incoherent reflection from the air/1st layer
% interface is taken into account before the coherent interference is
% calculated in the remaining layers. If there is no thick superstrate, 
% input 'Air' as the first layer so that the incoherent reflection calculates
% to 0. 
% The program
% also returns the calculated short circuit current for the device in 
% mA/cm^2 if the device had an IQE of 100%.

% The procedure was adapted from J. Appl. Phys Vol 86, No. 1 (1999) p.487
% and JAP 93 No. 7 p. 3693.

% George Burkhard and Eric Hoke February 5, 2010
% When citing this work, please refer to:
%
% G. F. Burkhard, E. T. Hoke, M. D. McGehee, Adv. Mater., 22, 3293.
% Accounting for Interference, Scattering, and Electrode Absorption to Make
% Accurate Internal Quantum Efficiency Measurements in Organic and Other 
% Thin Solar Cells

% Modifications:
% 3/3/11 Parastic absorption (parasitic_abs) is calculated and made
% accessable outside of script.  
% 3/5/15 Improved precision of fundamental constants


function [Gx1, LasGx1] = TMPC1(layers, thicknesses, actlyr1, actlyr2, laserlambda, pulsepow)
%------------BEGIN USER INPUT PARAMETERS SPECITIFCATION---------------
%
lambda=300:767; % Wavelengths over which field patterns are calculated
stepsize = 1;   % The electric field is calculated at a latice of points (nm)
                % in the device cross section seperated by this distance

% plotWavelengths specifies which wavelengths to plot when plotting E-field
% intensity distributions (figure 1). Specify values by adding wavelength
% values to the array. Values must be within the range of calculated
% wavelenths (ie. must be an element of the lambda array). All wavelengths
% are in nanometers.
plotWavelengths = laserlambda;

thicknesses = thicknesses*1e7;      % Convert cm to nm

% Specify Layers in device (an arbitrary number of layers is permitted) and 
% thicknesses.
%
% Change these arrays to change the order or number of layers and/or
% thickness of layers. List the layers in the order that they appear in the
% device starting with the side the light is incident on.  THE NAMES OF THE
% LAYERS MUST CORRESPOND TO THE NAMES OF THE MATERIALS IN THE INDEX OF 
% REFRACTION LIBRARY FILE, 'Index_of_Refraction_library.xls'. The first
% layer must be the transparent substrate (glass) or 'Air' if the active
% layers are on the reflective electrode (rather than transparent electrode) side 
% of the device.  The layer thicknesses are in nanometers.

%Test Values

%layers = {'Air' 'P3HT'}; % Names of layers of materials starting from side light is incident from
%thicknesses = [100 2000];  % thickness of each corresponding layer in nm (thickness of the first layer is irrelivant)
%actlyr1 = 2; % index of material layer where photocurrent is generated
%actlyr2 = 2; % Mod to allow gen profile for mpTiO2 layer- set to same as activelayer 1 for single layer. PC

%Standard Bottom-Cat Planar Device Values

%layers = {'SiO2' 'ITOsorizon' 'TiO2' 'MAPICl' 'Spiro' 'Au_TF'}; % Names of layers of materials starting from side light is incident from
%thicknesses = [100 100 50 400 100 100];  % thickness of each corresponding layer in nm (thickness of the first layer is irrelivant)
%actlyr1 = 4; % index of material layer where photocurrent is generated
%actlyr2 = 4; % Mod to allow gen profile for mpTiO2 layer- set to same as activelayer 1 for single layer. PC

%Standard Bottom-Cat MP Device Values

%layers = {'SiO2' 'ITOsorizon' 'TiO2' 'mpTi02/MAPICl' 'MAPICl' 'Spiro' 'Au_TF'}; % Names of layers of materials starting from side light is incident from
%thicknesses = [100 110 50 300 200 100 100];  % thickness of each corresponding layer in nm (thickness of the first layer is irrelivant)%
%actlyr1 = 4; % index of material layer where photocurrent is generated
%actlyr2 = 5; % Mod to allow gen profile for mpTiO2 layer- set to same as activelayer 1 for single layer. PC


%Standard Top-Cat Device Values

%layers = {'SiO2' 'ITOsorizon' 'PEDOT' 'MAPICl' 'PCBM' 'Au_TF'}; % Names of layers of materials starting from side light is incident from
%thicknesses = [100 110 30 350 100 100];  % thickness of each corresponding layer in nm (thickness of the first layer is irrelivant)
%actlyr1 = 4; % index of material layer where photocurrent is generated
%actlyr2 = 4; % Mod to allow gen profile for mpTiO2 layer- set to same as activelayer 1 for single layer. PC


%Double-Transparent electrode laminate Bottom-Cat Device Values

%Glass Side
%layers = {'SiO2' 'ITOsorizon' 'TiO2' 'mpTi02/MAPICl' 'MAPICl' 'Spiro' 'PEDOT'}; % Names of layers of materials starting from side light is incident from
%thicknesses = [100 100 50 200 200 100 100];  % thickness of each corresponding layer in nm (thickness of the first layer is irrelivant)
%actlyr1 = 4; % index of material layer where photocurrent is generated
%actlyr2 = 5; % Mod to allow gen profile for mpTiO2 layer- set to same as activelayer 1 for single layer. PC

%Laminate Side
%layers = {'PEDOT' 'Spiro' 'MAPICl' 'mpTi02/MAPICl' 'TiO2' 'ITOsorizon' 'SiO2'};
%thicknesses = [100 100 200 200 50 100 100];
%actlyr1 = 3; % index of material layer where photocurrent is generated
%actlyr2 = 4; % Mod to allow gen profile for mpTiO2 layer- set to same as activelayer 1 for single layer. PC

% Set plotGeneration to 'true' if you want to plot generation rate as a
% function of position in the device and output the calculated short circuit current
% under AM1.5G illumination (assuming 100% internal quantum efficiency)
plotGeneration = true;
plotGenLaser = true;

%
%------------END USER INPUT PARAMETERS SPECIFICATION-------------------

% Load in index of refraction for each material
n = zeros(size(layers,2),size(lambda,2));
for index = 1:size(layers,2)
    n(index,:) = LoadRefrIndex(layers{index},lambda);
end
t = thicknesses;

% Constants
h = 6.62606957e-34; 	% Js Planck's constant
c = 2.99792458e8;	% m/s speed of light
q = 1.60217657e-19;	% C electric charge

% Calculate Incoherent power transmission through substrate
% See Griffiths "Intro to Electrodynamics 3rd Ed. Eq. 9.86 & 9.87
T_glass=abs(4*1*n(1,:)./(1+n(1,:)).^2); 
R_glass=abs((1-n(1,:))./(1+n(1,:))).^2;

% Calculate transfer matrices, and field at each wavelength and position
t(1)=0;
t_cumsum=cumsum(t);
x_pos=(stepsize/2):stepsize:sum(t); %positions to evaluate field

%x_mat specifies what layer number the corresponding point in x_pos is in:
x_mat= sum(repmat(x_pos,length(t),1)>repmat(t_cumsum',1,length(x_pos)),1)+1;

%Mod to allow 2 active layer PC
x_mat1= sum(repmat(x_pos,length(t),1)>repmat(t_cumsum',1,length(x_pos)),1)+1; 
x_mat2= sum(repmat(x_pos,length(t),1)>repmat(t_cumsum',1,length(x_pos)),1)+1; 

R=lambda*0;
E=zeros(length(x_pos),length(lambda));
ESqu = zeros(length(x_pos),length(lambda));

for l = 1:length(lambda)
    % Calculate transfer matrices for incoherent reflection and transmission at the first interface
    
    S=I_mat(n(1,l),n(2,l));
    for matindex=2:(length(t)-1)
        S=S*L_mat(n(matindex,l),t(matindex),lambda(l))*I_mat(n(matindex,l),n(matindex+1,l));
    end
    R(l)=abs(S(2,1)/S(1,1))^2; %JAP Vol 86 p.487 Eq 9 Power Reflection from layers other than substrate
    T(l)=abs(2/(1+n(1,l)))/sqrt(1-R_glass(l)*R(l)); %Transmission of field through glass substrate Griffiths 9.85 + multiple reflection geometric series

    % Calculate all other transfer matrices
    for material = 2:length(t) 
        xi=2*pi*n(material,l)/lambda(l);
        dj=t(material);
        x_indices=find(x_mat == material); %indices of points which are in the material layer considered
        x=x_pos(x_indices)-t_cumsum(material-1); %distance from interface with previous layer
        
        % Calculate S matrices (JAP Vol 86 p.487 Eq 12 and 13)
        S_prime=I_mat(n(1,l),n(2,l));
        for matindex=3:material
            S_prime=S_prime*L_mat(n(matindex-1,l),t(matindex-1),lambda(l))*I_mat(n(matindex-1,l),n(matindex,l));
        end
        S_doubleprime=eye(2);
        for matindex=material:(length(t)-1)
            S_doubleprime=S_doubleprime*I_mat(n(matindex,l),n(matindex+1,l))*L_mat(n(matindex+1,l),t(matindex+1),lambda(l));
        end
        % Normalized Field profile (JAP Vol 86 p.487 Eq 22)
        E(x_indices,l)=T(l)*(S_doubleprime(1,1)*exp(-1i*xi*(dj-x))+S_doubleprime(2,1)*exp(1i*xi*(dj-x))) ./(S_prime(1,1)*S_doubleprime(1,1)*exp(-1i*xi*dj)+S_prime(1,2)*S_doubleprime(2,1)*exp(1i*xi*dj));
        
    end 
end

%Get the field intensity
ESqu = abs(E).^2;

% Overall Reflection from device with incoherent reflections at first
% interface (typically air-glass)
Reflection=R_glass+T_glass.^2.*R./(1-R_glass.*R);

% Plots electric field intensity |E|^2 vs position in device for
% wavelengths specified in the initial array, plotWavelengths. 

    figure(101)
    plotString = '';
    legendString = cell(1,size(plotWavelengths,2));

ESquLaser = zeros(length(x_pos),size(plotWavelengths,2));

for index=1:size(plotWavelengths,2)

    plotString = strcat(plotString, ['x_pos,abs(E(:,', num2str(find(lambda == plotWavelengths(index))), ').^2),']);
    legendString{index} = [num2str(plotWavelengths(index)), ' nm'];
    
    ESquLaser(:,index) = abs(E(:,(find(lambda == plotWavelengths(index))))).^2;
  
end

    eval(['plot(', plotString, '''LineWidth'',2)'])
    axislimit1=axis;

% Draws vertical lines at each material boundary in the stack and labels
% each layer
for matindex=2:length(t)
    line([sum(t(1:matindex)) sum(t(1:matindex))],[0 axislimit1(4)]);
    text((t_cumsum(matindex)+t_cumsum(matindex-1))/2,0,layers{matindex},'HorizontalAlignment','center','VerticalAlignment','bottom')
end

title('E-field instensity in device');
xlabel('Position in Device (nm)');
ylabel('Normalized Electric field intensity |E|^2');
legend(legendString);


% Absorption coefficient in cm^-1 (JAP Vol 86 p.487 Eq 23)
a=zeros(length(t),length(lambda));
for matindex=2:length(t)
     a(matindex,:)=4*pi*imag(n(matindex,:))./(lambda*1e-7);
end

% Plots normalized intensity absorbed /cm3-nm at each position and
% wavelength as well as the total reflection expected from the device
% (useful for comparing with experimentally measured reflection spectrum)

figure(101)

Absorption=zeros(length(t),length(lambda));
plotString = '';
for matindex=2:length(t)
    Pos=find(x_mat == matindex);
    AbsRate=repmat(a(matindex,:).*real(n(matindex,:)),length(Pos),1).*(abs(E(Pos,:)).^2);
    Absorption(matindex,:)=sum(AbsRate,1)*stepsize*1e-7;
    plotString = strcat(plotString, ['lambda,Absorption(', num2str(matindex), ',:),']);
end

    eval(['plot(', plotString, 'lambda,Reflection,''LineWidth'',2)'])
    title('Fraction of Light absorbed or reflected');
    xlabel('Wavelength (nm)');
    ylabel('Light Intensity Fraction');
    legend(layers{2:size(layers,2)}, 'Reflectance');


% Plots generation rate as a function of position in the device and
% calculates Jsc
if plotGeneration == true
    
    % Load in 1sun AM 1.5 solar spectrum in mW/cm2nm
    AM15_data=xlsread('AM15.xls');
    AM15=interp1(AM15_data(:,1), AM15_data(:,2), lambda, 'linear', 'extrap');
    
    figure(101)
    % Energy dissipation mW/cm3-nm at each position and wavelength (JAP Vol
    % 86 p.487 Eq 22)
    ActivePos1=find(x_mat1 == actlyr1);
    ActivePos2=find(x_mat2 == actlyr2);
    
    Q1=repmat(a(actlyr1,:).*real(n(actlyr1,:)).*AM15,length(ActivePos1),1).*(abs(E(ActivePos1,:)).^2);
    Q2=repmat(a(actlyr2,:).*real(n(actlyr2,:)).*AM15,length(ActivePos2),1).*(abs(E(ActivePos2,:)).^2);

    % Exciton generation rate per second-cm3-nm at each position and wavelength
    Gxl1=(Q1*1e-3).*repmat(lambda*1e-9,length(ActivePos1),1)/(h*c);
    Gxl2=(Q2*1e-3).*repmat(lambda*1e-9,length(ActivePos2),1)/(h*c);   
    
    if length(lambda)==1
        lambdastep= 1;
    else
        lambdastep=(max(lambda)-min(lambda))/(length(lambda)-1);
    end
       
    Gx1=sum(Gxl1,2)*lambdastep; % Exciton generation rate as a function of position/(sec-cm^3), active layer 1
    Gx2=sum(Gxl2,2)*lambdastep; % Exciton generation rate as a function of position/(sec-cm^3), active layer 2
    
    plot(x_pos(ActivePos1),Gx1,x_pos(ActivePos2),Gx2,'LineWidth',2)
    %plot(x_pos(ActivePos2),Gx2,'LineWidth',2)
    axislimit3=axis;
    axis([axislimit1(1:2) axislimit3(3:4)])
    
    % inserts vertical lines at material boundaries
    for matindex=2:length(t)
        line([sum(t(1:matindex)) sum(t(1:matindex))],[0 axislimit3(4)]);
        text((t_cumsum(matindex)+t_cumsum(matindex-1))/2,0,layers{matindex},'HorizontalAlignment','center','VerticalAlignment','bottom')
    end

    title('Generation Rate in Device') 
    xlabel('Position in Device (nm)');
    ylabel('Generation rate /(sec-cm^3)');
    
    figure(102)
    surf(Gxl2)
    
    figure(104)
    surf(Q2)
    % outputs predicted Jsc under AM1.5 illumination assuming 100% internal
    % quantum efficiency at all wavelengths
    Jsc1=sum(Gx1)*stepsize*1e-7*q*1e3; %in mA/cm^2
    Jsc2=sum(Gx2)*stepsize*1e-7*q*1e3; %in mA/cm^2

    % Calculate parasitic absorption
    parasitic_abs=(1-Reflection-Absorption(actlyr1,:))';
    parasitic_abs=(1-Reflection-Absorption(actlyr2,:))';
     
    
   %Tranpsose lambda for future use- seems to be necessary if code is
   %script and not function- PC
   
   %lambda = transpose(lambda)
end

% Mod Plots generation rate as a function of position in the device and
% calculates Jsc for single wavelength/Laser- may be easier ways to do
% this!

LasPow_data = (280:1:2000)';
LasPow_data(:, 2) = 0;
LasPow_data(laserlambda - LasPow_data(1), 2) = pulsepow;

if plotGenLaser == true
      
    LasPow=interp1(LasPow_data(:,1), LasPow_data(:,2), lambda, 'linear', 'extrap');

    figure(101)
    % Energy dissipation mW/cm3-nm at each position and wavelength (JAP Vol
    % 86 p.487 Eq 22)
    ActivePos1=find(x_mat1 == actlyr1);
    ActivePos2=find(x_mat2 == actlyr2);
    
    LasQ1=repmat(a(actlyr1,:).*real(n(actlyr1,:)).*LasPow,length(ActivePos1),1).*(abs(E(ActivePos1,:)).^2);
    LasQ2=repmat(a(actlyr2,:).*real(n(actlyr2,:)).*LasPow,length(ActivePos2),1).*(abs(E(ActivePos2,:)).^2);

    % Exciton generation rate per second-cm3-nm at each position and wavelength
    LasGxl1=(LasQ1*1e-3).*repmat(lambda*1e-9,length(ActivePos1),1)/(h*c);
    LasGxl2=(LasQ2*1e-3).*repmat(lambda*1e-9,length(ActivePos2),1)/(h*c);

    if length(lambda)==1
        lambdastep= 1;
    else
        lambdastep=(max(lambda)-min(lambda))/(length(lambda)-1);
    end
    
    LasGx1=sum(LasGxl1,2)*lambdastep; % Exciton generation rate as a function of position/(sec-cm^3), active layer 1
    LasGx2=sum(LasGxl2,2)*lambdastep; % Exciton generation rate as a function of position/(sec-cm^3), active layer 2

    plot(x_pos(ActivePos1),LasGx1,x_pos(ActivePos2),LasGx2,'LineWidth',2);
    %plot(x_pos(ActivePos2),Gx2,'LineWidth',2)
    
    axislimit3=axis;
    axis([axislimit1(1:2) axislimit3(3:4)]);

    % inserts vertical lines at material boundaries
    for matindex=2:length(t)
        line([sum(t(1:matindex)) sum(t(1:matindex))],[0 axislimit3(4)]);
        text((t_cumsum(matindex)+t_cumsum(matindex-1))/2,0,layers{matindex},'HorizontalAlignment','center','VerticalAlignment','bottom');
    end

    title('Generation Rate in Device'); 
    xlabel('Position (points)');
    ylabel('Generation rate /(sec-cm^3)');
    grid on

    % outputs predicted Jsc under AM1.5 illumination assuming 100% internal
    % quantum efficiency at all wavelengths
    LasJsc1=sum(LasGx1)*stepsize*1e-7*q*1e3; %in mA/cm^2
    LasJsc2=sum(LasGx2)*stepsize*1e-7*q*1e3; %in mA/cm^2

    % Calculate parasitic absorption
    Las_parasitic_abs=(1-Reflection-Absorption(actlyr1,:))';
    Las_parasitic_abs=(1-Reflection-Absorption(actlyr2,:))';

    % sends absorption, reflection, and wavelength data to the workspace
    TMsol.LasGx1 = LasGx1;
    TMsol.LasGx2 = LasGx2;
    TMsol.AM15Gx1 = Gx1;
    TMsol.AM15Gx2 = Gx2;
    assignin('base','TMsol',TMsol);
    
    %assignin('base','absorption',Absorption');
    %assignin('base','reflection',Reflection');
    %assignin('base','Las_parasitic_abs',Las_parasitic_abs);
    %assignin('base','lambda',lambda');
    %assignin('base','LasPow',LasPow);
    %assignin('base','Las_Gen_Rate_Active_Layer1',LasGx1);
    %assignin('base','Las_Gen_Rate_Active_Layer2',LasGx2);
    %assignin('base','E',E);
    %assignin('base','ESqu', ESqu);
    
end  

%------------------- Helper Functions ------------------------------------
% Function I_mat
% This function calculates the transfer matrix, I, for reflection and
% transmission at an interface between materials with complex dielectric 
% constant n1 and n2.
function I = I_mat(n1,n2)
r=(n1-n2)/(n1+n2);
t=2*n1/(n1+n2);
I=[1 r; r 1]/t;

% Function L_mat
% This function calculates the propagation matrix, L, through a material of
% complex dielectric constant n and thickness d for the wavelength lambda.
function L = L_mat(n,d,lambda)
xi=2*pi*n/lambda;
L=[exp(-1i*xi*d) 0; 0 exp(1i*xi*d)];

% Function LoadRefrIndex
% This function returns the complex index of refraction spectra, ntotal, for the
% material called 'name' for each wavelength value in the wavelength vector
% 'wavelengths'.  The material must be present in the index of refraction
% library 'Index_of_Refraction_library.xls'.  The program uses linear
% interpolation/extrapolation to determine the index of refraction for
% wavelengths not listed in the library.
function ntotal = LoadRefrIndex(name,wavelengths)

%Data in IndRefr, Column names in IndRefr_names
[IndRefr,IndRefr_names]=xlsread('Index_of_Refraction_library.xls');

% Load index of refraction data in spread sheet, will crash if misspelled
file_wavelengths=IndRefr(:,strmatch('Wavelength',IndRefr_names));
n=IndRefr(:,strmatch(strcat(name,'_n'),IndRefr_names));
k=IndRefr(:,strmatch(strcat(name,'_k'),IndRefr_names));

% Interpolate/Extrapolate data linearly to desired wavelengths
n_interp=interp1(file_wavelengths, n, wavelengths, 'linear', 'extrap');
k_interp=interp1(file_wavelengths, k, wavelengths, 'linear', 'extrap');

%Return interpolated complex index of refraction data
ntotal = n_interp+1i*k_interp; 
