function Gentot = beerlambert(par, x, source_type, laserlambda, figson)

%% Input arguments
% PAR - Parameters object
% SOURCE_TYPE - Currently either 'AM15' or 'laser'
% LASERLAMBDA - pulse wavelength only applied when 'laser' is chosen as the
% source. must be between 300 - 767 nm (owing to n & k data limits)
% FIGSON - logical toggling figures

%% Output arguments
% GENTOT - Generation rate integrated across the spectrum

% Author: Philip Calado, Imperial College London, 2018
% n & k loader taken from McGeehee group Stanford Transfer Matrix Code
% NOTE: In this version n & k change abruptly at interfaces

% Constants
h = 6.626e-34;
c = 2.998e8;

lambda = 300:767;  % wavelength range - currently limited by database MUST BE IN nm

switch par.side
    case 1
        layers = par.stack;
        layer_type = par.layer_type;
        pcum1 = [1,par.pcum+1];
        pcell = par.pcell;
    case 2
        % Flip properties for side 2
        layers = fliplr(par.stack);
        layer_type = fliplr(par.layer_type);
        x = x(end)-fliplr(par.xx);
        pcum1 = fliplr([1,par.pcum+1]);
        pcum1 = (pcum1(1)+1) - pcum1;
        pcell = fliplr(par.pcell);
end

switch source_type
    case 'AM15'      
        I0 = lightsource('AM15', lambda);
    case 'laser'
        % Throw up error if LASERLAMBDA is not within the acceptable range
        if laserlambda < lambda(1) || laserlambda > lambda(end)
            msg = ['Error in BEERLAMBERT - laser wavelength must be in the range ', num2str(lambda(1)), ' - ', num2str(lambda(end)), ' nm'];
            error(msg);
        end
        I0 = zeros(length(lambda));
        I0(laserlambda-lambda(1)) = 1e-3*par.pulsepow;   % convert to wcm-2
    otherwise
        warning('Light source unknown. Usinng AM1.5 as default')
        I0 = lightsource('AM15', lambda);
end

Eph = h*c./(lambda*1e-9);    % Photon energy in Joules

% Load in index of refraction for each material
n_layer = zeros(1,size(lambda,2));
k_layer = zeros(1,size(lambda,2));
n_layer_next = zeros(1,size(lambda,2));
k_layer_next = zeros(1,size(lambda,2));
n = zeros(length(x),size(lambda,2));
k = zeros(length(x),size(lambda,2));

try
    for i = 1:length(layers)    %size(layers,2)
        switch layer_type{1,i}
            case {'layer', 'active'}
                [n_layer, k_layer] = LoadRefrIndex(layers{i},lambda);
                n(pcum1(i):pcum1(i+1),:) = repmat(n_layer,pcell(i)+1, 1);
                k(pcum1(i):pcum1(i+1),:) = repmat(k_layer,pcell(i)+1, 1);       
            case 'junction'
                [n_layer, k_layer] = LoadRefrIndex(layers{i-1},lambda);
                [n_layer_next, k_layer_next] = LoadRefrIndex(layers{i+1},lambda);
                xprime = x(pcum1(i):pcum1(i+1)) - x(pcum1(i));
                xprime = xprime';
                k(pcum1(i):pcum1(i+1),:) = zeros(length(xprime),length(lambda));
%                 switch par.intgradfun
%                 	case 'linear'
%                         n(pcum1(i):pcum1(i+1),:) = n_layer + (n_layer_next - n_layer).*(xprime/xprime(end));
%                         k(pcum1(i):pcum1(i+1),:) = k_layer + (k_layer_next - k_layer).*(xprime/xprime(end));
%                     case 'erf'
%                         myerf = ((erf(2*pi*(xprime-par.d(i)/2)/(par.d(i)))+1)/2);
%                         n(pcum1(i):pcum1(i+1),:) = n_layer + (n_layer_next - n_layer).*myerf;
%                         k(pcum1(i):pcum1(i+1),:) = k_layer + (k_layer_next - k_layer).*myerf;
%                 end
        end
    end
catch
    error('Material name in stack does not match that in the Index of Refraction Library. Please check the names contained in the ''stack'' cell of the parameters object are correct.')
end

xcum = par.pcum;

alpha = zeros(length(x),size(lambda,2));
% calculate absorption coefficient from k
for i = 1:length(x)
    alpha(i,:) = 4*pi*k(i,:)./(lambda*1e-7);    % alpha in cm-1
end

I = zeros(length(x), length(lambda));

for i = 1:length(x)
for j = 1:length(lambda)
    if i == 1
        I(1,j) = I0(j);
        Gen(1, j) = I(1, j)*alpha(1, j)/Eph(j);
    else
        I(i,j) = I0(j)*exp(-alpha(i,j)*(x(i)-x(i-1)));
        Gen(i, j) = I(i, j)*alpha(i, j)/Eph(j);
        I0(j) = I(i, j);  % Set I0 for new layer
    end
end
end
Gentot = trapz(lambda, Gen, 2);

% Flip profiles if side 2
if par.side == 2
    I = flipud(I);
    Gen = flipud(Gen);
    Gentot = flipud(Gentot);
end

if figson == 1

figure(31)
surf(lambda, par.xx*1e7, I)
xlabel('Wavelength [nm]')
ylabel('Position [nm]')
zlabel('Intensity [Wcm-3]')

figure(32)
surf(lambda, par.xx*1e7, Gen)
xlabel('Wavelength [nm]')
ylabel('Position [nm]')
zlabel('Gen rate [cm-3s-1]')

figure(33)
plot(par.xx*1e7, Gentot)
xlabel('Position [nm]')
ylabel('Gen rate [cm-2s-1]')
xlim([0, par.xx(end)*1e7])
end

end