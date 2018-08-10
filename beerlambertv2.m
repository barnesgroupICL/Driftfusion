function gx = beerlambertv2(params, x, source_type)
% n & k loader taken from Stanford Transfer Matrix Code

% Constants
h = 6.626e-34;
c = 2.998e8;

layers = params.stack;
lambda = 300:767;          % wavelength range - currently limited by database MUST BE IN nm
Eph = h*c./(lambda*1e-9);    % Photon energy in Joules

% Load in index of refraction for each material
n = zeros(size(layers,2),size(lambda,2));
k = zeros(size(layers,2),size(lambda,2));
try
    for index = 1:size(layers,2)
        [n(index,:), k(index,:)] = LoadRefrIndex(layers{index},lambda);
    end
catch
    error('Material name in stack does not match that in the Index of Refraction Library. Please check the names contained in the ''stack'' cell of the parameters object are correct.')
end
% temporarily set I0 to one for all wavelengths
%I0 = ones(1, length(lambda))*10e-6;

I0 = lightsource('AM15', lambda);

xarr = [0, params.d];
xcum = [0, params.dcum];

% calculate absorption coefficient from k
for i = 1:length(layers)

    alpha(i,:) = 4*pi*k(i,:)./(lambda*1e-7);    % alpha in cm-1
end

I = zeros(length(x), length(lambda));

% Iterate across layers
for i =1:length(xarr)-1
  
    % Iterate across wavelengths
    for j = 1:length(lambda)

    
        x1 = sum(xarr(1:i));
        x2 = sum(xarr(1:(i+1)));
        
        % Set p1 to 1 on first loop to avoid problems with mesh not
        % starting at 0
        if i == 1
            
            p1 = 1;
            
        else
            
            p1 = find(x <= x1);
            p1 = p1(end);
        
        end
        
        p2 = find(x <= x2);
        p2 = p2(end);       
                       
        I(p1:p2, j) = I0(j).*exp(-alpha(i, j).*(x(p1:p2)-x(p1)));
        
        Gen(p1:p2, j) = I(p1:p2, j)*alpha(i, j)/Eph(j);
end
        I0 = I(p2, :);  % Set I0 for new layer
end


Gentot = trapz(lambda, Gen, 2);

figure(1)
surf(lambda, x, I)
xlabel('Wavelength [nm]')
ylabel('Position [nm]')
zlabel('Intensity [Wcm-3]')

figure(2)
surf(lambda, x, Gen)
xlabel('Wavelength [nm]')
ylabel('Position [nm]')
zlabel('Gen rate [cm-3s-1]')

figure(3)
plot(x, Gentot)
xlabel('Position [nm]')
ylabel('Gen rate [cm-2s-1]')

end