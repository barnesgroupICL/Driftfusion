function spvdat = spvana(spvsol)
% A function for analysing the solution SPVSOL from DOSPV2

%% Extract initial parameters and meshes
par = spvsol.ill.par;
t = spvsol.ill.t;
x = spvsol.ill.x;

%% Illuminated step
% integrated charge density as a function of time
rho.ill = dfana.calcrho(spvsol.ill);
sigmat.ill = trapz(x, rho.ill, 2);

% Get illuminated field
F.ill = dfana.calcF(spvsol.ill);
% Field at final position x=d as a function of time
Ft.ill = F.ill(:,end);
% Ionic component of field
Fion.ill = dfana.calcFion(spvsol.ill);
Fiont.ill = Fion.ill(:,end);

% Get illuminated potential
V.ill = (spvsol.ill.u(:,:,4));
% Direct readout of potential at x=d as a fucntion of time
Vt.ill = V.ill(:,end);

% Calculate Vion
Vion.ill = dfana.calcVion(spvsol.ill);
% Vion at x=d
Viont.ill = Vion.ill(:,end);

%% Dark decay time step
% integrated charge density as a function of time
rho.dk = dfana.calcrho(spvsol.dk);
sigmat.dk = trapz(x, rho.dk, 2);

% Get dark field
F.dk = dfana.calcF(spvsol.dk);
% Field at final position x=d as a function of time
Ft.dk = F.dk(:,end);
% Ionic component of field
Fion.dk = dfana.calcFion(spvsol.dk);
Fiont.dk = Fion.dk(:,end);

% Get dark potential
V.dk = (spvsol.dk.u(:,:,4));
% Direct readout of potential at x=d as a fucntion of time
Vt.dk = V.dk(:,end);

% Calculate Vion
Vion.dk = dfana.calcVion(spvsol.dk);
% Vion at x=d
Viont.dk = Vion.dk(:,end);

%% Concatenate vectors and write to output structure SPVDAT
spvdat.t = [t, t+t(end)];
spvdat.sigmat = [sigmat.ill; sigmat.dk];
spvdat.Ft = [Ft.ill; Ft.dk];
spvdat.Vt = [Vt.ill; Vt.dk];
spvdat.Fiont = [Fiont.ill; Fiont.dk];
spvdat.Viont = [Viont.ill; Viont.dk];
spvdat.Velt = spvdat.Vt - spvdat.Viont;

%% Get the deltas
spvdat.deltaFt = spvdat.Ft - spvdat.Ft(1);
spvdat.deltaVt = spvdat.Vt - spvdat.Vt(1);
spvdat.deltaVelt = spvdat.Velt - spvdat.Velt(1);
spvdat.deltaFiont = spvdat.Fiont - spvdat.Fiont(1);
spvdat.deltaViont = spvdat.Viont - spvdat.Viont(1);
spvdat.SPV = spvdat.Velt - spvdat.Velt(1);

%% Plot the outputs
figure(300)
plot(spvdat.t, spvdat.sigmat)
xlabel('Time [s]')
ylabel('Integrated charge desnity, \sigma [cm-2]')

figure(301)
plot(spvdat.t, spvdat.deltaFt)
xlabel('Time [s]')
ylabel('\Delta Field at x=d [Vcm-1]')

figure(302)
plot(spvdat.t, spvdat.deltaVt, spvdat.t, spvdat.deltaViont, spvdat.t, spvdat.SPV)
xlabel('Time [s]')
ylabel('\Delta Electric potential at x=d [V]')
legend('Vtotal', 'Vion', 'Vtotal - Vion')

end