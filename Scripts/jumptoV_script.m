par = pc('./PEM_workshop_Input_files/doped_Shottky_junction.csv');

soleq = equilibrate(par);
dfplot.ELx(soleq.el);

%% Jump to V
% sol_relax = jumptoV(sol_ini, Vjump, tdwell, mobseti, Int, stabilise, accelerate)
sol_relax = jumptoV(soleq.el, 0.4, 1e-3, 0, 0, 0, 0);

%% Plot potential before and after jump
dfplot.Vx(soleq.el, [0, 1e-3]);
hold on
dfplot.Vx(sol_relax, [0, 1e-3]);
