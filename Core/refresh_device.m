function par = refresh_device(par)

par.xx = meshgen_x(par);
par.x_ihalf = getvarihalf(par.xx);
par.dev = build_device(par, 'iwhole');
par.dev_ihalf = build_device(par, 'ihalf');
% Get generation profiles
par.gx1 = generation(par, par.light_source1, par.laser_lambda1);
par.gx2 = generation(par, par.light_source2, par.laser_lambda2);

end