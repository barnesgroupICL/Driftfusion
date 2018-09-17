function sol = changetau(sol_in, krad, taun, taup)

p = sol_in.p;
p.krad = krad;
p.taun = taun;
p.taup = taup;

sol = pindrift(sol_in, p);

end