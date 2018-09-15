function sol = changetau(sol_in, taun, taup)

p = sol_in.p;
p.taun = taun;
p.taup = taup;

sol = pindrift(sol_in, p);

end