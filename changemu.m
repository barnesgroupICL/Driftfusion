function sol = changemu(sol_in, mue, muh)

p = sol_in.p;
p.mue = mue;
p.muh = muh;

sol = pindrift(sol_in, p);

end