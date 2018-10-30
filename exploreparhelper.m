function par = exploreparhelper(par, parname, parvalue)

% takes parameter set and sets parname to parvalue- workaround for parallel
% computing loops
eval(['par.',parname,'=parvalue']);

end