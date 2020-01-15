function ppos = getpointpos(xpos, xmesh)
% Obtain point position PPOS from x-position XPOS for mesh defined by XMESH

if xpos <= xmesh(1)
    ppos = 1;
else
    ppos = find(xmesh <= xpos);
    ppos = ppos(end);
end

end
