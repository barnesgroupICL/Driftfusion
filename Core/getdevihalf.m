function devihalf = getdevihalf(par)
dev = par.dev;
fnames = fieldnames(dev);

for i = 1:length(fnames)
    prop = char(fnames(i));
    eval(['devihalf.',prop,' = dev.',prop,'(1:end-1) + (dev.',prop,'(2:end) - dev.',prop,'(1:end-1))/2;']);
end



end