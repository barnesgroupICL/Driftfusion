function stats = JVstats(JV)

if isfield(JV, 'ill')
    
    if isfield(JV.ill, 'f')
        
        p = find(JV.ill.f.Vapp >= 0);
        p = p(1);
        stats.Jsc_f = JV.ill.f.Jtotr(p, end);
        
        p = find(JV.ill.f.Jtotr(:, end) >= 0);
        p = p(1);
        stats.Voc_f = JV.ill.f.Vapp(p);
        
        pow_f = JV.ill.f.Jtotr(:,end).*JV.ill.f.Vapp';
        stats.mpp_f = min(pow_f);
        stats.FF_f = stats.mpp_f/(stats.Jsc_f*stats.Voc_f);
    
    else
        stats.Jsc_f = nan;
        stats.Voc_f = nan;
        stats.mpp_f = nan;
        stats.FF_f = nan;
    end
    
    if isfield(JV.ill, 'r')
        
        p = find(JV.ill.r.Vapp <= 0);
        p = p(1);
        stats.Jsc_r = JV.ill.r.Jtotr(p, end);
        
        p = find(JV.ill.r.Jtotr(:, end) <= 0);
        p = p(1);
        stats.Voc_r = JV.ill.r.Vapp(p);
        
        pow_r = JV.ill.r.Jtotr(:,end).*JV.ill.r.Vapp';
        stats.mpp_r = min(pow_r);
        stats.FF_r = stats.mpp_r/(stats.Jsc_r*stats.Voc_r);
    
    else
        stats.Jsc_r = nan;
        stats.Voc_r = nan;
        stats.mpp_r = nan;
        stats.FF_r = nan;
    end
else
    
    
end

end

