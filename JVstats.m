function stats = JVstats(JV)
% A function to pull statistics from a JV sweep using DOJV
% JV - a solution from DOJV

if isfield(JV, 'ill')
    
    if isfield(JV.ill, 'f')
        
        try
            p1 = find(JV.ill.f.Vapp >= 0);
            p1 = p1(1);
            stats.Jsc_f = JV.ill.f.Jtotr(p1, end);
            
        catch
            
            warning('No Jsc available- Vapp must pass through 0 to obtain Jsc')
            stats.Jsc_f = 0;
            
        end
        
        try
            p2 = find(JV.ill.f.Jtotr(:, end) >= 0);
            p2 = p2(1);
            stats.Voc_f = JV.ill.f.Vapp(p2);
            
        catch
            
            warning('No Voc available- try increasing applied voltage range')
            stats.Voc_f = 0;
            
        end
        
        if stats.Jsc_f ~= 0 && stats.Voc_f ~= 0
            pow_f = JV.ill.f.Jtotr(:,end).*JV.ill.f.Vapp';
            stats.mpp_f = min(pow_f);
            stats.FF_f = stats.mpp_f/(stats.Jsc_f*stats.Voc_f);
        end
        
    else
        stats.Jsc_f = nan;
        stats.Voc_f = nan;
        stats.mpp_f = nan;
        stats.FF_f = nan;
    end
    
    if isfield(JV.ill, 'r')
        
        try
            p1 = find(JV.ill.r.Vapp <= 0);
            p1 = p1(1);
            stats.Jsc_r = JV.ill.r.Jtotr(p1, end);
            
        catch
            
            warning('No Jsc available- Vapp must pass through 0')
            stats.Jsc_r = 0;
            
        end
        
        try
            p2 = find(JV.ill.r.Jtotr(:, end) <= 0);
            p2 = p2(1);
            stats.Voc_r = JV.ill.r.Vapp(p2);
            
        catch
            
            warning('No Voc available- try increasing applied voltage range')
            stats.Voc_r = 0;
            
        end
        
        if stats.Jsc_r ~= 0 && stats.Voc_r ~= 0
            pow_r = JV.ill.r.Jtotr(:,end).*JV.ill.r.Vapp';
            stats.mpp_r = min(pow_r);
            stats.FF_r = stats.mpp_r/(stats.Jsc_r*stats.Voc_r);
        end
        
    else
        stats.Jsc_r = nan;
        stats.Voc_r = nan;
        stats.mpp_r = nan;
        stats.FF_r = nan;
    end
else
    
    
end

end

