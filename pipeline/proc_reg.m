function [breathing_reg,all5_reg] = proc_reg(regfilename)

    regmat = load(regfilename);
    reg = regmat.reg;

    [v,idcs] = max(corr(reg(:,1:16),reg(:,17:end)),'all');
    [v2, idc_reg] = max(v);

    breathing_reg = squeeze(reg(:,16 + idc_reg));
    all5_reg = reg(:,17:end);

    filenamebase = strsplit(regfilename,'.');
    filenamebase_base = filenamebase{1};

    breath_regmat.reg = normalize(breathing_reg,1);
    save(strjoin({filenamebase_base,'breath_reg.mat'},'_'),'-struct','breath_regmat');
    
    

end

