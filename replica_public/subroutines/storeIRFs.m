%-------------------------------------------------------------------------%
% Store IRFs for group analysis
%-------------------------------------------------------------------------%

if modelOpt.dummies
    stoBVAR_d.(dataOpt.units{1}) = BVAR_d.IRFs;
    if miscOpt.SSAirf
        for i = 1:numel(ssaOpt.impulse)
            impulse0 = ssaOpt.impulse(i);
        stoBVARssa_d.(dataOpt.units{1}) = BVAR_d.SSA.(impulse0{:}).chIRFs;
        end
    end
    if plotOpt.homogeneity
        stoBVARhomogeneity_d.(dataOpt.units{1}) = BVAR_d.IRFsulm;
    end
end

if modelOpt.triangular
    stoBVAR_t.(dataOpt.units{1}) = BVAR_t.IRFs;
    if miscOpt.SSAirf
        stoBVARssa_t.(dataOpt.units{1}) = BVARssa_t.chIRFs;
    end
    if plotOpt.homogeneity
        stoBVARhomogeneity_t.(dataOpt.units{1}) = BVAR_t.IRFsulm;
    end
end

if z == 1

    if modelOpt.dummies
    save([stoFolder,'/stoBVAR_d.mat'],'-struct','stoBVAR_d',char(dataOpt.units));
        if miscOpt.SSAirf
            save([stoFolder,'/stoBVARssa_d.mat'],'-struct','stoBVARssa_d',char(dataOpt.units));
        end
        if plotOpt.homogeneity
            save([stoFolder,'/stoBVARhomogeneity_d.mat'],'-struct','stoBVARhomogeneity_d',char(dataOpt.units));
        end
    end

    if modelOpt.triangular
    save([stoFolder,'/stoBVAR_t.mat'],'-struct','stoBVAR_t',char(dataOpt.units));
        if miscOpt.SSAirf
            save([stoFolder,'/stoBVARssa_t.mat'],'-struct','stoBVARssa_t',char(dataOpt.units));
        end
        if plotOpt.homogeneity
            save([stoFolder,'/stoBVARhomogeneity_t.mat'],'-struct','stoBVARhomogeneity_t',char(dataOpt.units));
        end
    end
    
    sample = [dataOpt.units{1},data.dates(1),data.dates(end)];
    save([stoFolder,'/sample.mat'],'sample');

else

    if modelOpt.dummies
    save([stoFolder,'/stoBVAR_d.mat'],'-struct','stoBVAR_d',char(dataOpt.units),'-append');
        if miscOpt.SSAirf
            save([stoFolder,'/stoBVARssa_d.mat'],'-struct','stoBVARssa_d',char(dataOpt.units),'-append');
        end
        if plotOpt.homogeneity
            save([stoFolder,'/stoBVARhomogeneity_d.mat'],'-struct','stoBVARhomogeneity_d',char(dataOpt.units),'-append');
        end
    end

    if modelOpt.triangular
    save([stoFolder,'/stoBVAR_t.mat'],'-struct','stoBVAR_t',char(dataOpt.units),'-append');
        if miscOpt.SSAirf
            save([stoFolder,'/stoBVARssa_t.mat'],'-struct','stoBVARssa_t',char(dataOpt.units),'-append');
        end
        if plotOpt.homogeneity
            save([stoFolder,'/stoBVARhomogeneity_t.mat'],'-struct','stoBVARhomogeneity_t',char(dataOpt.units),'-append');
        end
    end

    load([stoFolder,'/sample.mat'])
    sample{end+1,1} = dataOpt.units{1};
    sample{end,2} = data.dates(1);
    sample{end,3} = data.dates(end);
    save([stoFolder,'/sample.mat'],'sample');
    
end

% Save model details (do it once)
if z == 1
    save([stoFolder,'/modelOpt'],'-struct','modelOpt')
    save([stoFolder,'/dataOpt'],'-struct','dataOpt')
    save([stoFolder,'/plotOpt'],'-struct','plotOpt')
    save([stoFolder,'/miscOpt'],'-struct','miscOpt')
    save([stoFolder,'/ssaOpt'],'-struct','ssaOpt')
end