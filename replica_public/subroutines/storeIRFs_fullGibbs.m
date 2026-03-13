%-------------------------------------------------------------------------%
% Store IRFs for group analysis
%-------------------------------------------------------------------------%

impulse0 = ssaOpt.impulse(1);
stoSSA.(dataOpt.units{1}) = SSAout.SSAgroup.(impulse0{:}).chIRFs;

if z == 1

    save([stoFolder,'/stoSSA.mat'],'-struct','stoSSA',char(dataOpt.units));
    sample = [dataOpt.units{1},data.dates(1),data.dates(end)];
    save([stoFolder,'/sample.mat'],'sample');

else

    save([stoFolder,'/stoSSA.mat'],'-struct','stoSSA',char(dataOpt.units),'-append');

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