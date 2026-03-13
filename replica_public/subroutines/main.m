%-------------------------------------------------------------------------%
% MATLAB TOOLBOX                                                          %
% By Riccardo Degasperi                                                   %
%                                                                         %
% riccardo.degasperi@bancaditalia.it                                      %
%-------------------------------------------------------------------------%

addpath([pwd '/subroutines/'])                   %path to child functions
addpath([pwd '/data/'])                          %path to data folder

%-------------------------------------------------------------------------%
% Data import
disp('--> Preliminary operations...')

[data,dataOpt] = dataLoaderBVAR(modelOpt,dataOpt,miscOpt);

%-------------------------------------------------------------------------%
% Checks and patches
checks

%-------------------------------------------------------------------------%
% Plot series
if plotOpt.plotTimeSeries

    plotSeries(data,dataOpt,plotOpt)
    
end


%-------------------------------------------------------------------------%
% Make X and Y matrices
[data,dims,hp,dataOpt] = makeXY(data,hp,modelOpt,dataOpt);
plotOpt.ylab = dataOpt.ylab;


%-------------------------------------------------------------------------%
% Full Gibbs sampler for Bayesian VAR (dummies or triangular)
if modelOpt.fullGibbs

out = fullGibbsSampler(data,dims,hp,modelOpt,dataOpt,miscOpt);

% Display priors and save options to a text file
filename = 'fullGibbs_triangular';
printDetails(out,filename,data,hp,modelOpt,dataOpt,miscOpt,plotOpt)

% Plot IRFs
chartName = [char(dataOpt.units),'_FullGibbs'];
BVARplots(out.IRFsulm,chartName,dataOpt,plotOpt);

% Forecast Error Variance Decomposition
if miscOpt.FEVD
    out.FEVD = getFEVD(out,dims,modelOpt);

    if exist('group')
        prefix = [group{z},'_'];
    else
        prefix = '';
    end

    if strcmp(plotOpt.plotType,'allShocks')
        customLegend = {};
        chartName = [prefix,'FullGibbs_FEVDall'];
        plotAllFEVD(out.FEVD.FEVDmean,dataOpt,plotOpt,customLegend,chartName);
    elseif strcmp(plotOpt.plotType,'singleShock')
        chartName = [prefix,'FullGibbs_FEVDchart'];
        plotSingleFEVD(out.FEVD,dataOpt,plotOpt,chartName);
    else
        error('ERROR: FEVD is available for either all or one shocks.')
    end

end

% Historical Decomposition
if miscOpt.HD
    out.HD = getHD(out,data,dims,modelOpt);

    if exist('group')
        prefix = [group{z},'_'];
    else
        prefix = '';
    end

    if strcmp(plotOpt.plotType,'allShocks')
        customLegend = {};
        chartName = [prefix,'FullGibbs_HDall'];
        plotAllHD(out.HD,data,dataOpt,plotOpt,customLegend,chartName);
    elseif strcmp(plotOpt.plotType,'singleShock')
        chartName = [prefix,'FullGibbs_HDchart'];
        plotSingleHD(out.HD,data,dataOpt,plotOpt,chartName);
    else
        error('ERROR: FEVD is available for either all or one shocks.')
    end
end 

%-------------------------------------------------------------------------%
% Full Gibbs sampler for in-sample structural scenario analysis
% (Antolin-Diaz, Petrella, Rubio-Ramirez, 2021; Breitenlechner, Georgiadis, Schumann, 2022)
elseif modelOpt.SSAfullIRF

    SSAout = SSAfullGibbsIRF(data,dims,hp,modelOpt,dataOpt,miscOpt,plotOpt,ssaOpt);

else
%-------------------------------------------------------------------------%
% Estimating BVAR via dummy observations
% (Banbura, Giannone, Reichlin, 2010; Giannone, Lenza, Primiceri, 2015)
if modelOpt.dummies
    
    % Estimate BVAR
    BVAR_d = BVARestimate_dummies(data,dims,hp,modelOpt,dataOpt,miscOpt,plotOpt);
    
    % Structural identification
    BVAR_d = BVARidentification(BVAR_d,data,dims,modelOpt,dataOpt,miscOpt);
    
    % Plot IRFs
    chartName = [char(dataOpt.units),'_DummyObs'];
    BVARplots(BVAR_d.IRFsulm,chartName,dataOpt,plotOpt);

    % Display priors and save options to a text file
    printDetails(BVAR_d,'DummyObs_details',data,hp,modelOpt,dataOpt,miscOpt,plotOpt)

    % Forecast Error Variance Decomposition
    if miscOpt.FEVD
        BVAR_d.FEVD = getFEVD(BVAR_d,dims,modelOpt);

        if strcmp(plotOpt.plotType,'allShocks')
            customLegend = {};
            chartName = 'DummyObs_FEVDall';
            plotAllFEVD(BVAR_d.FEVD.FEVDmean,dataOpt,plotOpt,customLegend,chartName);
        elseif strcmp(plotOpt.plotType,'singleShock')
            chartName = 'DummyObs_FEVDchart';
            plotSingleFEVD(BVAR_d.FEVD,dataOpt,plotOpt,chartName);
        else
            error('ERROR: FEVD is available for either all or one shocks.')
        end
    end

    % Historical Decomposition
    if miscOpt.HD
        BVAR_d.HD = getHD(BVAR_d,data,dims,modelOpt);

        if strcmp(plotOpt.plotType,'allShocks')
            customLegend = {};
            chartName = 'DummyObs_HDall';
            plotAllHD(BVAR_d.HD,data,dataOpt,plotOpt,customLegend,chartName);
        elseif strcmp(plotOpt.plotType,'singleShock')
            chartName = 'DummyObs_HDchart';
            plotSingleHD(BVAR_d.HD,data,dataOpt,plotOpt,chartName);
        else
            error('ERROR: FEVD is available for either all or one shocks.')
        end
    end 

    % Structural Scenario Analysis (IRF)
    if miscOpt.SSAirf
        plotOpt.customLabels = {};
        plotOpt.estimation = 'DummyObs';
        BVAR_d.SSA = getSSA_irf(BVAR_d,data,dims,dataOpt,plotOpt,ssaOpt);
    end
    
end


%-------------------------------------------------------------------------%
% Estimating BVAR via triangularisation algorithm
% (Chan, 2019; Carriero, Clark, Marcellino, 2019)
if modelOpt.triangular

    % Estimate BVAR
    BVAR_t = BVARestimate_triangular(data,dims,hp,modelOpt,dataOpt,miscOpt,plotOpt);
    
    % Structural identification
    BVAR_t = BVARidentification(BVAR_t,data,dims,modelOpt,dataOpt,miscOpt);
    
    % Plot IRFs
    chartName = [char(dataOpt.units),'_Triangular_',miscOpt.caseTriangular];
    BVARplots(BVAR_t.IRFsulm,chartName,dataOpt,plotOpt);

    % Display priors and save options to a text file
    printDetails(BVAR_t,['Triangular_',miscOpt.caseTriangular,'_details'],data,hp,modelOpt,dataOpt,miscOpt,plotOpt)

    % Forecast Error Variance Decomposition
    if miscOpt.FEVD
        BVAR_t.FEVD = getFEVD(BVAR_t,dims,modelOpt);

        if strcmp(plotOpt.plotType,'allShocks')
            customLegend = {};
            chartName = ['Triangular_',miscOpt.caseTriangular,'_FEVDall'];
            plotAllFEVD(BVAR_t.FEVD.FEVDmean,dataOpt,plotOpt,customLegend,chartName);
        elseif strcmp(plotOpt.plotType,'singleShock')
            chartName = ['Triangular_',miscOpt.caseTriangular,'_FEVDchart'];
            plotSingleFEVD(BVAR_t.FEVD,dataOpt,plotOpt,chartName);
        else
            error('ERROR: FEVD is available for either all or one shocks.')
        end
    end

    % Historical Decomposition
    if miscOpt.HD
        BVAR_t.HD = getHD(BVAR_t,data,dims,modelOpt);

        if strcmp(plotOpt.plotType,'allShocks')
            customLegend = {};
            chartName = ['Triangular_',miscOpt.caseTriangular,'_HDall'];
            plotAllHD(BVAR_t.HD,data,dataOpt,plotOpt,customLegend,chartName);
        elseif strcmp(plotOpt.plotType,'singleShock')
            chartName = ['Triangular_',miscOpt.caseTriangular,'_HDchart'];
            plotSingleHD(BVAR_t.HD,data,dataOpt,plotOpt,chartName);
        else
            error('ERROR: FEVD is available for either all or one shocks.')
        end
    end 

    % Structural Scenario Analysis
    if miscOpt.SSAirf
        plotOpt.customLabels = {};
        plotOpt.estimation = ['Triangular_',miscOpt.caseTriangular];
        BVAR_t.SSA = getSSA_irf(BVAR_t,data,dims,dataOpt,plotOpt,ssaOpt);
    end
    
end



end %full Gibbs samplers


