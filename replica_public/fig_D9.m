%-------------------------------------------------------------------------%
% The Global Transmission of U.S. Monetary Policy                         %
% Riccardo Degasperi, Simon Hong, Giovanni Ricco                          %
%                                                                         %
% Replica of Figure D.9                                                   %
%                                                                         %
% URL to the file:                                                        %
% Fig.D.9 :  ./plots/fig_D9/AE_oil_proxy_comparison.pdf                   %
%                                                                         %
% last updated: 03/03/2026                                                %
%-------------------------------------------------------------------------%

clear
close all
clc

addpath(genpath(pwd))  %path to child functions

% Load options
base   = load('AEs/z_oilout');
MARext = load('AEs_MARext/z_oilout');
JK     = load('AEs_JK/z_oilout');
BS     = load('AEs_BS/z_oilout');

plotOpt  = base.plotOpt;
dataOpt  = base.dataOpt;
modelOpt = base.modelOpt;
miscOpt  = base.miscOpt;


% Get varplot indices
if sum(ismember(dataOpt.shockVar,dataOpt.varplot)) < numel(dataOpt.shockVar)
    iS = ismember(dataOpt.shockVar,dataOpt.varplot)==0;
    dataOpt.varplot = [dataOpt.varplot,dataOpt.shockVar(iS)];
end
idx = ismember(dataOpt.varendo,dataOpt.varplot);


% Pack IRF for plot
stoIRF.baseline = base.IRFmedian(idx,idx);
stoIRF.MARext   = MARext.IRFmedian(idx,idx);
stoIRF.JK       = JK.IRFmedian(idx,idx);
stoIRF.BS       = BS.IRFmedian(idx,idx);

% Plot figure
plotOpt.folder = 'fig_D9/';
legNames = {'Baseline','MAR extended','Jarocinski-Karadi','Bauer-Swanson'};
chartName = 'AE_oil_proxy_comparison';
longnames = dataOpt.varendolong(idx);
style = {'k-','k--','k-.','k:','ko'};
IRFdisplay_structure2(stoIRF,dataOpt.varplot,longnames,dataOpt.shockVar,legNames,style,'varplot',plotOpt,chartName)