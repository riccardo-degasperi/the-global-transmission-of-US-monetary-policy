%-------------------------------------------------------------------------%
% The Global Transmission of U.S. Monetary Policy                         %
% Riccardo Degasperi, Simon Hong, Giovanni Ricco                          %
%                                                                         %
% Replica of Figure D.4.a                                                 %
%                                                                         %
% URL to the file:                                                        %
% Fig.D.4.a :  ./plots/fig_D4/overlay/AE_MARextended_comparison.pdf       %
%                                                                         %
% last updated: 03/03/2026                                                %
%-------------------------------------------------------------------------%

clear
close all
clc

addpath(genpath(pwd))  %path to child functions

% Get results for advanced economies
run('runAE_MARext.m')

%% Load options
AEmodelOpt = load('sto_AEs/modelOpt');
AEdataOpt  = load('sto_AEs/dataOpt');
AEplotOpt  = load('sto_AEs/plotOpt');
AEmiscOpt  = load('sto_AEs/miscOpt');

EMmodelOpt = load('sto_AEs_MARext/modelOpt');
EMdataOpt  = load('sto_AEs_MARext/dataOpt');
EMplotOpt  = load('sto_AEs_MARext/plotOpt');
EMmiscOpt  = load('sto_AEs_MARext/miscOpt');

% Load IRFs
load('sto_AEs/sto_groupBaseline_BVAR_t');
AEIRF = IRFmedian;
load('sto_AEs_MARext/sto_groupBaseline_BVAR_t');
EMIRF = IRFmedian;

% Get varplot indices
if sum(ismember(AEdataOpt.shockVar,AEdataOpt.varplot)) < numel(AEdataOpt.shockVar)
    iS = ismember(AEdataOpt.shockVar,AEdataOpt.varplot)==0;
    AEdataOpt.varplot = [AEdataOpt.varplot,AEdataOpt.shockVar(iS)];
end
iAE = ismember(AEdataOpt.varendo,AEdataOpt.varplot);
if sum(ismember(EMdataOpt.shockVar,EMdataOpt.varplot)) < numel(EMdataOpt.shockVar)
    iS = ismember(EMdataOpt.shockVar,EMdataOpt.varplot)==0;
    EMdataOpt.varplot = [EMdataOpt.varplot,EMdataOpt.shockVar(iS)];
end
iEM = ismember(EMdataOpt.varendo,EMdataOpt.varplot);

% Index useful to skip plot of variables available only for AE
iskip = ismember(AEdataOpt.varplot,EMdataOpt.varplot);

% Pack IRF for plot
stoIRF.EMIRF = EMIRF(iAE,iAE);
stoIRF.AEIRF = AEIRF(iAE,iAE);

% Plot figure
AEplotOpt.folder = 'fig_D4/overlay/';
legNames = {'MAR extended','Baseline'};
chartName = 'AE_MARextended_comparison';
longnames = AEdataOpt.varendolong(iAE);
plotOpt.ylab = AEdataOpt.ylab(iAE);
IRFdisplay_overlay(stoIRF,AEdataOpt.varplot,iskip,longnames,AEdataOpt.shockVar,legNames,'varplot',AEplotOpt,chartName)

