%-------------------------------------------------------------------------%
% The Global Transmission of U.S. Monetary Policy                         %
% Riccardo Degasperi, Simon Hong, Giovanni Ricco                          %
%                                                                         %
% Replica of Figures 3, 4, C.2, C.3, E.1, H.2                             %
%                                                                         %
% URL to the file:                                                        %
% Fig.3     : ./plots/fig_3_4_C2_C3_E1_H2/overlay/overlay_AE_EM_1.pdf     %
% Fig.4     : ./plots/fig_3_4_C2_C3_E1_H2/AEs/AA_commodity_median-median_varplot_1.pdf
% Fig.C.2.a : ./plots/fig_3_4_C2_C3_E1_H2/AEs/AA_Triangular_median-median_varplot_1.pdf
% Fig.C.2.b : ./plots/fig_3_4_C2_C3_E1_H2/EMs/AA_Triangular_median-median_varplot_1.pdf
% Fig.C.3   : ./plots/fig_3_4_C2_C3_E1_H2/EMs/AA_commodity_median-median_varplot_1.pdf
% Fig.E.1.a : ./plots/fig_3_4_C2_C3_E1_H2/AEs/AA_Triangular_median-median_homogeneity_1.pdf
% Fig.E.1.b : ./plots/fig_3_4_C2_C3_E1_H2/EMs/AA_Triangular_median-median_homogeneity_1.pdf
% Fig.H.2   : ./plots/fig_3_4_C2_C3_E1_H2/EMs/AA_Triangular_median-median_homogeneity_1.pdf
%                                                                         %
% last updated: 02/03/2026                                                %
%-------------------------------------------------------------------------%

clear
close all
clc

addpath(genpath(pwd))  %path to child functions

% Get results for advanced economies
run('runAE.m')

% Get results for emerging economies
run('runEM.m')


%% Load options
AEmodelOpt = load('sto_AEs/modelOpt');
AEdataOpt  = load('sto_AEs/dataOpt');
AEplotOpt  = load('sto_AEs/plotOpt');
AEmiscOpt  = load('sto_AEs/miscOpt');

EMmodelOpt = load('sto_EMs/modelOpt');
EMdataOpt  = load('sto_EMs/dataOpt');
EMplotOpt  = load('sto_EMs/plotOpt');
EMmiscOpt  = load('sto_EMs/miscOpt');

% Load IRFs
load('sto_AEs/sto_groupBaseline_BVAR_t');
AEIRF = IRFmedian;
load('sto_EMs/sto_groupBaseline_BVAR_t');
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
stoIRF.AEIRF = AEIRF(iAE,iAE);
stoIRF.EMIRF = EMIRF(iEM,iEM);

% Plot figure
AEplotOpt.folder = 'fig_3_4_C2_C3_E1_H2/overlay/';
legNames = {'AE','EM'};
chartName = 'overlay_AE_EM';
longnames = AEdataOpt.varendolong(iAE);
plotOpt.ylab = AEdataOpt.ylab(iAE);
IRFdisplay_overlay_AE_EM(stoIRF,AEdataOpt.varplot,iskip,longnames,AEdataOpt.shockVar,legNames,'varplot',AEplotOpt,chartName)

