%-------------------------------------------------------------------------%
% The Global Transmission of U.S. Monetary Policy                         %
% Riccardo Degasperi, Simon Hong, Giovanni Ricco                          %
%                                                                         %
% Replica of all figures of the paper                                     %
%                                                                         %
% last updated: 03/03/2026                                                %
%-------------------------------------------------------------------------%


tic

% Produce all figures in the main text
run('fig_1_C1.m')            % Global and domestic effects of US monetary policy
run('fig_2_C4_C5.m')         % Forecast Error Variance Decomposition
run('fig_3_4_C2_C3_E1_H2.m') % Median responses of Advanced and Emerging Economies
run('fig_5_C6_C7.m')         % The commodity price channel
run('fig_6a_C8.m')           % Disentangling the channels -- Global
run('fig_6b_C9.m')           % Disentangling the channels -- Advanced
run('fig_6c_C10.m')          % Disentangling the channels -- Emerging

% Produce all remaining figures in the appendix

run('fig_D1.m')
run('fig_D2.m')
run('fig_D3.m')
run('fig_D4a.m')
run('fig_D4b.m')
run('fig_D5.m')
run('fig_D6a.m')
run('fig_D6b.m')
run('fig_D7.m')
run('fig_D8a.m')
run('fig_D8b.m')
run('fig_D9.m')
run('fig_D10.m')
run('fig_D11.m')
run('fig_D12.m')
run('fig_D13_D14_D15.m')
run('fig_D16.m')
run('fig_D17.m')
run('fig_D18.m')
run('fig_D19.m')

toc

% Tested with Matlab R2024a
