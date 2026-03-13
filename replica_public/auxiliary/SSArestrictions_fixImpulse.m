function out = SSArestrictions_fixImpulse(impulse,dims,dataOpt,ssaOpt)
% SSArestrictions defines the restrictions on the path of observables and
% on the innovations/structural shocks of the VAR that deliver the chosen
% path for the observables.
% Adapted from Antolin-Diaz, Petrella, Rubio-Ramirez (2021)
%
% Please customise the 'PART TO CUSTOMISE'
%
% Last updated: 28/02/2024
% riccardo.degasperi@bancaditalia.it
%-------------------------------------------------------------------------%

% Unpack
h        = dims.h+1;           %+1 for consistency with IRFbuild function
N        = dims.N;
varendo  = dataOpt.varendo;
SSAy     = ssaOpt.SSAy;
SSAe     = ssaOpt.SSAe;

% Find positions
i_SSAy   = ismember(varendo,SSAy);
i_SSAe   = ismember(varendo,SSAe);
if ~isempty(impulse)
iImpulse = ismember(varendo,impulse);
end

%-------------------------------------------------------------------------%
% PART TO CUSTOMISE

% Impose restriction that IRF of policy variable is the same for cf and bl
fixImpulse = 1;    %1:on, 0:off
% Note: When active, fixPolicy automatically substitutes the median IRF of
% the impulse variable in yCondition(:,iImpulse), overwriting whatever is
% specified in this script.

% Choose restrictions
choose = 'default_irf';

switch choose

    case 'default_irf'

    % Conditions on observables
    yCondition = nan(h,N);                                 %nan means no restriction
    yCondition(:,i_SSAy) = zeros(h,sum(i_SSAy));           %set variables SSAy to a specific path (zero is steady state)
    Omega_f = zeros(sum(sum(isfinite(yCondition))));       %"hard" restrictions
    
    % Conditions on innovations/shocks
    eCondition = zeros(h,N);                               %zero means a shock is set to its unconditional distribution
    eCondition(:,i_SSAe) = nan;                            %SSAe shocks can deviate from their unconditional distribution (you can specify for which horizon)
    eCondition(1,iImpulse) = 1;                            %shockVar shocks set to one on impact
    %Omega_g = 1*diag(ones(sum(sum(isfinite(eCondition))),1)); %Non-driving shocks restricted to their unconditional variance
    Omega_g = zeros(sum(sum(isfinite(eCondition))));       %Non-driving shocks restricted to zero ("hard" restrictions)


    case 'test_ts'

    % Conditions on observables
    yCondition = nan(h,N);                                 %nan means no restriction
    yCondition(:,i_SSAy) = [-253.18648].*ones(h,sum(i_SSAy)); %set variables SSAy to a specific path (zero is steady state)
    Omega_f = zeros(sum(sum(isfinite(yCondition))));       %"hard" restrictions
    %Omega_f = [];                                         %"loose" restrictions
    
    % Conditions on innovations/shocks
    eCondition = zeros(h,N);                               %zero means a shock is set to its unconditional distribution
    eCondition(:,i_SSAe) = nan;                            %SSAe shocks can deviate from their unconditional distribution (you can specify for which horizon)
    Omega_g = 1*diag(ones(sum(sum(isfinite(eCondition))),1)); %Non-driving shocks restricted to their unconditional variance
    %Omega_g = zeros(sum(sum(isfinite(eCondition))));       %Non-driving shocks restricted to zero ("hard" restrictions)

    case 'test_ts2'

    % Conditions on observables
    yCondition = nan(h,N);                                 %nan means no restriction
    yCondition(:,i_SSAy) = [0.33 -253.1865].*ones(h,sum(i_SSAy)); %set variables SSAy to a specific path (zero is steady state)
    Omega_f = zeros(sum(sum(isfinite(yCondition))));       %"hard" restrictions
    %Omega_f = [];                                         %"loose" restrictions
    
    % Conditions on innovations/shocks
    eCondition = zeros(h,N);                               %zero means a shock is set to its unconditional distribution
    eCondition(:,i_SSAe) = nan;                            %SSAe shocks can deviate from their unconditional distribution (you can specify for which horizon)
    Omega_g = 1*diag(ones(sum(sum(isfinite(eCondition))),1)); %Non-driving shocks restricted to their unconditional variance
    %Omega_g = zeros(sum(sum(isfinite(eCondition))));       %Non-driving shocks restricted to zero ("hard" restrictions)








    case 'default_cf'

    % Conditions on observables
    yCondition = nan(h,N);
    yCondition(:,i_SSAy) = zeros(h,sum(i_SSAy));
    Omega_f = zeros(sum(sum(isfinite(yCondition))));

    % Conditions on innovations/shocks
    eCondition = zeros(h,N);
    eCondition(1,iImpulse) = 1;
    Omega_g = zeros(sum(sum(isfinite(eCondition))));

    case 'default_cf2'

    % Conditions on observables
    yCondition = nan(h,N);
    yCondition(:,i_SSAy) = zeros(h,sum(i_SSAy));
    Omega_f = zeros(sum(sum(isfinite(yCondition))));

    % Conditions on innovations/shocks
    eCondition = nan(h,N);
    eCondition(1,iImpulse) = 1;
    Omega_g = zeros(sum(sum(isfinite(eCondition))));

    case 'default_cf3'

    % Conditions on observables
    yCondition = nan(h,N);
    yCondition(:,i_SSAy) = zeros(h,sum(i_SSAy));
    Omega_f = zeros(sum(sum(isfinite(yCondition))));

    % Conditions on innovations/shocks
    eCondition = nan(h,N);
    eCondition(1,:) = 0;
    eCondition(1,iImpulse) = 1;
    Omega_g = zeros(sum(sum(isfinite(eCondition))));

    case 'default_cf4'

    % Conditions on observables
    yCondition = nan(h,N);
    yCondition(:,i_SSAy) = zeros(h,sum(i_SSAy));
    Omega_f = zeros(sum(sum(isfinite(yCondition))));

    % Conditions on innovations/shocks
    eCondition = zeros(h,N);
    eCondition(:,iImpulse) = nan;
    eCondition(1,iImpulse) = 1;
    Omega_g = zeros(sum(sum(isfinite(eCondition))));


    case 'demo1_houthi'

    % Conditions on observables
    yCondition = nan(h,N);
    yCondition(:,i_SSAy) = 20*ones(h,1);
    Omega_f = zeros(sum(sum(isfinite(yCondition))));

    % Conditions on innovations/shocks
    eCondition = zeros(h,N);
    eCondition(:,i_SSAe) = nan;
    eCondition(1,iImpulse) = 1;
    Omega_g = zeros(sum(sum(isfinite(eCondition))));

    case 'demo2_houthi'

    % Conditions on observables
    yCondition = nan(h,N);
    yCondition(:,i_SSAy) = [20*ones(8,1); zeros(h-8,1)];
    Omega_f = zeros(sum(sum(isfinite(yCondition))));

    % Conditions on innovations/shocks
    eCondition = zeros(h,N);
    eCondition(:,i_SSAe) = nan;
    eCondition(1,iImpulse) = 1;
    Omega_g = zeros(sum(sum(isfinite(eCondition))));

    case 'demo3_houthi'  %conditional forecast

    % Conditions on observables
    yCondition = nan(h,N);
    yCondition(:,i_SSAy) = [20*ones(8,sum(i_SSAy)); zeros(h-8,sum(i_SSAy))];
    Omega_f = zeros(sum(sum(isfinite(yCondition))));

    % Conditions on innovations/shocks
    eCondition = zeros(h,N);
    eCondition(1,iImpulse) = 1;
    Omega_g = zeros(sum(sum(isfinite(eCondition))));
            

end
%-------------------------------------------------------------------------%

% Checks
if fixImpulse
   check = ismember(SSAy,impulse);
   if check == 0
       error('ERROR: include "impulse" among the restricted variables in ssaOpt.SSAy if the option fixImpulse is active.')
   end
end

% Pack output
out.yCondition = yCondition;
out.Omega_f    = Omega_f;
out.eCondition = eCondition;
out.Omega_g    = Omega_g;
out.chooseSS   = choose;
out.fixImpulse = fixImpulse;