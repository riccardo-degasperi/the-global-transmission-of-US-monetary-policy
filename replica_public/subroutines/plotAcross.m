function [] = plotAcross(in,chartName,modelOpt,dataOpt,plotOpt)
% Plot responses of specific variables across units


% Unpacking

h           = modelOpt.h;

varacross   = dataOpt.varacross;
group       = dataOpt.group;
varendo     = dataOpt.varendo;
shockVar    = dataOpt.shockVar;

% Number of units
U = numel(group);

% Position of shockVar
iShock = find(strcmp(varendo,shockVar));


for ii = 1:size(varacross,2)

    % Position of variable ii
    idx = find(strcmp(varendo,varacross{ii}));

    % Get responses for variable ii across units
    for i = 1:U
        tmp = in.(group{i});
        stack(i,1) = tmp(idx,iShock);           %stack unit responses
    end

    % Initialise container of IRF median, upper and lower bounds
    IRFulm = cell(U,1);

    c1 = 0.9;   %90% confidence bands
    c2 = 0.68;  %68% confidence bands

    % Compute bands
    for z = 1:U
      for j = 1:h
      % Lower bound 90%
      IRFulm{z,1}(1,j) = quantile(stack{z,1}(j,:),(1-c1)/2);
      % Lower bound 68%
      IRFulm{z,1}(2,j) = quantile(stack{z,1}(j,:),(1-c2)/2);
      % Median
      IRFulm{z,1}(3,j) = quantile(stack{z,1}(j,:),0.5);
      % Upper bound 68%
      IRFulm{z,1}(4,j) = quantile(stack{z,1}(j,:),1-(1-c2)/2);
      % Upper bound 90%
      IRFulm{z,1}(5,j) = quantile(stack{z,1}(j,:),1-(1-c1)/2);
      end
    end

    % Display IRFs
    ylab = plotOpt.ylab(ii);
    IRFdisplay_across(IRFulm,varacross{ii},group,ylab,chartName,plotOpt)
end

