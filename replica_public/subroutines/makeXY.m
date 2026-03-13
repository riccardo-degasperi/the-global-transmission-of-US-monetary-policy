function [data,dims,hp,dataOpt] = makeXY(data,hp,modelOpt,dataOpt)
% makeXY prepares the model matrices and performs some preliminary
% housekeeping operations.
%
% last modified: 26/03/2024
% riccardo.degasperi@bancaditalia.it
%-------------------------------------------------------------------------%
disp('    Prepare data matrices')

% Unpacking
p                     = modelOpt.p;
pexo                  = modelOpt.pexo;
h                     = modelOpt.h;
constant              = modelOpt.constant;
endo                  = data.endo;
endo_global           = data.endo_global;
exo                   = data.exo;
iRW_endo              = data.iRW_endo;
iRW_global            = data.iRW_global;
iLogsEn               = data.iLogsEn;
iLogsEn_global        = data.iLogsEn_global;
shockVar              = dataOpt.shockVar;
varendo               = dataOpt.varendo;
varendo_global        = dataOpt.varendo_global;
varplot               = dataOpt.varplot;
varendolong           = dataOpt.varendolong;
varendolong_global    = dataOpt.varendolong_global;
ylab_endo             = dataOpt.ylab_endo;
ylab_global           = dataOpt.ylab_global;
covid                 = hp.covid;

%=========================================================================%
% Generate pandemic time dummies (Cascaldi-Garcia, 2024)
if covid

    dates = data.dates;
    lTc = dataOpt.lT_covid;
    uTc = dataOpt.uT_covid;

    ilT = find(strcmp(dates,lTc));
    iuT = find(strcmp(dates,uTc));

    if isempty(ilT) || isempty(iuT)
        error('ERROR: pandemic time dummies misspecified. Control settings.')
    end

    Xcovid = zeros(size(dates,1),iuT-ilT+1);
    Xcovid(ilT:iuT,:) = eye(iuT-ilT+1);

else

    Xcovid = [];

end
%=========================================================================%

% Concatenate endo2 (global endogenous variables)
endo = [endo endo_global];

% Generate matrix Y
Y = endo(p+1:end,:);                                %(T-p)xN matrix of values at time t

% Generate matrix X
if constant == 0
    X = lagX(endo,p);                               %(T-p)x(NUp) matrix of lagged values
    Xex = [exo(pexo+1:end,:) lagX(exo,pexo)];
    X = [X Xex(p-pexo+1:end,:) Xcovid(p+1:end,:)];  %(T-p)x(NUp+m) add exogenous
elseif constant == 1
    X = lagX(endo,p);
    Xex = [exo(pexo+1:end,:) lagX(exo,pexo)];
    X = [X ones(size(X,1),1) Xex(p-pexo+1:end,:) Xcovid(p+1:end,:)];
elseif constant == 2
    X = lagX(endo,p);
    Xex = [exo(pexo+1:end,:) lagX(exo,pexo)];
    trend = [1:size(X,1)]';
    X = [X ones(size(X,1),1) trend Xex(p-pexo+1:end,:) Xcovid(p+1:end,:)];
end

% Determine dimensions of the new system
[T,N] = size(endo);                        %this is the reshaped endo
if constant == 0
    m = size(Xex,2);                       %no constant
elseif constant == 1
    m = size(Xex,2) + 1;                   %plus 1 for the constant
elseif constant == 2
    m = size(Xex,2) + 2;                   %plus 2 for constant & trend
end
[~,Nx] = size(exo);
mc = size(Xcovid,2);

% If global endogenous are there, modify iRW accordingly
iRW = [iRW_endo; iRW_global];

% If global endogenous are there, modify iLogs accordingly
iLogs = [iLogsEn; iLogsEn_global];

% If global endogenous are there, modify ylab accordingly
ylab = [ylab_endo; ylab_global];

% This is needed when printing the results
sto_varendolong = varendolong;

% If global endogenous are there, update varendo and varendolong
varendo     = [varendo varendo_global];
varendolong = [varendolong varendolong_global];


% Pack output
dims.T               = T;
dims.N               = N;
dims.Nx              = Nx;
dims.m               = m;
dims.mc              = mc;
dims.p               = p;
dims.pexo            = pexo;
dims.h               = h;

data.endo            = endo;
data.X               = X;
data.Y               = Y;
data.Xex             = Xex;
data.Xcovid          = Xcovid;
data.iLogs           = iLogs;           %update index for RW prior

hp.iRW               = iRW;             %update index for RW prior

dataOpt.varplot      = varplot;         %update names
dataOpt.varendo      = varendo;
dataOpt.varendolong  = varendolong;
dataOpt.sto_varendolong = sto_varendolong;
dataOpt.shockVar     = shockVar;
dataOpt.ylab         = ylab;




