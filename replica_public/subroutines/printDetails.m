function [] = printDetails(in,nameFile,data,hp,modelOpt,dataOpt,miscOpt,plotOpt)
% Prints a text file containing some relevant estimation details.

% INPUT:
% - in       : output from estimateBVAR or estimateBLP;
% - nameFile : either 'BVAR_details' or 'BLP_details';
% - data     : output from makeXY;
% - hp       : hyperparameters and related options;
% - modelOpt : structure containing model options;
% - dataOpt  : structure containing data options; 
% - miscOpt  : structure containing misc options;

% OUTPUT:
% - text file named [units,nameFile].txt in folder "plots".

% r.degasperi@warwick.ac.uk
%-------------------------------------------------------------------------%

% Determine if BLP or BVAR
if strcmp('BLP_details',nameFile)
    triangular = 0;
elseif strcmp('DummyObs_details',nameFile)
    triangular = 0;
elseif contains(nameFile,'riangular')
    triangular = 1;
elseif strcmp(nameFile,'TVAR_details')
    triangular = 0;
else
    error('ERROR: nameFile can either be "BVAR_details", "BLP_details", "TVAR_details" or "Triangular_[algorithm]".')
end


%-------------------------------------------------------------------------%
% Unpacking
dates              = data.dates;
iRW_endo           = data.iRW_endo;
iRW_global         = data.iRW_global;
iLogsEn            = data.iLogsEn;
iLogsEn_global     = data.iLogsEn_global;
iLogsEx            = data.iLogsEx;
p                  = modelOpt.p;
h                  = modelOpt.h;
identification     = modelOpt.identification;
sim                = modelOpt.sim;
burnin             = modelOpt.burnin;
jump               = modelOpt.jump;
shutVar            = modelOpt.shutVar;
shutEq             = modelOpt.shutEq;
caseTriangular     = miscOpt.caseTriangular;
varendo            = dataOpt.varendo;
varexo             = dataOpt.varexo;
units              = dataOpt.units;
shockVar           = dataOpt.shockVar;
shockUnit          = dataOpt.shockUnit;
shockSign          = dataOpt.shockSign;
variv              = dataOpt.variv;
sto_varendolong    = dataOpt.sto_varendolong;
varendolong_global = dataOpt.varendolong_global;
varexolong         = dataOpt.varexolong;
dataset            = dataOpt.dataset;
folder             = plotOpt.folder;
hpVAR              = in.hpVAR;

if strcmp(identification,'iv')
    Fstat          = in.Fstat;
    Fdistr         = in.Fdistr;
    Frobust        = in.Frobust;
    Frdistr        = in.Frdistr;
    RMdistr         = in.RMdistr;
end

neig           = in.neig;



% Patch for chartNames when units are more than 1
if numel(units) > 1
    units0 = join(string(units),'_');
else
    units0 = units;
end


%-------------------------------------------------------------------------%
% Display on-screen details

% Display estimation sample
disp(' ')
disp(['Common sample from ',char(dates(1)),' to ',char(dates(end))])
disp(' ')

% Display priors
disp('--------- Optimal priors for VAR ---------')
disp(['NIW tightness          : ',num2str(hpVAR.l1)])
if triangular
disp(['Cross-var shrinkage    : ',num2str(hpVAR.l6)])
end
if hpVAR.lag; tmp = num2str(hpVAR.l2); else tmp = 'off'; end
disp(['Lag decay              : ',tmp])
if hpVAR.soc; tmp = num2str(hpVAR.l3); else tmp = 'off'; end
disp(['Sum-of-coefficients    : ',tmp])
if hpVAR.cop; tmp = num2str(hpVAR.l4); else tmp = 'off'; end
disp(['Co-persistence         : ',tmp])
disp(['Constant and exogenous : ',num2str(hpVAR.l5)])
if hp.covid
disp(['Pandemic priors        : ',num2str(hpVAR.l7)])
end
if hp.GLP ~= 0 && ~triangular
disp(['Log-posterior          : ',num2str(hpVAR.fh)])
end
disp('------------------------------------------')


% Display instrument relevance
if strcmp(identification,'iv')
    disp(' ')
    disp('-------- Relevance of instruments --------')
    for kk = 1:size(RMdistr,1)
    disp(['Shock ',num2str(kk),':'])
    disp(['First-stage F stat     : ',num2str(Fstat(kk)),' [',num2str(Fdistr(kk,1)),',',num2str(Fdistr(kk,2)),',',num2str(Fdistr(kk,3)),']'])
    disp(['Robust F stat          : ',num2str(Frobust(kk)),' [',num2str(Frdistr(kk,1)),',',num2str(Frdistr(kk,2)),',',num2str(Frdistr(kk,3)),']'])
    disp(['Reliability (eig ',num2str(kk),')    : ',num2str(RMdistr(kk,2)),' [',num2str(RMdistr(kk,1)),',',num2str(RMdistr(kk,2)),',',num2str(RMdistr(kk,3)),']'])
    disp('------------------------------------------')
    if size(RMdistr,1) > 1
    disp('Note: shocks numbered according to their position in varendo.')
    end
    end
end
disp(' ')
disp(' ')

% Define filename
fileName = [char(units0),'_',nameFile,'.txt'];

% Save model details
if ~isempty(folder)
    folder_ = ['./plots/',folder];
    if not(isfolder(folder_))
        mkdir(folder_)
    end
else
    folder_ = './plots';
end
    
oldcd = cd(folder_);                 %change directory to plots
fileID = fopen(fileName,'w');        %open text file
cd(oldcd)                            %restore original directory


%-------------------------------------------------------------------------%
% Print details

% Title & datetime
fprintf(fileID,'ESTIMATION DETAILS:\r\n');  %overall title
fprintf(fileID,'Date:');
fprintf(fileID,' %s',string(datetime));     %print date
fprintf(fileID,'\r\n\r\n\r\n');

% Relevance diagnostics for Proxy-SVAR
if strcmp(identification,'iv')
    fprintf(fileID,'RELEVANCE DIAGNOSTICS FOR PROXY-SVAR:\r\n\r\n');
    for kk = 1:size(RMdistr,1)
    t1 = ['Shock ',num2str(kk),':'];
    t2 = ['First-stage F stat     : ',num2str(Fstat(kk)),' [',num2str(Fdistr(kk,1)),',',num2str(Fdistr(kk,2)),',',num2str(Fdistr(kk,3)),']'];
    t3 = ['Robust F stat          : ',num2str(Frobust(kk)),' [',num2str(Frdistr(kk,1)),',',num2str(Frdistr(kk,2)),',',num2str(Frdistr(kk,3)),']'];
    t4 = ['Reliability (eig ',num2str(kk),')    : ',num2str(RMdistr(kk,2)),' [',num2str(RMdistr(kk,1)),',',num2str(RMdistr(kk,2)),',',num2str(RMdistr(kk,3)),']'];
    t5 = ('------------------------------------------');
    fprintf(fileID,'%s \r\n',t1,t2,t3,t4,t5);
    end
    fprintf(fileID,'\r\n\r\n');
end


% Number of explosive draws
fprintf(fileID,'# explosive draws for VAR: %s\r\n',num2str(neig(1)));
fprintf(fileID,'\r\n\r\n');

% Priors
t1 = ['Normal-Inverse-Wishart : ',num2str(hpVAR.l1)];
t2 = ['Lag decay              : ',num2str(hpVAR.l2)];
t3 = ['Sum-of-coefficients    : ',num2str(hpVAR.l3)];
t4 = ['Co-persistence         : ',num2str(hpVAR.l4)];
t5 = ['Constant and exogenous : ',num2str(hpVAR.l5)];
t6 = ['Cross-variable         : ',num2str(hpVAR.l6)];
t7 = ['Panemic priors         : ',num2str(hpVAR.l7)];

fprintf(fileID,'PRIORS FOR VAR:\r\n\r\n');
if triangular
    fprintf(fileID,'%s \r\n',t1,t6,t2,t3,t4,t5);
else
    fprintf(fileID,'%s \r\n',t1,t2,t3,t4,t5,t7);
end
fprintf(fileID,'\r\n\r\n');

% Priors for triangular
if triangular
    fprintf(fileID,'DETAILS OF ASYMMETRIC PRIORS:\r\n\r\n');
    fprintf(fileID,'Shut equations         :');
    fprintf(fileID,' %s',string(shutEq));
    fprintf(fileID,'\r\n');
    fprintf(fileID,'Shut coefficients      :');
    fprintf(fileID,' %s',string(shutVar));
    fprintf(fileID,'\r\n');
    fprintf(fileID,'Algorithm              :');
    fprintf(fileID,' %s',string(caseTriangular));
    fprintf(fileID,'\r\n\r\n\r\n');
end

% Prior optimisation
t1 = ['Prior selection        : ',num2str(hp.GLP)];
t2 = ['Hyperpriors            : ',num2str(hp.hyperpriors)];
t3 = ['NIW on/off             : ',num2str(hp.niw)];
t4 = ['Lag decay on/off       : ',num2str(hp.lag)];
t5 = ['SoC on/off             : ',num2str(hp.soc)];
t6 = ['Cop on/off             : ',num2str(hp.cop)];
t7 = ['Std on/off             : ',num2str(hp.std)];
t8 = ['Cross on/off           : ',num2str(hp.cross)];


fprintf(fileID,'PRIOR OPTIMISATION:\r\n\r\n');
fprintf(fileID,'%s \r\n',t1,t2);
fprintf(fileID,'\r\n');
if triangular
fprintf(fileID,'%s \r\n',t3,t8,t4,t5,t6,t7);  
else
fprintf(fileID,'%s \r\n',t3,t4,t5,t6,t7);
end
fprintf(fileID,'\r\n\r\n');


% Model options
t1  = ['Lags VAR               : ',num2str(p)];

t3  = ['Horizons               : ',num2str(h)];
t4  = ['Identification         : ',identification];
t6  = 'Sign restricitons file : N/A';
t8  = ['Sampler iterations     : ',num2str(sim)];
t9  = ['Burnin                 : ',num2str(burnin)];
t10 = ['Jump                   : ',num2str(jump)];

fprintf(fileID,'MODEL OPTIONS:\r\n\r\n');
fprintf(fileID,'%s \r\n',t1,t3,t4,t6,t8,t9,t10);
fprintf(fileID,'\r\n\r\n');

% Data options
fprintf(fileID,'DATA OPTIONS:\r\n\r\n');
fprintf(fileID,'Dataset                :');
fprintf(fileID,' %s',dataset);
fprintf(fileID,'\r\n');
fprintf(fileID,'Initial period         :');
fprintf(fileID,' %s',char(dates(1)));
fprintf(fileID,'\r\n');
fprintf(fileID,'Final period           :');
fprintf(fileID,' %s',char(dates(end)));
fprintf(fileID,'\r\n');
fprintf(fileID,'Shock variable         :');
fprintf(fileID,' %s',string(shockVar));
fprintf(fileID,'\r\n');
if ~isempty(shockUnit)
    fprintf(fileID,'Shock unit             :');
    fprintf(fileID,' %s',char(shockUnit));
    fprintf(fileID,'\r\n');
end
fprintf(fileID,'Sign of the shock      :');
fprintf(fileID,' %s',shockSign);
fprintf(fileID,'\r\n');
fprintf(fileID,'External instrument    :');
fprintf(fileID,' %s',string(variv));
fprintf(fileID,'\r\n');

% Variables
fprintf(fileID,'Endogenous variables   :');
fprintf(fileID,' %s',string(varendo));
fprintf(fileID,'\r\n');
fprintf(fileID,'Exogenous variables    :');
fprintf(fileID,' %s',string(varexo));
fprintf(fileID,'\r\n');
fprintf(fileID,'Units                  :');
fprintf(fileID,' %s',string(units));
fprintf(fileID,'\r\n');

% Transformations
fprintf(fileID,'Transformations  : \r\n');
fprintf(fileID,' %35s %s %s \r\n','Long names','Logs','RW');
p1 = [sto_varendolong,varendolong_global,varexolong]';
p2 = string([iLogsEn;iLogsEn_global;iLogsEx]);
p3 = repmat('n/a',size(varexo,2),1);
p4 = [string([iRW_endo;iRW_global]);p3];
t1 = [p1,p2,p4]';
fprintf(fileID,' %35s %4s %2s\r\n',t1);

fclose(fileID);
