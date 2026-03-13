function [Yd,Xd] = dummyObs(Ybar,Xbar,hp,dims)

% Retrieve hyperpriors
niw   = hp.niw;
soc   = hp.soc;
cop   = hp.cop;
covid = hp.covid;
l1    = hp.l1;
l2    = hp.l2;
l3    = hp.l3;
l4    = hp.l4;
l5    = hp.l5;
l7    = hp.l7; %pandemic prior
iRW   = hp.iRW;

% Retrieve dimensions
m  = dims.m;
mc = dims.mc;
p  = dims.p;
N  = dims.N;

% Use diagonal elements of the scale matrix of the IW prior on the residual variance
SS = hp.SS;

% Dummy observations for NIW prior ---------------------------------------%
% (Normal Inverse Wishart)
if niw

    Yrw  = [diag(SS.*iRW)/l1;
            zeros(N*(p-1),N);
            diag(SS)        ;
            zeros(m,N)     ];

    Xrw  = [kron(diag((1:p).^l2),diag(SS)/l1) zeros(N*p,m) ;
            zeros(N,N*p)            zeros(N,m)             ;
            zeros(m,N*p)            kron(1/(l1*l5),eye(m))];
      
% NOTE: relative position of second and third block does not matter.
        
else
    Yrw = [];  Xrw = [];
end

% Dummy observations for SoC prior ---------------------------------------%
% (Sum of Coefficients, or No Cointegration)
if soc

    Ysoc = diag(Ybar.*iRW)/l3;

    Xsoc = [kron(ones(1,p),Ysoc) zeros(N,m)];

else
    Ysoc = [];  Xsoc = [];
end

% Dummy observations for CoP prior ---------------------------------------%
% (Co-persistence, or Dummy Initial Observation)
if cop

    Ycop = Ybar'/l4;

    Xcop = [kron(ones(1,p),Ycop) 1/l4 Xbar/l4];

else
    Ycop = [];  Xcop = [];
end

% Combine dummy observations
Yd = [Yrw; Ysoc; Ycop];
Xd = [Xrw; Xsoc; Xcop];


% Dummy observations for pandemic prior ----------------------------------%
% (Cascaldi-Garcia, 2024)
if covid

    Ycovid = zeros(mc,N);

    Xcovid = [zeros(mc,N*p) zeros(mc,1) l7.*eye(mc)];

    Yd = [Yd;Ycovid];
    Xd = [Xd zeros(size(Xd,1),mc);Xcovid];

end

