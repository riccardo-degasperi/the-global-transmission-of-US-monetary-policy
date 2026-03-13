function r=GammaCoef(mode,sd,plotit)
% Giannone, Lenza, Primiceri 2015,
% Modified by Riccardo Degasperi, 28/11/2018

% Compute scale (theta) and shape (k) parameters of the gamma pdf
r.k=(2+mode^2/sd^2+sqrt((4+mode^2/sd^2)*mode^2/sd^2))/2;
r.theta=sqrt(sd^2/r.k);

% NOTE: knowing mode and std, one can easily derive the k and theta
% coefficients using the formulas for variance and mode. k has to be
% greater than zero, so only one solution of the resulting quadratic
% equation is valid.

% Plot hyperprior distribution
if plotit==1
xxx=[0:.000001:mode+5*sd];
plot(xxx,[xxx.^(r.k-1).*exp(-xxx./r.theta)*r.theta^-r.k/gamma(r.k)],'k--','LineWidth',2)
end