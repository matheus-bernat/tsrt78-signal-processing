function [th,P,lam,epsi,yHat]=sig2linmod(y,Phi);
% [thHat,P,lam,epsi] = sig2linmod(y,phi);
% Inputs
%        - y is the whole output signal to be estimated
%        - phi is: y(t) = theta * phi(t)' + e(t)
% Outpputs
%        - th = estimated model parameters. The order is decided by the
%        size of phi
%        - P = covariance of model parameters
%        - lam = ?
%        - epsi = the accumulated error between the real data y and the
%        estimated signal yHat
% See page 226
% Estimate parameters in a linear model
th=Phi\y;
epsi=y-Phi*th;
R=Phi'*Phi;
lam=epsi'*epsi/length(y);
P=lam*inv(R);

% Matheus added:
yHat = Phi*th;
