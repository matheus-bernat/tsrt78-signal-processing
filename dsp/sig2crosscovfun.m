function [R,k,varR]=sig2crosscovfun(x,y,kmax)
%Cross covariance function estimation
if nargin<3; kmax=30; end
N=length(x);
Rxy=fastconv(x,y(end:-1:1))/N;
if nargout==3
   Rxx=fastconv(x,x(end:-1:1))/N;
   Ryy=fastconv(y,y(end:-1:1))/N;
   tmp=fastconv(Rxy,Rxy);
   varR=1/N*(sum(Rxx.*Ryy)+tmp(2*N-1:2:2*N-1+2*kmax));
end
k=-kmax:kmax;
R=Rxy(N-kmax:N+kmax);
