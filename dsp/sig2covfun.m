function [R,k,varR]=sig2covfun(x,kmax)
%Covariance function estimation
if nargin<2; kmax=30; end
N=length(x);
R=fastconv(x,x(end:-1:1))/N;
if nargout==3
   tmp=fastconv(R,R);
   varR=1/N*(sum(R.^2)+tmp(2*N-1:2:2*N-1+2*kmax));
end
k=0:kmax;
R=R(N:N+kmax);
if nargout==0, plot(k,R), end