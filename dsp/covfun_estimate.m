function [R,k,varR]=covfun_estimate(x,kmax)
N=length(x);
R=conv(x,x(end:-1:1))/N;
if nargout==3
   tmp=conv(R,R);
   varR=1/N*(sum(R.^2)+tmp(2*N-1:2:2*N-1+2*kmax));
end
k=0:kmax;
R=R(N:N+kmax);
