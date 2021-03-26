function [n,W]=arorderaic(y,nmax);
%Select AR model order using AIC
N=length(y);
for n=1:nmax
    [th,P,lam]=sig2ar(y,n);
    W(n)=N*lam/(y'*y)*(1+2*n/N);
end
if nargout==0
   plot(1:nmax,W,'-o')
else
   [dum,n]=min(W);
end



