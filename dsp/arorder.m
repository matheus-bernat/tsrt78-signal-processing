function [W,Uaic,Ubic]=arorder(y,nmax);
%Select AR model order 
N=length(y);
for n=1:nmax
    [th,P,lam]=sig2ar(y,n);
    W(n)=N*lam/(y'*y);
    Uaic(n)=W(n)*(1+2*n/N);
    Ubic(n)=W(n)*(1+n*log(N)/N);
end
if nargout==0
   plot(1:nmax,[W' Uaic' Ubic'],'-o')
end



