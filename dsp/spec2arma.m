function [c,a,lambda]=spec2arma(num,den)
%ARMA to spectrum conversion
%  num and den polynomials in cos(wk)
%  c and a polynomials in z
n=length(den); m=length(num);
dnum=[0.5*num(end:-1:2) num(1) 0.5*num(2:end)];
r=roots(dnum);
ind=find(abs(r)<1);
c=[zeros(1,n-m) poly(r(ind))];
dden=[0.5*den(end:-1:2) den(1) 0.5*den(2:end)];
r=roots(dden);
ind=find(abs(r)<1);
a=poly(r(ind));
lambda=polyval(num,1)/polyval(den,1)/...
       (polyval(c,1)/polyval(a,1))^2
