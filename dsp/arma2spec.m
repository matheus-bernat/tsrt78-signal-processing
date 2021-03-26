function [num,den]=arma2spec(c,a,lambda)
%ARMA to spectrum conversion
%  num and den polynomials in cos(wk)
%  c and a polynomials in z
if nargin<3; lambda=1; end
m=length(c);
d=conv(c,c(end:-1:1));
num=lambda*[d(m) 2*d(m+1:2*m-1)];
n=length(a);
d=conv(a,a(end:-1:1));
den=[d(n) 2*d(n+1:2*n-1)];


