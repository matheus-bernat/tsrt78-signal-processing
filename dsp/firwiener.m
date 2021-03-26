function [h,lam]=firwiener(Ryy,Rss,n);
%FIR Wiener filter
A=toeplitz(Ryy(1:n),Ryy(1:n));
b=Rss(1:n);
h=A\b;
lam=Rss(1)-h'*Rss(1:n);