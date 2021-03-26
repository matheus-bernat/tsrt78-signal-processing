function [Phi,f]=sig2blackmantukey(x,gamma,T)
%Spectral analysis using Blackman-Tukey
if nargin<3, T=1; end
gamma=round(pi*gamma);
R=sig2covfun(x,gamma);
R=[R(end:-1:2);R];
w=getwindow(2*gamma+1,'hamming');
R=w.*R;
Phi=T*abs(fft(R));
Phi=Phi(1:gamma);
f=(0:gamma-1)/gamma/T;
