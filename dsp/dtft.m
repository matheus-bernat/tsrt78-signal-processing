function [X,f]=dtft(x,T,N)
%Approximation of the DTFT
if nargin<2, T=1; end
if nargin<3, N=1024; end
X=fft([x(:); zeros(N-length(x),1)]);
f=(0:N-1)'/N/T;