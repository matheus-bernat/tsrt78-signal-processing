function [X,f]=tdft(x,T,N)
   % DTFT approximation
   if nargin<2, T=1; end
   if nargin<3, N=2048; end
   X=fft([x(:); zeros(N-length(x),1)]);
   f=(0:N-1)'/N/T;
end
