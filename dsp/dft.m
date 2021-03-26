function [X,f]=dft(x,T,f);
  N=length(x);
  if nargin<2; T=1; end
  if nargin<3; f=[0:N-1]'/N/T; end
  W=exp(i*2*pi*f(:)*[0:N-1]*T);  
  X=W*x(:);                      
end
