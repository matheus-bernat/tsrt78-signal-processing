function [x,f]=idft(X,T,f);
  N=length(X);
  if nargin<2; T=1; end
  if nargin<3; f=[0:N-1]'/N/T; end
  W=exp(i*2*pi*f(:)*[0:N-1]*T);
  x=1/N*W'*X(:);               
end
