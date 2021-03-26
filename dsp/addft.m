function [X,f]=addft(x,lambda,f)
%Adaptive DFT
N=length(x);
if nargin<2, lambda=0.95; end
if nargin<3, f=1/N*(0:N/2-1)'; end
W=exp(-i*2*pi*f);
X(:,1)=W*x(1);
for k=2:N
   W=exp(-i*2*pi*f*k);
   X(:,k)=lambda*X(:,k-1)+(1-lambda)*W*x(k);
end
if nargout==0
    image(f,1:N,256/max(abs(X(:)))*abs(X))
    set(gca,'Ydir','normal')
end
