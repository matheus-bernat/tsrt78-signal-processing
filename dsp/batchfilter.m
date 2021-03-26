function y=batchfilter(c,a,e,M,xi)
%Filter implementation where y is computed in batches
if nargin<5, xi=zeros(length(a)-1,1); end
if nargin<4, M=1; end
R=floor(length(e)/M); 
y=zeros(size(e));
for i=1:R
    [y(i*M-M+1:i*M), xi]=filter(c,a,e(i*M-M+1:i*M),xi);
end
y(R*M+1:end)=filter(c,a,e(i*M+1:end),xi);
