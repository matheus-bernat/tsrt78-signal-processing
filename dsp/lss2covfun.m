function [R,k]=lss2covfun(m,kmax)
%State space model to covariance function conversion
%if ~islss(m), error('m must be a state space model'), end
k=0:kmax;
R=zeros(size(k));
Pibar=dlyap(m.A,m.B*m.Q*m.B');
R(1)=m.C*Pibar*m.C'+m.R; 
tmp=m.C;
S=zeros(size(m.A,1),size(m.C,1));
if isfield(m,'S'), S=m.S; end
for l=1:kmax
   R(l+1)=tmp*m.A*Pibar*m.C'+tmp*m.B*S;
   tmp=tmp*m.A;
end
if nargout==0, plot(k,R), end