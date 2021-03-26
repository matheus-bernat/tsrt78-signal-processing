function [tht,Pt]=arkf(y,na,Q)
%RLS for adaptive AR parameter estimation
if nargin<3; Q=0.01*eye(na); end
N=length(y);
th=zeros(na,1);
P=100*eye(na);
for t=na+1:N
    phi=-y(t-1:-1:t-na,:); 
    K=P*phi./(1+phi'*P*phi);
    P=P-K*phi'*P+Q;
    epsi(:,t)=y(t)-phi'*th;
    th=th+K*epsi(:,t);
    tht(:,t)=th; Pt(:,:,t)=P;
end