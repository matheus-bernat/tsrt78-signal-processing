function [th,P,lam,epsi]=sig2ar(y,na);
%Estimate parameters in an auto-regressive model
N=length(y);
ind=(na:N-1)'*ones(1,na)-ones(N-na,1)*(0:na-1);
Phi=-reshape(y(ind),size(ind));
Y=y(na+1:N);
[th,P,lam,epsi]=sig2linmod(Y,Phi);
if nargout==0
    disp(['  a = ',num2str([th'])])
    disp(['std = ',num2str(diag(P)')])
    disp(['lam = ',num2str(lam)])
end
