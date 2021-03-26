function [th,P,lam,epsi]=sig2arx(y,u,na,nb,nk)
%Estimate parameters in an ARX model
if nargin<5; nk=0; end
n=max([na nb+nk]);
N=length(y);
PhiT=[ -toeplitz(y(n:end-1),y(n:-1:n-na+1)),...
        toeplitz(u(n+1-nk:end),u(n+1-nk:-1:n-nb-nk+2)) ];
th=PhiT\y(n+1:end);
epsi=y(n+1:end)-PhiT*th;
lam=epsi'*epsi/length(epsi);
P=lam*inv(PhiT'*PhiT);
if nargout==0
    disp(['[a,b] = ',num2str([th'])])
    disp([' std  = ',num2str(diag(P)')])
    disp([' lam  = ',num2str(lam)])
end
