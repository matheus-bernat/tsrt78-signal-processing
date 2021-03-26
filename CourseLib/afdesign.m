function [V,thhat]=afdesign(th,sigma,adm,adg,MC);

%AFDESIGN  Illustrates performance of adaptive AR-parameter estimation
%
%          [V,thhat] = afdesign(th,sigm,adm,adg,MC)
%
% Simulates an AR signal with parameters th structured with columns as 
% time index and rows as parameter vector. V is the parameter error and 
% thhat the parameter estimate, both averaged over MC Monte Carlo 
% realizations. The adaptive algorithm parameters adm and adg are 
% forwarded to the rarx function in System Identification toolbox.

if nargin<4,
  MC=1;
end
[n,N]=size(th);
if length(sigma)==1,
  sigma=sigma*ones(N,1);
end
thhat=zeros(size(th));
V=0;

for i=1:MC,
  randn('seed',i+1);
  e=randn(N,1);
  y=zeros(N,1);
  for t=n+1:N,
    y(t)=-y(t-1:-1:t-n)'*th(:,t)+sigma(t)*e(t);
  end
  thhatc=rarx(y,n,adm,adg);
  thhatc=thhatc';
  thhat=thhat+thhatc;
  V=V+sum(sum((th(:,41:N)-thhatc(:,41:N)).^2))/(N-40);
end
thhat=thhat/MC;
V=V/MC;
disp(['Parameter error V: ',num2str(V)])
plot(th','--'),hold on,plot(thhat','-'),hold off
title('Real parameter (dashed) and MC averaged estimate (solid)')
axis([0 N -2 2])

