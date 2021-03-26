function [th,P,lam]=sig2arpp(y,na);
%Estimate parameters in an auto-regressive model
R=sig2covfun(y,na+1);
[th,PN,lam]=covfun2ar(R,na);
P=PN/length(y);


