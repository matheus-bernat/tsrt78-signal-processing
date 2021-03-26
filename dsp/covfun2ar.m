function [th,PN,lam]=covfun2ar(R,n);
%Estimate parameters in an auto-regressive model
Rbar=toeplitz(R(1:n),R(1:n));
fbar=-R(2:n+1);   % Rhatyy(1) to Rhatyy(n)
th=Rbar\fbar;     % Solution to Yule-Walker equations
lam=R(1);         % By definition
PN=lam*inv(Rbar); % Note, N*P, N unknown in the function