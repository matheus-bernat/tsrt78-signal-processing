function [c,a]=lss2arma(A,B,C,D)
%ARMA to state space conversion
if nargin==1; m=A; A=m.A; B=m.B; C=m.C; D=m.D; end
a=real(poly(A));
c=real(poly(A-B*C)) + (D-1)*a;

