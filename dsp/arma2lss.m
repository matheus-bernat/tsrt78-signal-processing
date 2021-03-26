function [A,B,C,D]=arma2lss(c,a)
%State space to ARMA conversion
n=length(a);
A = [-a(2:end); eye(n-2,n-1)];
B = eye(n-1,1);
C = c(2:n) - c(1) * a(2:end);
D = c(1);
if nargout==1 % Struct form of model
    m.A=A; m.B=B; m.C=C; m.D=D; A=m;
end
