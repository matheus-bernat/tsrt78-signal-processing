function mpred=lss2predictor(m,kstep)
%k-step predictor for state space model on innovation form
if any(m.Q~=m.R) | any(m.Q~=m.S), 
    error('Model not on innovation form'), 
end
if nargin<2; kstep=1; end
mpred.A=m.A-m.B*m.C;  
mpred.B=m.B; 
mpred.C=m.C*m.A^kstep;; 
mpred.D=0;

