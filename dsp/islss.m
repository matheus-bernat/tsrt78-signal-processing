function [b,str]=islss(m)
%Check if all matrices are feasible
b=true; str='';
if ~isfield(m,'A') |~isfield(m,'B')  |~isfield(m,'C') 
   b=false, str='A,B,C are mandatory fields'; return
end
if ~isfield(m,'Q') | ~isfield(m,'R')
   b=false, str='Q,R are mandatory fields'; return
end
try  
   m.C*m.A*m.B; m.B*m.Q*m.B'; m.C'*m.R*m.C;
catch 
   b=false, str='A, B, C, Q, R are not compatible'; return
end
if ~all(m.Q==m.Q') | ~all(eig(m.Q)>0)
   b=false; str='Q must be positive definite'; return
end
if ~all(m.Q==m.Q') | ~all(eig(m.Q)>0)
   b=false; str='R must be positive definite'; return
end
if isfield(m,'S')
  try m.B*m.S*m.R; 
     catch b=false; str='S is not compatible'; return
  end
end
