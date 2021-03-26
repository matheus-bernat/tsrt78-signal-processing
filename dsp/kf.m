function [xf,Pf,xp,Pp]=kf(m,y)
%Kalman filter and predictor
N=length(y);
xhat=m.x0;
P=m.P0;
xp(1,:)=xhat; Pp(1,:,:)=P;
for t=1:N
   K=P*m.C'*inv(m.C*P*m.C'+m.R);
   xhat=xhat+K*(y(t)-m.C*xhat);
   P=P-K*m.C*P;
   xf(t,:)=xhat; Pf(t,:,:)=P;
   xhat=m.A*xhat;
   P=m.A*P*m.A'+m.Q;
   xp(t+1,:)=xhat; Pp(t+1,:,:)=P;
end

