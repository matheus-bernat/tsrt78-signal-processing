function [yhat,xhat,Pp]=kalmanfilt(y,A,B,C,Qtil,R,P0,x0)

%SYNTAX:  [yhat,xhat,P]=kalmanfilt(y,A,B,C,Qtil,R,P0,x0)
%	y=measurements
%	A,B,C=System matrices
%	Qtil=Process noise covariance
%	R=Measurement noise covariance
%	P0=Initial value for the error covariance
%	x0=Initial value for the state estimate
%
%       If P0 and x0 are excluded, the stationary Kalman filter 
%       is computed, otherwise, the time variable Kalman filter
%       is employed. Note that for the time variable Kalman filter,
%       the sequence of covariance matrices are stored after each
%       other in P, hence P(k+1:k+n,1:n)= the covariance matrix 
%       after k-1 time steps (dim A = nxn)
  
if nargin< 7, 
	[L,P,Z]=dlqe(A,B,C,Qtil,R);
	K=real(A*L);
	Pp=real(P);
	D=0*C*C';
	[yhat,xhat]=dlsim(A-K*C,K,C,D,y);
end

if nargin>6,
	n=length(A);
	Ndata=max(size(y));
	xhat=[];
	yhat=[];
	xt=zeros(size(A,1),1);
	Q=B*Qtil*B';
	if nargin > 7, xt=x0;end;
	P=P0;
	for t=1:Ndata,
		xhat=[xhat;xt'];
        yhat=[yhat;(C*xt)'];
        S=inv(C*P*C'+R);
		K=A*P*C'*S;
		P=A*P*A'+Q-A*P*C'*S*C*P*A';
		xt=(A-K*C)*xt+K*y(t,:)';
		Pp(t*n-n+1:t*n,:)=P;
	end
end		


