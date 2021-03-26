function th=arlms(y,na,mu)
%LMS for adaptive AR parameter estimation
if nargin<3; mu=0.01; end
N=length(y);
th=zeros(na,N);
for t=na+1:N
    phi=-y(t-1:-1:t-na,:); 
    th(:,t)=th(:,t-1)+mu*phi*(y(t)-phi'*th(:,t-1));
end