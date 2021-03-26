function th=sig2armainit(y,na,nc)
%Initialize parameters in an ARMA model
th=sig2ar(y,10*(na+nc));  % Estimate high order AR model
u=randn(1000,1);          % Deterministic input
y=filter(1,[1;th],u);     % Simulation of the AR model
th=sig2arx(y,u,na,nc+1);  % Estimatin of ARX model
th(na+1:end)=th(na+1:end)/th(na+1); % Normalize c0
th(na+1)=[];              % Remove c0
if nargout==0
    disp(['[a,c] = ',num2str([th'])])
end
