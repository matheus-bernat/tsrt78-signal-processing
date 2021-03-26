function [th,s_hat,err]=MyLMS(y,u,nb,nk,mu,lambda)
% [th,s_hat]=MyLMS(y,u,nb,nk,mu,lambda)
%
% A leaky LMS algorithm
%
% Inputs:       y       N x 1 vector with signal measurements (y=s+n)
%               u       N x 1 vector with input values
%               nb      model order
%               nk      number of samples the input should be delayed
%               mu      step length
%               lambda  leakage factor (actually 'gamma' in the book)
%
% Outputs       th      N x (nb+1) matrix with the estimated filter parameters
%               s_hat   N x 1 vector with the filtered signal
%               err     N x 1 vector with the innovation, i.e., err(k)=y(k)-phi(k)'*theta(k-1)
%
% Authors: Caspian Süsskind & Matheus Bernat
% Date: 2020-12-9
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Place your code here
    N = length(y);
    th = zeros(nb+1, N+1);
    s_hat = zeros(N, 1);
    err = zeros(N, 1);
    
    for t = 1:N
        if(t-nk > 0)
            phi = [u(t-nk:-1:max(t-nk-nb, 1)); zeros(nk+nb+1-t, 1)];
        else
            phi = zeros(nb+1, 1);
        end
        th(:, t+1) = (1-lambda)*th(:, t) + mu*phi*(y(t) - phi'*th(:, t));
        err(t) = y(t) - phi'*th(:, t);
        s_hat(t) = y(t) - phi'*th(:, t+1); % shat = y - y_hat
    end
    th = th(:, 2:end)';









