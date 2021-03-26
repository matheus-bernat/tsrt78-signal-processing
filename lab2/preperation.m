

addpath('../dsp')
addpath('../CourseLib')

%% Preparation 1

N = 1e4;
n = (0:N-1);
h = [0 0 0 0.9 -0.4 0 0 0.2];
lambda = 1;
u = sqrt(lambda)*randn(N, 1); 
y = filter(h, 1, u);
figure;
subplot(2,1,1);
plot(n, u);
subplot(2,1,2);
plot(n, y);

%% Preparation 2 (Exercise 9.6)
% Functions at the bottom of the file.

%% Preparation 3
% Code in MyLMS.m file


%% Preparation 4
[th, s_hat, err] = MyLMS(y, u, 7, 0, 0.015, 0.005);
figure;
plot(th)
xlabel('Iteration');
figure; 
stem(1:8, th(5000, :));
xlabel('parameter');

% Seems to work well

%% Preparation 5
% Done!

%% Questions

% ------------------------------ 1) ------------------------------ %
% It is possible to estimate the time delay for different paths but
% it's basically impossible to find all time delayed signals coming in
% using only geometry.

% It's impossible to find the time delays by comparing the observed signal
% and the real one.

% It is possible to find the time delays by looking at the peaks of the
% cross-correlation function. E.g if Ryu(i) is large then it is likely that
% u(t-i) is a component in the received signal y. One drawback is that you
% can't know what the coefficients in front of the delayed signals are.

% It is possible, look at preparation exercise 4 for example. There we knew
% the model order needed to be 7 but if we were unsure we could have chosen
% a larger value and gotten the same result, with some extra parameter 
% being close to 0 of course. A drawback of this method is that we don't
% know the model order so to be sure we might have to use really large
% model order. 


% ------------------------------ 2) ------------------------------ %
% The step length affects how fast the parameter convergences and the
% variance once it has converged. If it's small the convergence will
% take longer but the variance will be smaller, and vice versa. LMS doesn't
% work for an arbitrary step length since once the equation 
% |eig(I-µ*phi*phi')| >= 1 the system becomes unstable (poles outside unit
% circle).

% ------------------------------ 3) ------------------------------ %
% To evaluate a model one can look at epsilon(k) = y(k)-phi(k)'*theta(k-1).
% This gives an indication of how the real model have changed since the
% last time instance. If it's high the current parameters are no good
% anymore.


% ------------------------------ 4) ------------------------------ %
% Block diagram on paper

%% Functions
% For exercise 9.6 b)
function [theta] = rekid_lms(mu, order, y, u)
    N = length(y);
    theta = zeros(order+1, N+1);
    yhat = zeros(N, 1);
    for t = 1:N
        phi = [u(t:-1:max(t-order,1)); zeros(order+1-t, 1)];
        K = mu*phi;
        theta(:, t+1) = theta(:,t) + K*(y(t) - phi'*theta(:, t)); 
        yhat(t) = phi'*theta(:, t+1);
    end
    theta = theta(:, 2:end)';
end 


% For exercise 9.6 c)
function [theta] = rekid_rls(lambda, order, y, u)
    N = length(y);
    theta = zeros(order+1, N+1);
    P = eye(order+1, order+1);
    yhat = zeros(N, 1);
    for t = 1:N
        phi = [u(t:-1:max(t-order,1)); zeros(order+1-t, 1)];
        P = (P - (P*(phi*phi')*P)/(lambda + phi'*P*phi))/lambda;
        K = P*phi;
        theta(:, t+1) = theta(:,t) + K*(y(t) - phi'*theta(:, t)); 
        yhat(t) = phi'*theta(:, t+1);
    end
    theta = theta(:, 2:end)';
end 


