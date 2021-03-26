% Lessons in TSRT78 - Signal Processing

% Add paths to files used in lessons
addpath C:/Users/mathe/liu/tsrt78-signal-processing/CourseLib
addpath C:/Users/mathe/liu/tsrt78-signal-processing/dsp
addpath C:/Users/mathe/liu/tsrt78-signal-processing/
addpath C:/Users/mathe/liu/tsrt78-signal-processing/lab2/

%% Old exams

% 2020-08-28
% Question 5
load polydata % x, y
N = length(y);

maxOrder = 10;
phi = ones(N, 1);
AIC = zeros(maxOrder, 1);
for m = 1:maxOrder
    phi = [phi x.^m];
    [thHat,P,lam(m),epsi] = sig2linmod(y,phi);
    AIC(m) = 2*m + N*log(sum(epsi.^2)); % not really like the book's definition in page 256 but it's like this in the answers
end

figure(1);
clf;
plot(AIC); xlabel('Order'); ylabel('AIC value'); title('AIC');
% Model order 3 is good. A "knee" is observed then.

% Get the 4 model parameters for a model of order 3
phi = [ones(N,1) x x.^2 x.^3];
[thHat,P,lam,epsi] = sig2linmod(y,phi);
thHat 

% Visualize the simulated model at values -1...1
xe = (-1:0.01:1)';
N = size(xe,1);
yHat = zeros(N, 1);
for i = 1:N
    yHat(i) = thHat(1) + thHat(2)*xe(i) + thHat(3)*xe(i)^2 + thHat(4)*xe(i)^3;
end
plot(xe, yHat, 'k.'); title('Estimated model');
hold on;
scatter(x, y, 'r');






%% -------------- Lesson 10: Adaptive filters ----------------

%% 9.4
signal = load('sig92.mat');
y = signal.y;
th_real = signal.th;
order = 1;

% LMS - least mean square (pg 336) - minimizes expected value of squared prediction error 
% th(t) = th(t-1) + mu*phi(t) *(y(t)-phi'(t)*th(t-1))

mu = 0.02;  
% 0.01 is a normal start value step-size. The greater the step-size is, 
% the faster the parameters change when the prediction error is larger. Too 
% large step-sizes may give unstable systems.
th_lms = arlms(y, order, mu);

% RLS - recursive least squares (pg 340) - minimizes sum over all squared
% errors, where each error is valued differently depending on how old it is
% through the so called "forgetting factor"
lambda = 0.95; % forgetting factor, between 0 and 1. The smaller it is the less we care about old errors -> the faster the system adapts to new errors
[th_rls,Pt,epsi] = arrls(y, order, lambda);

% KF (pg 345)
Q = 0.001;
[th_kf,P_kf] = arkf(y, order, Q);

figure(1); clf;
plot(th_real, '-', 'LineWidth', 1.5); hold on;
plot(th_lms, '-', 'LineWidth', 1.2); hold on;
plot(th_rls, '-', 'LineWidth', 1.2); hold on;
plot(th_kf, '-', 'LineWidth', 1.2); 
legend('Real', 'LMS', 'RLS', 'KF');
title('Adaptive model parameters estimation');

%% 9.7
% y(t) = s(t) + d(t), where d(t) = b0*u(t) + b1*u(t-T) +...+ bnb*u(t-nb*T)

% ---------- Setup variables in problem
signal = load('ecgdata.mat');
y = signal.ekg;             % disturbed electrocardiogram signal: y = s + d, where d is noise
s = signal.signal;          % true ECG signal
fs = signal.fsamp;          % Sampling frequency in Hz
Ts = 1/fs;
N = length(y);
Y = fft(y); 
freqHz = (0:N-1)*fs/N;

% a) Suppress the noise (that has frequency 50Hz) with a BP-filter

% ----------- Design BP-filter
nyquistFreq = fs/2;
Wp = [45 55]/nyquistFreq;          % passband cutoff frequency, normalized between 0 and 1, where 1 is the Nyquist frequency pi radians/sample
Ws = [40 60]/nyquistFreq;          % stopband cutoff frenquency

passbandRipple = 3;                % [db]
stopbandAttenuation = 40;          % [db]
[minOrder,Wn] = buttord(Wp, Ws, passbandRipple, stopbandAttenuation); % n is the lowest order of a Butter filter that achieves the specifications, Wn are the achieved cutoff frequencies
[b,a] = butter(minOrder, Wn, 'stop'); % returns the coefficients of the BP-filter with the specified order and cutoff frequencies

ekgFiltered = filter(b, a, y);
EKG = fft(ekgFiltered);

% ----------- Plot filtered signal in time and frequency
figure(1);
subplot(2,1,1);
plot((0:N-1)*Ts, y); xlabel('Time [s]'); hold on;
plot((0:N-1)*Ts, ekgFiltered); 
title('EKG in time');
legend('Original', 'Filtered');

subplot(2,1,2);
plot(freqHz(1:(N-1)/2), Y(1:(N-1)/2)); hold on; % only plot up to Nyquist frequency
plot(freqHz(1:(N-1)/2), EKG(1:(N-1)/2));
xlabel('Frequency [Hz]'); 
title('EKG in frequency');
legend('Original', 'Filtered');

% b

% ----------- Generate noise with frequency 50 [Hz] and random phases
f0 = 50;                    % sinusoid frequency, in [Hz]
amp = 10;                   % sinusoid amplitude
noise = (amp*sin(f0*[0:N-1]*Ts + randn))';

% ----------- Adapt a filter with LMS
mu = 0.001;                 % step-size

order = 3;
leakageFactor = 0;
[th,yHat] = MyLMS(y, noise, order, 0, mu, leakageFactor);

% How to know which order is necessary? How to know if the disturbance was
% eliminated?

figure(101);
plot((0:N-1)*Ts, s); hold on;
plot((0:N-1)*Ts, y-yHat, 'r'); % y - yHat = estimate of ekg, since yHat is an estimate of the noise component
xlabel('t [s]');
legend('True signal', 'Estimate of true signal');

%% 9.13
signal = load('volvo.mat');
s = signal.s;   % wheel slip between tire and road
M = signal.M;   % normalized engine torque

z = [s [M ones(313,1)]]; % from 9.13a


nk = 0;
nn = [na, nb, nk];

[th, yhat] = rarx(z, nn, adm, adg); 



%% ------------------ Lesson 9: Kalman filter -----------------------
%% 8.3

% b
A = -0.5; 
B = 1; 
C = 1; 
D = 1; 
Q = 0.25; 
R = 0.25;
[P,K] = stat_kalman(A,B,C,D,Q,R);

%% 8.17

w = 0.1*pi;
A = [1 0; 0 1]; 
B = 0;
C = [1 0];

% Design parameters Q and R:
R = eye(2);
P0 = eye(2);

% See function kf_1_step_predictor in the bottom of this file
ym = load('sig80').ym;

[xHat, phiHat] = phaseEstimator(ym, A, B, C, R, P0);



%% -------------- Lesson 8: Wiener and Kalman filter ----------------
%% 8.15
clear;
signal = load("lunarmodule.mat");
trueX = signal.x;                       % True state vector
trueY = signal.y;                       % True position (= x(:,1))
ylm = signal.ylm;                       % Position observations
yhm = signal.yhm;                       % Velocity observations

% Model parameters
A = [1 1; 0 1];                
% System matrix
B = [0.5; 1];                           % Noise gain matrix
C = [1 0];                              % Observation matrix

Qtilde = 1;                             % Normalized variance of process noise u(t). 
Q = B * Qtilde * B';                        % Cov matrix of process noise Bu(t)
% Note: only the ratio between Q and R matter in the Kalman filter. That's
% why we set Qtilde to 1 and so only R needs to be tuned.
R = 500; % what to set R to ???????????????????????????????????????????????????????????????

% Kalman filter

[y_hat,x_hat,P] = kalmanfilt(ylm, A, B, C, Qtilde, R);
% Now that P0 and x0 are omitted, the stationary KF is computed!


clf;figure(1);
subplot(2,1,1);
plot(trueX(:,1));xlabel('time');ylabel('Position');
hold on;
plot(x_hat(:,1));
legend('True position','Kalman estimate');

subplot(2,1,2);
plot(trueX(:,2));xlabel('time');ylabel('Velocity');
hold on;
plot(x_hat(:,2));
legend('True velocity','Kalman estimate');

% ------------- b

P0 = diag([1e5,5]);         % Since x0 is just a guess,  let P0 be large
x0 = zeros(2,1);            % Initial guess for x0
[y_hat_b,x_hat_b,P_b] = kalmanfilt(ylm, A, B, C, Qtilde, R, P0, x0);

figure(2);
subplot(2,1,1);
plot(trueX(:,1));xlabel('time');ylabel('Position');
hold on;
plot(x_hat_b(:,1));
legend('True position','Kalman estimate');

subplot(2,1,2);
plot(trueX(:,2));xlabel('time');ylabel('Velocity');
hold on;
plot(x_hat_b(:,2));
legend('True velocity','Kalman estimate');

% ------------- c: Start state is known
P0 = diag([5,5]);           % Since x0 is known, let P0 (prediction error variance) be zero
x0 = trueX(1,:)';            
[y_hat_c,x_hat_c,P_c] = kalmanfilt(ylm, A, B, C, Qtilde, R, P0, x0);

figure(3);
subplot(2,1,1);
plot(trueX(:,1));xlabel('time');ylabel('Position');
hold on;
plot(x_hat_c(:,1));
legend('True position','Kalman estimate');

subplot(2,1,2);
plot(trueX(:,2));xlabel('time');ylabel('Velocity');
hold on;
plot(x_hat_c(:,2));
legend('True velocity','Kalman estimate');

% ANSWER: Known initial states (x0) enables us to use relatively low
% initial prediction error covariance P0, for example P = [0 0; 0 0].


% -------------- d: Velocity is also observed
C_d = [1 0; 0 1];               % Both position and velocity is also observed. Note that C = [1 1] gives the wrong result (y = position + velocity instead of y = [position; velocity])
P0 = diag([0.1,0.1]);           % Since x0 is known, let P0 (prediction error variance) be zero
x0 = trueX(1,:)'; 
R = diag([500,2]);              % FORGOT that R should also be changed if now y is a 2x1 vector (containing position AND velocity!)
[y_hat_d,x_hat_d,P_d] = kalmanfilt([ylm yhm], A, B, C_d, Qtilde, R, P0, x0);

figure(4);
subplot(2,1,1);
plot(trueX(:,1));xlabel('time');ylabel('Position');
hold on;
plot(x_hat_d(:,1));
legend('True position','Kalman estimate');

subplot(2,1,2);
plot(trueX(:,2));xlabel('time');ylabel('Velocity');
hold on;
plot(x_hat_d(:,2));
legend('True velocity','Kalman estimate');

% COMMENT: the velocity estimate is now a lot better, which is natural
% given that we can measure the velocity at each time step.


% Kalman filter
%N = size(trueY, 1);                     
%x_hat = zeros(2, N);                    % Allocate memory for states (1st row: positions in each time, 2nd row: velocities)
%x_hat(:, 1) = zeros(2,1);               % Initial state estimate

%P = zeros(2,2,N);                       % Allocate memory for state error covariance
%P(:,:,1) = diag([1, 0.1]);              % Initial state covariance P0
% Not necessary, s stationary KF has only one value of state error
% covariance that is P_bar.

% Stationary prediction error variance P_bar
%syms P_bar;
%equation = P_bar == (A*P_bar*A') - (A*P_bar*C'*C*P_bar*A')/(C*P_bar*C' + R) + Q;
%P_bar = solve(equation, P_bar);



% QUESTIONS:
% - how to choose P0 in the case of stationary KF? Calculate myself with
% formula? 
% ==========> ANSWER: guess a value according to how good you
% think that the initial estimate for the initial state x0 is.
% - why to start the state as the estimate with Cholesky thing?
% - what is "y_hat" that is returned by "kalmanfilt"?

%% -------------- Lesson 7: Wiener and Kalman filter ----------------
%% 7.8

% QUESTIONS:
% - when to use filter and filtfilt? 
%      (It seems like filtfilt is used for noncausal and filter for causal
%      filters.)
% - what is the theoretical difference between filtfilt and filter
%   functions in MATLAB?
% - when taking variance of 1-step smoothing, why are boundaries 1:end-1 
%   and 2:end?

% OBSERVATIONS:
% - be careful not to confuse numerator with denominator in filter and
%   filtfilt
% - when using filtfilt, first factorize H(z) = H(z)*H(z^-1) and send only
%   the causal part to filtfilt. See last slide in lecture 7 and page 162. 

%%
clear;
signal = load("sig70.mat");
y = signal.y;
s = signal.s;
%% shat = y
variance = var(s - y); % 1.0601

%% Non-causal WF
num = sqrt(0.45/2);
den = [1 -0.5];
shat = filtfilt(num, den, y);
variance = var(s - shat); % 0.2663

%% Causal WF
% Here we use filter instead of filtfilt as it is a causal filter (???)
num = 0.375;
den = [1 -0.5];
shat = filter(num, den, y);
variance = var(s - shat); % 0.3542

%% Causal WF 1-step ahead predictor
num = [0 0.3];
den = [1 -0.5];
shat = filter(num, den, y);
variance = var(s - shat); % 0.5567

%% Causal WF 1-step smoothing
num = [0.188 0.225];
den = [1 -0.5]; % my coefficients are actually equal to [0 1 -0.5]
shat = filter(num, den, y);
variance = var(s(1:end-1) - shat(2:end)); % 0.2879 =========> WHY 1:end-1 and 2:end ?

%% Causal FIR WF of 1st order
num = [0.4048 0.2381];
den = [1 0];
shat = filter(num, den, y);
variance = var(s - shat); % 0.385
% OBSERVATION: careful not to exchange places between NUM and DEN!!!!

%% Causal FIR WF of 0th order
num = [0.5];
den = 1;
shat = filter(num, den, y);
variance = var(s - shat); % 0.4696

%% -------------- Lesson 6: Wiener filters ----------------
clear;
% Good MATLAB functions: 
% residue() - finds partial-fraction expansion
% roots() - finds roots of equation
% toeplitz() - constructs matrix to use in linear equation system

% Partial-fraction expansion in exercises 7.2 nd 7.7
%[r,p,k] = residue([0 1],[1 3 1]);
%[r,p,k] = residue([0 1],[1 0.5+2.618 0.5*2.618]);
%[r,p,k] = residue([0 1],[1 -2.5 1]);
%[r,p,k] = residue([0 1],[1 -2.8 0.8*2]);
[r,p,k] = residue([1 0],[1 -2.8 .8*2]);
A = [2 .8 ; .8 2];
b = [1; .8];
filterCoeffs = A\b;


%% -------------- Lesson 5: Signal models & estimation ----------------
%% 6.5 - AR-models in noisy signals
clear;clf;
% ------ a 
% Generate 300 samples of sinusoid with angular frequency 0.1 [rad/s]
w0 = 0.1; % rad/s
N = 300;
n = 0:N-1;
x = sin(w0*n);
figure(1);subplot(1,2,1);
plot(n, x);xlabel("n");ylabel("x[n]");

% Estimate second model order AR-model for sinusoid
modelOrder = 2;
arMod = ar(x, modelOrder);
subplot(1,2,2);
zplane(1,arMod.A); 
% OBS: arMod.A contains the polynomial in the denominator of the transfer 
% function. Therefore, its roots are the poles of the system
poles = roots(arMod.A);
dist = abs(poles)
angleRad = angle(poles)

% ANSWER: Relation between angle and angular frequency of sinusoid is equal
% to the angular frequency of the signal! ?

% ------ b
% Add white noise to signal with variance 0.01
% Generate 300 samples of sinusoid with angular frequency 0.1 [rad/s]
varianceNoise = 0.01;
xNoise = sin(w0*n) + varianceNoise*randn(1,N);
figure(2);
subplot(1,2,1); plot(n, xNoise);xlabel("n");ylabel("xNoise[n]");

% Estimate second model order AR-model for sinusoid
modelOrder = 2;
arMod = ar(xNoise, modelOrder);
subplot(1,2,2);
zplane(1,arMod.A); 
% OBS: arMod.A contains the polynomial in the denominator of the transfer 
% function. Therefore, its roots are the poles of the system
poles = roots(arMod.A);
dist = abs(poles);
angleRad = angle(poles); % angles at which the poles are


% ANSWER: The angle of the poles estimated by the AR-model for the noisy
% signal is 0.0971, which differs a little from the angle of the poles
% estimated by the AR-model by the non-noisy signal where the angle is 0.1.

%% 6.18 - validation of model order for AR-models
clear;clf;
% Load signal
signal = load("sig60.mat");
y = signal.y;
figure(1);plot(y);
N = size(y,1);
% Divide data in training and validation data
trainData = y(1:round(2*N/3));
validData = y((round(2*N/3)+1):end);

% -------- a. Get appropriate model order
loss = [];
stds = [];
for order = 1:50
    arMod = ar(trainData,order);
    predErr = pe(arMod, validData); % predErr = y_i - ^y_i for all i's
    totalErr = sum(predErr.^2);
    stds = [stds sum(std(arMod.A))];
    loss = [loss totalErr];
end
figure(2);plot(loss);xlabel("Model order");ylabel("Squared total loss");
figure(3);plot(stds);xlabel("Model order");ylabel("Sum of standard devs");

% As we see in the graph, a model order of 5 is good (since it gives less
% loss than the model orders 6 to 12. If we want an even lower loss, we can
% choose 19, which is less than 20 to 43.
chosenOrder = 2;
arMod2 = ar(trainData,chosenOrder);

% -------- b. Study standard deviations of parameter estimates
% ANSWER: it seems as: the biggest the model order => the smallest the
% standard deviation of the model parameters are.
% OBS: the estimated model parameters (with index greater than 2) seem to 
% have standard deviations that are of the same order of the model 
% parameter itself, which is weird. No 
present(arMod);

% -------- c. Plot models in freq domain together with estimate
clf; figure(4);
bode(etfe(iddata(validData)));
for i = 1:3
    hold on; bode(ar(trainData,i));
end
% ANSWER: model order 1 is bad in comparison to the rest. Choose between
% model order 2 and 3.

% -------- d. Check whether prediction errors are white with cov function

% Idea: if  residuals are white, the covariance function are 0 for time
% indices > 0 (see page 92), since white noise is independent (and has
% therefore covariance 0) of all other signals in different time instants.
arMod1 = ar(trainData,1);
arMod3 = ar(trainData,3);

predErr1 = pe(arMod1,validData);
predErr2 = pe(arMod2,validData);
predErr3 = pe(arMod3,validData);

figure(5);plot(1:size(predErr1),predErr1/predErr1(1));
hold on; plot(1:size(predErr2),predErr2/predErr2(1));
hold on; plot(1:size(predErr3),predErr3/predErr3(1));
title("Covariance function of residuals");
legend("model order 1","model Order 2","model Order 3");

% ANSWER: Model 3 gives white residuals, since the covariance of the
% residual with other time instant residual is equal to zero (definition of
% whiteness). Model order 1 is less white in comparison to 2 (that varies 
% closer to 0).

%% 6.7 - Determine frequency components using parametric method
clear;
signal = load("sig30.mat");y = signal.y;

% Parametric method
modelOrder = 15;
arMod = ar(y, modelOrder);
figure;bode(arMod);
zplane(1, arMod.A);
phi = angle(roots(arMod.A));

% QUESTION: how to know what the main frequencies are? Is it those with the
% greatest distance from the center of the circle? Or just the ones with
% smallest angle?

% Double check with non-parametric method
N = size(y); n = 0:N-1;
y1024 = [y' zeros(1,1024-size(y,1))]; Y1024=fft(y1024);
%figure;
%subplot(2,1,1);plot(1:1024,y1024);
%subplot(2,1,2);plot(1:1024,abs(Y1024));

% ANSWER: from the Bode plot given by the AR(2)-model, it seems as the two
% frequencies that compose the signal y are 0.201 [rad/s] and 0.251 [rad/2]

%% Testing how a prediction of a sinus looks in the zplane
Ts = 1e-3; fs = 1/Ts;
N = 1000;
n = 0:Ts:N-1;
f0 = 0.5; % [Hz]
amp = 0.5;
x = amp*sin(f0*n);
figure;subplot(2,1,1);
plot(n*Ts, x);xlabel("Time (s)");ylabel("x[t]");
fHz = n*fs/N;
subplot(2,1,2);
X = fft(x);
plot(fHz,abs(X));xlabel("Frequency (Hz)");ylabel("X[f]");

modelOrder = 20;
arMod = ar(x,modelOrder);
zplane(1,arMod.A);
abs(roots(arMod.A)) % 1, doesn't matter what the amplitude of sinus is

% CONCLUSION: it looks like a sinus curve, when modeled by a AR(i)-model
% has 2 of its poles ON the unit circle, while the rest is within the unit
% circle. The two poles ON the unit circle form the conjugate pair with the
% smallest angle.

% QUESTION: Why was "the closest the poles are to the unit circle" a purity 
% measure used in lab 1?

%% 6.4 - Construct own function that estimates an AR-model

%  
%function [theta,varMatrix,costFunctionVal] = arestimate(signal, modelOrder)
%    % According to the calculations made in the book page 220-221, the
%    % estimated parameters can be calculated as:
%    N = max(size(signal));
%    
%end

%% -------------- Lesson 3: Spectrum, filter ----------------
%% 2.25: decimation
clf;clear;

% Signal with 60 samples
N = 60;
ts = 1; % This is an assumption, as it's not written in the exercise
n60 = 0:N-1;
w0 = 2*pi/5;
x60 = sin(w0*n60);
X60 = fft(x60);
nr_plots = 4;

fs = 1/ts; % sampling rate
freqHz = n60*fs/N;

subplot(nr_plots,1,1);stem(n60*ts,x60);xlabel("time [s]");ylabel("x60[t]");
subplot(nr_plots,1,2);stem(freqHz,abs(X60)); xlabel("frequency [Hz]");

% Signal with 60/2 = 30 samples
factor = 2;
N30 = N/factor; n30 = 0:N30-1;
x30 = x60(1:factor:end);
X30 = fft(x30);

freqHz = n30*fs/N30;

subplot(nr_plots,1,3);stem(n30*ts,x30);xlabel("time [s]");ylabel("x30[t]");
subplot(nr_plots,1,4);stem(freqHz,abs(X30)); xlabel("frequency [Hz]");

% From the formula of the DFT: F{x[n]} = 1/T * sum .... Thefore, when the
% sampling time T is increased to 2T, the amplitude of the transform
% decreases with a factor 2.

% Let the original sampling frequency of the time signal be equal to f
% before decimation. After decimation (when every other sample is taken),
% the sampling frequency becomes then f/2. Notice that, in order to not
% have aliasing, the sampling frequecy must be equal or greater than 2
% times the greatest frequency component of the time signal: which in our
% case is 2*pi/5. Therefore sampling frequency must be at least 2 * 2pi/5 =
% 4pi/5.  

% QUESTIONS: how can one know the exact sampling frequency? How is the
% amplitude affected by the sampling frequency? How to see if there's
% aliasing or not just by looking at the spectrums?

%% 3.2 Get composing frequencies by decomposing (and zero-padding for better resolution)
clf; figure(1);
signal = load("sig30.mat");
y = signal.y;
N = size(y,1);
ts = 1; fs = 1/ts;
n = 1:N;
Y = fft(y);
freqHz = n*fs/N;
subplot(3,1,1);plot(n*ts,y);xlabel("time [s]");
subplot(3,1,2);stem(freqHz,abs(Y));xlabel("frequency [Hz]");

% Zero-pad:
Y256 = fft([y', zeros(1,256-length(y))]);
subplot(3,1,3);stem((0:255)*fs/256, abs(Y256));xlabel("frequency [Hz]");
title("Zero padded transform");

% ANSWER: Through zero-padding we can easier see that the searched
% frequencies that together compose the signal y are: 0.0297 Hz and 
% 0.04292 Hz.

%% 4.16 Decompose signals by constructing filters with "filter" & "filtfilt"
clf; 
signal = load("sig40.mat");

s1 = signal.s1; % sinus with freq 0.1 rad/s
s2 = signal.s2; % sinus with freq 0.5 rad/s
s = signal.s;   % sum of s1 and s2
N = size(s,2); ts = 1; fs = 1/ts;
n = 0:N-1;
nr_plots = 5;
figure(1);
subplot(nr_plots,1,1);plot(n*ts,s1);xlabel("time [s]");title("s1");
subplot(nr_plots,1,2);plot(n*ts,s2);xlabel("time [s]");title("s2");
subplot(nr_plots,1,3);plot(n*ts,s);xlabel("time [s]");title("s");

S = fft(s);
freqHz = n*fs/N;
subplot(nr_plots,1,4);stem(2*pi*freqHz,abs(S));xlabel("frequency [rad/s]");title("S");

% Separate the sinuses composing s with filters:
filter_order = 3;
[lowB,lowA] = butter(filter_order, 0.1,"low");

filtered_S = filter(lowB,lowA,S);
subplot(nr_plots,1,5);stem(2*pi*freqHz,abs(filtered_S));xlabel("frequency [rad/s]");title("S");

% ANSWER: notice that the higher the order of the filter used in "filter",
% the bigger is the time delay.

% OBS: didn't design filters with filtfilt

%% 4.18
clf; clear;
signal = load("soderasen.mat");
i4r = signal.i4r;
fs = 1000; % sample frequency [Hz]
ts = 1/fs;
N = size(i4r,1); 
n = 0:N-1;
freqHz = fs*n/N;
transf = fft(i4r);

figure(1); nr_plots = 5;
subplot(nr_plots,1,1);plot(n*ts,i4r);xlabel("time [s]");title("i4r");
subplot(nr_plots,1,2);plot(freqHz,abs(transf));xlabel("frequency [Hz]");title("transform");

% Extract fundamental frequency at 50 Hz with a BP Butter filter
filter_order = 3;
[b,a] = butter(filter_order, [40/(0.5*fs) 60/(0.5*fs)]); % cut-off frequencies must be between 0 and 1 (i.e. normalized with frequency fs)
%h = fvtool(b,a);
filtered_i4r = filtfilt(b,a,i4r);
filtered_transf = fft(filtered_i4r);
subplot(nr_plots,1,3);
plot(freqHz, abs(transf), freqHz, abs(filtered_transf));
xlabel("frequency [Hz]");title("original transf and filtered fundamental freq");
legend("original", "filtered");


% First harmonics at 100 Hz
filter_order = 3;
[b,a] = butter(filter_order, [90/(0.5*fs),110/(0.5*fs)]);
first_har = filtfilt(b,a,i4r);
first_har_transf = fft(first_har);
subplot(nr_plots,1,4);
plot(freqHz, abs(transf), freqHz, abs(first_har_transf));
xlabel("frequency [Hz]");title("original transf and 1st harmonic");
legend("original", "filtered");


% Second harmonics at 200 Hz
filter_order = 3;
[b,a] = butter(filter_order, [140/(0.5*fs),160/(0.5*fs)]);
second_har = filtfilt(b,a,i4r);
second_har_transf = fft(second_har);
subplot(nr_plots,1,5);
plot(freqHz, abs(transf), freqHz, abs(second_har_transf));
xlabel("frequency [Hz]");title("original transf and 2nd harmonic");
legend("original", "filtered");

% Power of the filtered signals and total power
power_fundamental_freq = bandpower(filtered_i4r);
power_first_har = bandpower(first_har);
power_second_har = bandpower(second_har);
sum = power_fundamental_freq + power_first_har + power_second_har;
power_total = bandpower(i4r);
diff = power_total - sum;


%% -------------- Lesson 2: DFT/DTFT Spectrum ----------------
%% 2.3
% read:http://paulbourke.net/miscellaneous/dft/#:~:text=More%20precisely%2C%20a%20continuous%20function,when%20digitising%20a%20continuous%20signal.
fs = 4000; % sample frequency [Hz]
ts = 1/fs;
fc = 10; % carrier frequency [Hz]
%t = (0:ts:stop_time); % time vector from 0 to 1 s with ts interval between values. In [s]
x_t = cos(2*pi*fc*t);
N = length(x_t);
time_vector = 1:N;
figure(1);
plot(time_vector * ts, x_t); % plot signal in time-domain
xlabel("Time [s]"); ylabel("x(t)"); title("Signal in time domain, x(t)");

Y_f = fft(x_t);

frequency_vector = ((0:N-1)/N)*fs;

figure(2);
plot(frequency_vector, abs(Y_f));
xlabel("Frequency [Hz]"); ylabel("X(f)"); title("Spectrum of x(t)");

% 3c) Reason for leakage: multiplication in the time domain between an 
% infinite cosine signal with a rectangular window corresponds to a convo-
% lution between the Diracs (trasform of cos) and sincs (transform of 
% rectangular function. The sincs accumulate and therefore we have leakage.
%% 2.5
% a -----------------------
% Take signal defined above as "u2".
N = length(u2);
ts = 3.9e-4;

% Transform u2 and plot it
subplot(2,1,1);
[U2,w] = dtft(u2,ts);
plot(w,abs(U2));

% Take every third sample of u2, transform and plot its transform.
interval = 2;
y = u2(1:interval:end); 
subplot(2,1,2);
[Y,w] = dtft(y,interval*ts);
plot(w,abs(Y));
% b --------------------------
% Transform u2 and plot it
subplot(2,1,1);
[U2,w] = dtft(u2,ts);
plot(w,abs(U2));

% Lowpass filter signal before transforming it
y2 = decimate(u2,3); % reduces sample rate of u2 by factor 3.
subplot(2,1,2);
[Y2,w] = dtft(y2,3*ts);
plot(w, abs(Y2));
% Notice that there is less aliasing when the signal is lowpass filtered
% before sampling.
%% 2.10
n = 0:15;
x0_n = cos((2*pi/8)*n);
x1_n = cos((2*pi/7)*n);
subplot(4,1,1); plot(n, x0_n); title("x0[n]");
X0 = fft(x0_n); % this is the dft of the signal x0_n
subplot(4,1,2); stem(2*pi*n/16, abs(X0)); title("DFT x0"); xlabel("k");
subplot(4,1,3); plot(n, x1_n); title("x1[n]");
X1 = fft(x1_n);
subplot(4,1,4); stem(2*pi*n/16, abs(X1)); title("DFT x1"); xlabel("k");

% In light of figure 2.1 and seeing the plotted dft's, we see that the
% dft of x0_n is zero at other values then at the top since the dtft is 
% sampled at times when it is equal to zero. FACIT: this occurs when
% the DFT frequency axis gradation is an integer factor of the signal an-
% gular frequency. 
% Notice from figure 2.1 that there's leakage in both dfts.

%% 2.22 Zero padding
% Calculate 16 point DFT of x_n:
clf; % Clear old graphs.
N = 16;
n = 0:N-1;
w0 = 2*pi/sqrt(17);
x_n = cos(w0*n);
X_k = fft(x_n);
subplot(5,1,1); stem(n,x_n); title("x[n]");
subplot(5,1,2); plot(2*pi*n/N, abs(X_k)); title("Not zero-padded X[k]");

% Zero-pad x[n] to lengths 32
X32 = fft([x_n, zeros(1,32-length(x_n))]);
subplot(5,1,3); plot(2*pi*(0:31)/32, abs(X32)); hold on; title("Zero-padded X[k] N=32");
plot(2*pi*n/N, abs(X_k), 'r*');

% Zero-pad x[n] to lengths 64
X64 = fft([x_n, zeros(1,64-length(x_n))]);
subplot(5,1,4); plot(2*pi*(0:63)/64, abs(X64)); hold on; title("Zero-padded X[k] N=64");
plot(2*pi*n/N, abs(X_k), 'r*');

% Zero-pad x[n] to lengths 256
X256 = fft([x_n, zeros(1,256-length(x_n))]);
subplot(5,1,5); plot(2*pi*(0:255)/256, abs(X256)); hold on; title("Zero-padded X[k] N=256");
plot(2*pi*n/N, abs(X_k), 'r*');


% Conclusion: by comparing the spectrums achieved before and after
% zero-padding, we notice that the one achieved after zero-padding had more
% points in X[k] and therefore needed LESS INTERPOLATION (i.e. guessing of
% how the signal looks like). More zero-padding -> more continuous
% spectrum.

% QUESTION: why does it get this way just by adding some zeros in the end
% of x[n]????
% ANSWER: zero-padding just increases the sampling interval of DTFT
% (getting then DFT) in the frequency domain. Increasing length N without
% adding more information just results in a denser sampling of the
% underlying DTFT of the signal. See link: https://dsp.stackexchange.com/questions/34919/how-does-zero-padding-affect-the-magnitude-of-the-dft
% and https://dsp.stackexchange.com/questions/34211/is-spectral-leakage-due-to-windowing-different-for-the-dtft-and-dft/34216#34216



%% ================== Function definitions from all lessons ===============

%% Lesson 10
% Exercise 9.6


%% Lesson 9:
% Exercise 8.3
function [Pp_bar,K_bar,x_hat,y_hat] = stat_kalman(A,B,C,D,Q,R,y)    
    % From page 327
    [K_bar,Pp_bar,Pf_bar] = dlqe(A,B,C,Q,R);  % Difference between Pp_bar and Pf_bar ?
    K_bar = real(K_bar);
    Pp_bar = real(Pp_bar);
    % Obs: numerical innacuracy can lead to imaginary gain)
    if nargin > 6
        [y_hat,x_hat] = dlsim(A - A*K_bar*C, C, D, y);
    end
end

% Exercise 8.17

% Kalman filter for one step ahead prediction

function [xHat, phiHat] = phaseEstimator(y,A,B,C,R,P0,w)
% [th,s_hat]=MyLMS(y,u,nb,nk,mu,lambda)
%
% 1 step ahead KF predictor. Two states.
%
% Inputs:       y       N x 1 vector with signal measurements 
%               A       2 x 2 x N state update matrix
%               B       2 x 2 x N noise gain matrix
%               C       1 x 2 x N measurement matrix
%               Q       2 x 2 x N process noise covariance matrix
%
% Outputs       x_hat   2 x N 
%               P       N x 1 state covariance 

    N = length(y);
    xHat = zeros(2, N);         % state
    yHat = zeros(2,N);
    P = P0;                     % initial state covariance
    phaseHat = zeros(1,N);      % predicted phase phi, given x1(t) & x2(t)
    
    for t = 1:N
        % ----- Time update
        xHat(:,t+1) = A*xHat(:,t);
        P = A*P*A';                 % Q is given as 0
        
        % ----- Measurement update
        C = [cos(w*t) sin(w*t)];
        % Define S
        S = inv(C*P*C' + R);
        % First get the gain K
        K = P*C'*S;
        xHat(:,t+1) = xHat(:,t+1) + K*(y(t+1) - C*xHat(:,t+1));
        P =  P - P*C*S*C'*P;
    end
end




