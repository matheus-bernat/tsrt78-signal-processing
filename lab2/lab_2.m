% Lab 2

addpath ../dsp
addpath ../CourseLib

%% Load recorded signals from microphones A and B
chirp_a = load('chirp_a.mat').chirp_a;
chirp_b = load('chirp_b.mat').chirp_b;
sine_a = load('sine_a.mat').sine_a;
sine_b = load('sine_b.mat').sine_b;
multisine_a = load('multisine_a.mat').multisine_a;
multisine_b = load('multisine_b.mat').multisine_b;
white_a = load('white_a.mat').white_a;
white_b = load('white_b.mat').white_b;

% Normalize signals
chirp_a = chirp_a/max(chirp_a);
chirp_b = chirp_b/max(chirp_b);
sine_a = sine_a/max(sine_a);
sine_b = sine_b/max(sine_b);
multisine_a = multisine_a/max(multisine_a);
multisine_b = multisine_b/max(multisine_b);
white_a = white_a/max(white_a);
white_b = white_b/max(white_b);

%% Define sampling frequency and time
fs = 8000; % [Hz]
Ts = 1/fs;

%% ============================== Determine nk ============================
%% With geometry: 
% ANSWER: nk = 8000*0.0012 = 9.6

%% With time-plots: 
% ANSWER: by looking at the multisine plots for mics A and B we
% could tell that the signal comes ~9 samples later to mic A => nk = 9

figure(2);
subplot(2,1,1);
plot(chirp_a); title('Chirp a');
subplot(2,1,2);
plot(chirp_b); title('Chirp b');

figure(3);
subplot(2,1,1);
plot(sine_a); title('Sine a');
subplot(2,1,2);
plot(sine_b); title('Sine b');

figure(4);
subplot(2,1,1);
plot(multisine_a); title('Multisine a');
subplot(2,1,2);
plot(multisine_b); title('Multisine b');

figure(5);
subplot(2,1,1);
plot(white_a); title('White a');
subplot(2,1,2);
plot(white_b); title('White b');

%% With xcorr
% ANSWER: the correlation is high for lag = 9, therefore nk = 9. This can
% be seen for all correlation plots besides for the sine (where the 
% correlation also is a sine). 

maxLag = 70;
[R,lags] = xcorr(chirp_a,chirp_b,maxLag);
figure(6);
subplot(4,1,1);
plot(lags,R); xlabel('Lags'); title('Chirp correlation for different lags'); 

[R,lags] = xcorr(sine_a,sine_b,maxLag);
subplot(4,1,2);
plot(lags,R); xlabel('Lags'); title('Sine correlation for different lags'); 

[R,lags] = xcorr(multisine_a,multisine_b,maxLag);
subplot(4,1,3);
plot(lags,R); xlabel('Lags'); title('Multisine correlation for different lags'); 

[R,lags] = xcorr(white_a,white_b,maxLag);
subplot(4,1,4);
plot(lags,R); xlabel('Lags'); title('White correlation for different lags'); 

%% Estimate model parameters with MyLMS

nk = floor(0.0012*fs);          % time-delay in sample
nb = 40;                        % nb + nk = model order 
mu = 0.015;                     % step-length
lambda = 0.005;                 % leakage factor
% nk = 0; % to estimate time-delay nk below
[th_chirp,sHat_chirp,err_chirp] = MyLMS(chirp_a, chirp_b, nb, nk, mu, lambda);
[th_sine,sHat_sine,err_sine] = MyLMS(sine_a, sine_b, nb, nk, mu, lambda);
[th_multisine,sHat_multisine,err_multisine] = MyLMS(multisine_a, multisine_b, nb, nk, mu, lambda);
[th_white,sHat_white,err_white] = MyLMS(white_a, white_b, nb, nk, mu, lambda);

figure(1);
subplot(4,1,1);
plot(th_chirp); xlabel('Iteration'); title('Chirp model parameters');
subplot(4,1,2);
plot(th_sine); xlabel('Iteration'); title('Sine model parameters');
subplot(4,1,3);
plot(th_multisine); xlabel('Iteration'); title('Multisine model parameters');
subplot(4,1,4);
plot(th_white); xlabel('Iteration'); title('White model parameters');

%% Having a small nk and a large model order:
% ANSWER: In white and multisine one can see that the first nonzero model
% parameter is parameter 9. Therefore nk = 9. In the other plots (sine and 
% chirp) we couldn't see the same effect as all model parameters starting
% from the first one are non-zero.

figure(7);
subplot(4,1,1); stem(1:nb+1,  th_chirp(16000,:)); title('Chirp');
subplot(4,1,2); stem(1:nb+1,  th_sine(16000,:)); title('Sine');
subplot(4,1,3); stem(1:nb+1,  th_multisine(16000,:)); title('Multisine');
subplot(4,1,4); stem(1:nb+1,  th_white(16000,:)); title('White');

%% ========================= Determine best nb ============================

%% ---------------- With AIC/BIC/loss plots
figure(8);
maxOrder = 50;
subplot(4,1,1);
arorder(chirp_a,maxOrder);
legend('Minimal loss function',...
       'Akaikes information criterion (AIC)',...
       'Akikes information criterion B (BIC)');
title('Chirp model order');

subplot(4,1,2);
arorder(sine_a,maxOrder);
legend('Minimal loss function',...
       'Akaikes information criterion (AIC)',...
       'Akikes information criterion B (BIC)');
title('Sine model order');

subplot(4,1,3);
arorder(multisine_a,maxOrder);
legend('Minimal loss function',...
       'Akaikes information criterion (AIC)',...
       'Akikes information criterion B (BIC)');
title('Multisine model order');

subplot(4,1,4);
arorder(white_a,maxOrder);
legend('Minimal loss function',...
       'Akaikes information criterion (AIC)',...
       'Akikes information criterion B (BIC)');
title('White model order');

% ANSWER: Chirp requires a model order of 5, sine requires model order 2, 
% multisine requires model order 40 and white is 10. For higher model 
% orders than those above, no significant decrease in the AIC/BIC is 
% noticeable. 

%% ======================= 3. Determine energy of sHat ====================

chirp_energy = sum(abs(chirp_a).^2);
sine_energy = sum(abs(sine_a).^2);
multisine_energy = sum(abs(multisine_a).^2);
white_energy = sum(abs(white_a).^2);

sHat_chirp_energy = sum(abs(sHat_chirp).^2);
sHat_sine_energy = sum(abs(sHat_sine).^2);
sHat_multisine_energy = sum(abs(sHat_multisine).^2);
sHat_white_energy = sum(abs(sHat_white).^2);

sHat_chirp_energy/chirp_energy             % 0.1621
sHat_sine_energy/sine_energy               % 0.1519
sHat_multisine_energy/multisine_energy     % 0.7462
sHat_white_energy/white_energy             % 0.6984

%% ============== 4. Compare real signal and its estimation ===============

yHat_chirp = chirp_a - sHat_chirp;
yHat_sine = sine_a - sHat_sine;
yHat_multisine = multisine_a - sHat_multisine;
yHat_white = white_a - sHat_white;

clf;
figure(9);
subplot(4,1,1);
plot(chirp_a); hold on; plot(yHat_chirp); title('Chirp: real and estimation');
legend('y','yHat');
subplot(4,1,2);
plot(sine_a); hold on; plot(yHat_sine); title('Sine: real and estimation');
legend('y','yHat');
subplot(4,1,3);
plot(multisine_a); hold on; plot(yHat_multisine); title('Multisine: real and estimation');
legend('y','yHat');
subplot(4,1,4);
plot(white_a); hold on; plot(yHat_white); title('White: real and estimation');
legend('y','yHat');

% ANSWER: For the chirp and sine signals, the estimation is closer to thhe
% real signal. This was expected, since the ratio of energies y-yhat/y was 
% low then, close to 15%. For the multisine and white signals, the
% estimation was bad, which also was expected since the energy ratio
% between the y-yhat/y was over 69%. Therefore the easiest to suppress the
% noise from are sine and chirp, since we can estimate that noise well.


%% ============= PART 2: Adaptive noise cancellation with LMS =============

playsound('multisine')

% 1. How does different step lengths and leakage factors affect the result?
% Can you hear the difference?

% ANSWER: With a step length of 2^-3 (the largest possible) and leakage
% factor on the low end between 0 and 0.03 good noise reducing was
% achieved. A relatively large value of step-length tells how fast the 
% values of theta are adjusted for errors. A low leakage factor tells that
% the earlier estimates of theta are good.

% ANSWER: Assuming the type of noise is unknown: choose nb to the largest
% required by the types of noise above, i.e. nb = 40, and nk = 0.

% ANSWER: The easiest to cancel out was sine by far. The others sucked.



