% Homework 3

clear;

%% Causal Wiener filter for 1-step ahead prediction

% Plot frequency and phase spectrum
fs = 5; % [Hz]
nrPoints = 256; % standard
b = 0.3;
a = [-0.5 1];
freqz(b,a,nrPoints,fs);
title('Frequency and phase response of 1-step ahead causal WF predictor');

%% Stationary 1-step ahead Kalman fitler predictor

% Plot frequency and phase spectrum
fs = 5; % [Hz]
nrPoints = 256; % standard
b = 0.3;
a = [-0.5 1];
freqz(b,a,nrPoints,fs);
title('Frequency and phase response of stationary 1-step ahead KF predictor');

% Reference:
% https://se.mathworks.com/help/signal/ref/freqz.html#bt8l9fo-1