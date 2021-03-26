% TSRT78 Signal Processing
% Useful code to be used in exam

%% ============================== General code ==========================
%% Read data, plot time and frequency

signal = load('y.mat');
y = signal.y;
N = length(y);
fs = 100; % [Hz]
T = 1/fs; % [s]
timeVec = (0:N-1)*T;
freqHz = (0:N-1)*fs/N;

figure(1);
subplot(2,1,1); plot(timeVec, y); xlabel('Time [s]');
Y = fft(y);
subplot(2,1,2); plot(freqHz,Y); xlabel('Frequency [Hz]');

% Frequency resolution
freqResolution = 1/(N*ts);

%% Energy of signal
% Calculate energy in the time domain
totalEnergy_t = 0;
for i = 1:length(y)
    totalEnergy_t = totalEnergy_t + abs(y(i))^2; 
end

% Calculate energy in the frequency domain
totalEnergy_f = 0;
for i = 1:length(Yf)
    totalEnergy_f = totalEnergy_f + abs(Yf(i)^2)/nSamp;
end

%% ========================= PART 1: Transforms =====================
%% Zero-padding signal for better frequency resolution (from exercise 2.22)
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

%% Decimation (from exercise 2.25)
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

%% Filtering with a BP (from exercise 9.7)

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

%% Filtering with BP (from lab1 ass4)
yFiltered = bandpass(y, [dominantFreq - sth, dominantFreq + sth], fSamp);



%% ======================= PART 2: Model estimation ===================

%% Estimation with AR-model (from lessons)
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
dist = abs(poles);
angleRad = angle(poles);

%% Estimation with AR-model (from lab 1 ass 5)
% AR models for vowel 'a' of order 2
estIdx_x = floor(2*Nx/3);
modelOrder_x = 2;
arMod_x_2 = ar(x(1:estIdx_x), modelOrder_x, 'Ts', Ts);

%% Validation of AR-model (from lab 1 ass 5)
% Validation method 1: compare power spectra
f = 0:0.05:1;
Phi1 = arMod_x_18.NoiseVariance*Ts*abs(freqz(1, arMod_x_18.a, pi*f)).^2;
Phi2 = arMod_x_9.NoiseVariance*Ts*abs(freqz(1, arMod_x_9.a, pi*f)).^2;
Phi3 = arMod_x_2.NoiseVariance*Ts*abs(freqz(1, arMod_x_2.a, pi*f)).^2;
[Phi4, f4] = sig2blackmantukey(x(estIdx_x+1:end), 30, Ts);
figure;
subplot(2,1,1);
semilogy(fs*f/2,Phi1, fs*f/2,Phi2,fs*f/2, Phi3, f4/2, Phi4);
legend('AR(18)', 'AR(9)', 'AR(2)', 'Blackman-Tukey Estimate');
xlabel('Frequency (Hz)'); title('Spectrum Estimates for a');

% Validation method 2: residual whiteness test
eps_x_18 = pe(arMod_x_18, x(estIdx_x+1:end)); 
[Rex_18, k] = sig2crosscovfun(eps_x_18, x(estIdx_x+1:end));
eps_x_9 = pe(arMod_x_9, x(estIdx_x+1:end)); 
Rex_9 = sig2crosscovfun(eps_x_9, x(estIdx_x+1:end));
eps_x_2 = pe(arMod_x_2, x(estIdx_x+1:end)); 
Rex_2 = sig2crosscovfun(eps_x_2, x(estIdx_x+1:end));

eps_y_10 = pe(arMod_y_10, y(estIdx_y+1:end)); 
Rey_10 = sig2crosscovfun(eps_y_10, y(estIdx_y+1:end));
eps_y_5 = pe(arMod_y_5, y(estIdx_y+1:end)); 
Rey_5 = sig2crosscovfun(eps_y_5, y(estIdx_y+1:end));
eps_y_3 = pe(arMod_y_3, y(estIdx_y+1:end)); 
Rey_3 = sig2crosscovfun(eps_y_3, y(estIdx_y+1:end));

figure;
subplot(2,1,1);
plot(k, Rex_18, 'xm', k, Rex_9, '--b', k, Rex_2, '.-r');
legend('Rea(k) for AR(18)', 'Rea(k) for AR(9)', 'Rea(k) for AR(2)');
title('Cross correlation of residuals for a');
xlabel('k'); ylabel('Rea(k)');
subplot(2,1,2);
plot(k, Rey_10, 'xm', k, Rey_5, '--b', k, Rey_3, '.-r');
legend('Reo(k) for AR(10)', 'Reo(k) for AR(5)', 'Reo(k) for AR(3)');
title('Cross correlation of residuals for o');
xlabel('k'); ylabel('Reo(k)');


%% AIC/BIC model order plot
% Select best model order
figure;
maxOrder = 50;
arorder(x, maxOrder);
legend('Minimal loss function',...
       'Akaikes information criterion (AIC)',...
       'Akikes information criterion B (BIC)');
title('Order selection analysis of vowel a');
xlabel('Model order');
