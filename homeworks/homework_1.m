% Homework 1 in TSRT78 - Signal Processing

addpath C:/Users/mathe/liu/tsrt78-signal-processing/homeworks

signal = load('signal.mat');
z = signal.z;
Ts = 0.1; fs = 1/Ts;
N = length(z);
n = (0:N-1);
t = Ts*n;
f = fs*n/N; 

figure(1);
subplot(3,1,1);stem(t, z); xlabel('time [s]'); ylabel('z[t]');
Z = fft(z);
subplot(3,1,2); stem(f,abs(Z));xlabel('frequency [Hz]'); ylabel('Z[f]');

% Zero pad the signal to better discern what the two main frequencies are:
zeroPads = 1024;
Zp = fft([z, zeros(1,zeroPads-length(z))]);
subplot(3,1,3);
plot(fs*(0:zeroPads-1)/zeroPads, abs(Zp));
xlabel('frequency [Hz]'); ylabel('Zero padded X[f]');

frequencyResolution = 1/(N*Ts); % 0.1639

% ANSWERS:
% The spectrum of x(t) is given in figure 1. The frequencies f0 and f21 are
% respectively 1.25 [Hz] and 1.65 [Hz]. it is clear that the first
% frequency component corresponds to the highest lobe. The second
% frequency, in turn, must be the second side lobe to the right of the mian
% lobe. This can be derived by looking at the graph and noticing that the
% absolute value of X does not go to zero at the frequency 1.53 [Hz] as
% much as it goes to zero at frequency 1.084 [Hz], which is explained by
% the fact that there is another frequency component at 1.65 [Hz] which
% brings up the value of the total spectrum.

% The frequency resolution is 0.1639.








