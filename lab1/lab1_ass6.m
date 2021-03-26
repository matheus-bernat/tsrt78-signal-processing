% TSRT78, Lab 1: Fundamental Signal Processing
% Matheus Bernat (matvi959) & Caspian Süsskind (cassu286)

% ================= 6 Assignment: Speech encoding as in GSM ===============


% ------------------ INTRO ------------------
clear;
% --------------- Read wav file: extract data and sampling frequency
[x,fs] = audioread('cleverFox.wav'); 
N = size(x,1);
segmentLength = 160;
% --------------- Plot signal in time and frequency axes
Ts = 1/fs;
t = Ts*(0:N-1); % time vector in seconds
numSegments = floor(N/segmentLength);
figure(1);clf();
subplot(2,1,1); plot(t, x), xlabel('time [s]'); ylabel('x[t]');
X = fft(x);
f = (0:N-1)*fs/N;

subplot(2,1,2); plot(f, abs(X)) 
xlabel('frequency [Hz]'); ylabel('X[f]');



% --------------- Create 115 AR models, one for each 160 point segment
modelOrder = 16;
sounds = zeros(size(x)); 

for row = 1:numSegments
    segment = detrend(x(1+(row-1)*segmentLength:row*segmentLength)); % Fetch 160 points from recording
    m = ar(segment, modelOrder, 'Ts', Ts);
    
    % Check for unstable poles and mirror
    poles = roots(m.a);
    if max(abs(poles)) > 1
        disp('Pole outside unit circle');
        for idx = 1:size(poles,1)
            if abs(poles(idx)) > 1 
                poles(idx) = 1/poles(idx);
            end 
        end
        m.a = poly(poles);
    end 
    
    e = filter(m.a, 1, x(1+(row-1)*segmentLength:row*segmentLength));
    r = covf(e, 100);
    [A, D] = max(r(20:end));
    D = D + 19; % Look at max from t>19, so add 19 to time lag
    amp = sqrt(A);
    ehat = amp*(mod(1:160, D) == 1);
    yhat = filter(1, m.a, ehat);
    sounds(1+(row-1)*segmentLength:row*segmentLength) = yhat;
end

% --------------- Play up sound by the reconstructed sound
sounds = reshape(sounds, N, 1);
sound(10*sounds, fs);
figure;clf();
plot(t, sounds), xlabel('time [s]'); ylabel('Simulation');


