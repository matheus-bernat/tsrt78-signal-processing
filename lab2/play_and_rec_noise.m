function [A_rec,B_rec]=play_and_rec_noise(port,noise_type)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that plays noise via the sound card of the computer and records
% the noise via the Arduino board. For the script to function properly
% the program "DSPLAB_part1" should be loaded into the Arduino board and
% the board connected to the computer via its USB-port. The name of the
% communaction port used can be found under the meny "tools" in the Arduino
% IDE program or in the list of COM-ports in the device manager.
%
% Sampling frequency: 8000 Hz.
%
% OBS! The noise is only played in one of the sound channels.
%
% Function call: [A_rec,B_rec]=play_and_rec_noise(port,noise_type)
%
%
% Inputs:   port        communication port, e.g., 'COM11'. The port to be
%                       used can be found via the Arduino IDE.
%           noise_type  Noise type to be played,  Valid noise types are:
%                       'white', 'chirp', 'sine', and 'multisine'.
%
% Output:   A_rec       16000x1 vector with the sound recorded via
%                       microphone A.
%
%           B_rec       16000x1 vector with the sound recorded via
%                       microphone B.
%
%
%
% Edit: Isaac Skog, isaac.skog@liu.se, 2018-11-14
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of samples to be recorded (Arduino sampling freq. = 8 kHz).
% To change this value, the Arduino code most be changed as well.
Nrec=16000;

% Sampling frequency matlab
Fs=44.1e3;

% Length in seconds of noise in matlab
T=5;

% Number of samples in matlab
N=T*Fs;


% Generate noise
switch lower(noise_type)

    case 'white'
        noise=2*(rand(N,1)-0.5);
    case 'chirp'
        noise=chirp((0:Fs-1)./Fs,200,1-1/Fs,2000);
        noise=[noise flip(noise)]';
        noise=repmat(noise,ceil(N/length(noise)),1);
        noise=noise(1:N);
    case 'sine'
        noise=sin(2*pi*400/Fs*(0:N-1))' ;
    case 'multisine'
        noise=zeros(N,1);
        for k=1:20
            noise=noise+sin(2*pi*200*k/Fs*(0:N-1))' ;
        end
    otherwise
        error('Unknown noise type. Valid noise types are: white, chirp, sine, and multisine.')
end

% Normalize the noise so it always has the same energy
noise=sqrt((1/3)/mean(noise.^2)).*noise;


% open serial port
com = serial(upper(port),'BaudRate', 115200,'InputBufferSize',100000,'ByteOrder','bigEndian');
fopen(com);
pause(1)

% Flush buffers
while(com.BytesAvailable>0)
    fread(com);
end

% Start playing sound
sound([noise zeros(size(noise))],Fs');

% To avoid tranisient effects of sound card etc.
pause(1);

% Send command to start recording data
fwrite(com,255,'uint8');

% Read recorded data
A_rec=zeros(Nrec,1);
B_rec=zeros(Nrec,1);
ctr=0;
while ctr<Nrec
    if com.BytesAvailable>3
        ctr=ctr+1;
        B_rec(ctr)=fread(com,1,'int16');
        A_rec(ctr)=fread(com,1,'int16');
    end
end

% Clean thing up
fclose(com);
delete(com);
clear com


% Plot the signals
matlab_version=version('-release');
t=(0:Nrec-1)./8000;
figure(1001)
clf
subplot(2,1,1);
plot(t,A_rec);
hold on;
plot([t(1) t(end)],[2^15-1 2^15-1],'r--');
if str2double(matlab_version(1:4))>2016
    legend('Signal','Saturation level','AutoUpdate','off')
else
    legend('Signal','Saturation level')
end
plot([t(1) t(end)],[-2^15 -2^15],'r--');
title('A-signal')
xlabel('Time [s]')
ylabel('Amplitude')

subplot(2,1,2)
plot(t,B_rec);
hold on;
plot([t(1) t(end)],[2^15-1 2^15-1],'r--');
if str2double(matlab_version(1:4))>2016
    legend('Signal','Saturation level','AutoUpdate','off')
else
    legend('Signal','Saturation level')
end
plot([t(1) t(end)],[-2^15 -2^15],'r--');
title('B-signal')
xlabel('Time [s]')
ylabel('Amplitude')
end