% This file comes with the simulink diagram recordSounds.mdl
% It requires a functioning microphone and the DSP systems toolbox.
% Tested on a Laptop with Matlab R2014a, Windows 7, October 2014. 

% The code records a sound file using simulink.
% Please run this file section by section!

%% first, we'll open simulink and record some sounds

tRec = 10; % will record tRec seconds
fSamp = 8000; % sampling frequency in Hz

open_system('recordSounds'); % opens the simulink diagram

% Some advice:
% - You can open the "From Audio Device" and adjust the recording device if
%   necessary. 
% - Start the recording with the play button.
% - After running the simulink diagram, the variable y will appear in your 
%   workspace. It contains the recorded sound signal.
% - Make sure that the volume of your recording device is reasonable. You
%   can listen to your recordings in matlab and also look at a plot of the
%   data to check the quality.

%% second, we'll listen to our recording

sound(y,fSamp); % make sure that the quality of the recording is okay

%% third, we'll have a look at the data

nSamp = size(y,1); % number of samples
t = (0:nSamp-1)/fSamp; % time vector in seconds

figure(1);clf();
plot(t, y)
xlabel('time in seconds')
ylabel('recorded signal') % axis description is important!