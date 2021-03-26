% This file comes with the file recording.wav
% Tested on a Laptop with Matlab R2014a, Windows 7, November 2014. 

% The code shows how to read in a wav file that you have recorded, for example,
% on your Android phone. 

%% first, we'll read in the pre-recorded wave file

% We assume that the wave file has been recorded with a sampling frequency
% of 8000Hz! 

[y, fSamp] = audioread('whistle.wav'); % extract data and sampling frequency

fSamp % check if 8000Hz
size(y) % if your recording is in stereo, you can remove one of the channels

% Some advice:
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