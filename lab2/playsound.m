function playsound(noise_type)
%
% Function that repeatedly plays a piece of music and noise via the sound
% card of the computer. To stop playing the sound press Ctrl + C
%
% Function call: playsound(noise_type)
%
%
% Inputs:   noise_type  Noise type to be played,  Valid noise types are:
%                       'white', 'chirp', 'sine', and 'multisine'.
%
% Edit: Isaac Skog, isaac.skog@liu.se, 2018-01-11
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Press Ctrl + C to stop playing the sound')

music=load('handel');
Fs=music.Fs;
signal=music.y;
N=length(signal);


% Generate noise
switch lower(noise_type)

    case 'white'
        noise=2*(rand(N,1)-0.5);
    case 'chirp'
        noise=chirp((0:4*Fs-1)./Fs,200,4-1/Fs,2000);
        noise=[noise flip(noise)]';
        noise=repmat(noise,ceil(N/length(noise)),1);
        noise=noise(1:N);
    case 'sine'
        noise=sin(2*pi*500/Fs*(0:N-1))' ;
    case 'multisine'
        noise=zeros(N,1);
        for k=1:5
            noise=noise+sin(2*pi*200*k/Fs*(0:N-1))' ;
        end
    otherwise
        error('Unknown noise type. Valid noise types are: white, chirp, sine, and multisine.')
end

% Normalize the noise so it always has the same energy
signal=sqrt((1/4)/mean(signal.^2)).*signal;
noise=sqrt((1/2)/mean(noise.^2)).*noise;


max(abs(signal))

while 1
   sound([noise signal],Fs);
   pause(1.03*N/Fs)
end