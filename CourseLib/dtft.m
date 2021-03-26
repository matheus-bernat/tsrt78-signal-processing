  function [X,w] = dtft(x,Ts,option)

%DTFT  Discrete Time Fourier Transform of truncated signals
%
%   Computes the DTFT of signals
%
%   SYNTAX:
%           [X,w] = dtft(x)
%
%   Sample interval can be provided:
%      [X,w] = dtft(x,Ts)
%

% Fredrik Gunnarsson, Linköpings universitet, 990201

if nargin < 1,
  error('No signal provided.')
end

if nargin < 2,
  Ts = 1;
end

N = length(x);
if nargin < 3,
  p=16;
  option='end';
  if N*p < 1024,
    p=ceil(1024/N);
  end
else
  p=16;
  if N*p < 4096,
    p=ceil(4096/N);
  end
end

if size(x,1) > size(x,2),
  tall = 1;
  x=x';
else
  tall = 0;
end


if strcmp(option,'middle'),

  if mod(N,2),
    x1 = [x(1+(0:(N-1)/2)),...
          zeros(1,N*(p-1)),...
          x(1+((N+1)/2:N-1))];
  else
    x1 = [x(1+(0:(N/2-1))),...  
          zeros(1,N*(p-1)/2),...
          x(1+N/2),...          
          zeros(1,N*(p-1)/2),...
          x(1+((N/2+1):N-1))]; 
  end

else
  x1 = [x,zeros(1,N*(p-1))];
end

w = 2*pi/N/p/Ts*(0:N*p-1);

if tall,
  x1=x1';
  w=w';
end

X = Ts*fft(x1);

if (nargout==0),
  subplot(2,1,1)
  plot(w,abs(X));
  xlabel('rad/s');
  title('DTFT');
  set(gca,'xlim',[min(w),max(w)]);
  ylabel('Magnitude');
  subplot(2,1,2)
  plot(w,angle(X));
  xlabel('rad/s');
  ylabel('Phase');
  set(gca,'xlim',[min(w),max(w)]);
end
