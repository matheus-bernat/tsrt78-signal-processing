function y=freqfilter(b,a,u)
%Frequency domain filter
N=length(u);
omega = 2*pi*(0:N-1)/N;
H = freqz(b,a,omega);
U = fft(u)
Y = U.* H;
y = ifft(Y);

