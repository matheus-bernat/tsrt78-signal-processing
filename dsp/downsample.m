function x=downsample(x,D)
%Downsampling a factor D
   N=length(x);
   if N/D ~= round(N/D), error('N/D not an integer'), end
   X=fft(x(:));
   X=[X(1:N/D/2+1);X(end-N/D/2+2:end)];
   x=real(ifft(X));

