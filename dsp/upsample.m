function x=upsample(x,U)
%Upsampling a factor U
   if U ~= round(U), error('U is not an integer'), end
   N=length(x);
   X=fft(x(:));
   X=[X(1:N/2+1);zeros((U-1)*N,1);X(N/2+2:N)];
   x=real(ifft(X));