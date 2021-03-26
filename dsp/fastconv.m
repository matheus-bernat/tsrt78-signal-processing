function z=fastconv(x,y)
%Fast convolution in the frequency domain
   Nx=length(x);
   Ny=length(y);
   xn = [x(:); zeros(Ny,1)];
   yn = [y(:); zeros(Nx,1)];
   z = real(ifft(fft(yn).*fft(xn)));
   z(end)=[];  % Last value is always 0
