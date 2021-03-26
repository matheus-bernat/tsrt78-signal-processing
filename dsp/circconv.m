function z=circconv(x,y)
%Circular convolution
   z = real(ifft(fft(y(:)).*fft(x(:))));
