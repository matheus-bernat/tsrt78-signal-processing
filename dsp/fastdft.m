function X=fastdft(x)
x=x(:);
N=length(x);
if log2(N)~=round(log2(N)), 
    error('N must be a power of two'), 
end
if N==2;
    X=[x(1)+x(2);x(1)-x(2)];
else
   Xe=fastdft(x(1:2:end));
   Xo=fastdft(x(2:2:end));
   X=[Xe;Xe]+exp(-i*2*pi/N*(0:N-1)').*[Xo;Xo];
end
