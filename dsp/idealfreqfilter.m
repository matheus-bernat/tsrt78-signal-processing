function y=idealfreqfilter(wc,u,type)
%Ideal frequency domain filter
% wc=[w1 w2], w1=0 and w2=pi are allowed
% type='stop' or 'pass' (default)
if nargin<3, type='pass'; end
if length(wc)==1; wc=[wc pi]; end
N=length(u);
U=fft(u);
w=2*pi*(0:N/2)/N;
if type=='pass'
   ind=find(wc(1)<=w & wc(2)>=w);
else
   ind=find(wc(1)>=w | wc(2)<=w);
end
Y=zeros(N/2+1,1);
Y(ind)=U(ind);
Y=[Y;conj(Y(N/2:-1:2))];
y=real(ifft(Y));
