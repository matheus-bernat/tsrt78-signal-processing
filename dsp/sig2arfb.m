function th=sig2arfb(y,na);
%Estimate parameters in an auto-regressive model
thf=sig2arpp(y,na);
thb=sig2arpp(y(end:-1:1),na);
th=0.5*(thf+thb);


