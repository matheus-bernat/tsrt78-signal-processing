function y=zpfilter(b,a,u)
%Zero-phase filtering ('filtfilt')
yr = filter(b,a,u(end:-1:1));
y  = filter(b,a,yr(end:-1:1));
