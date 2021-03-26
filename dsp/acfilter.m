function y=acfilter(b,a,u)
%Anti-causal filter
yr = filter(b(end:-1:1),a(end:-1:1),u(end:-1:1));
y  = yr(N:-1:1);
