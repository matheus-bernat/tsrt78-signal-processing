function y=fbfilter(bf,af,bb,ab,u)
%Forward-backward filtering
y1 = filter(bf,af,u);
yr = filter(bb,ab,y1(end:-1:1));
y = yr(end:-1:1);
