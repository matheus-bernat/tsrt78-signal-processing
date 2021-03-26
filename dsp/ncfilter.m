function y=ncfilter(b,a,u)
%Non-causal filter
ind=find(b~=0);
k1=b(ind(1));
z=roots(b);
ind=find(abs(z)<=1);
bf=poly(z(ind));
ind=find(abs(z)>1);
bb=poly(1./z(ind));
k2=bb(end);

z=roots(a);
ind=find(abs(z)<=1);
af=poly(z(ind));
ind=find(abs(z)>1);
ab=poly(1./z(ind));
k3=ab(end);
y=k1/k2*k3*fbfilter(bf,af,bb,ab,u);