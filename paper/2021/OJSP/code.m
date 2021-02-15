pad=’symmetric’;
N=(2*r+1)^2;
h=ones(2*r+1)/N;

%patch mean of I

mu=imfilter(I,h,pad);

%patch mean of G

nu=imfilter(G,h,pad);

%patch cov

phi=imfilter(I.*G,h,pad)-mu.*nu;

%patch var of G
vS=imfilter(G.*G,h,pad)-nu.*nu;
a=phi./(vS+Epsilon);
Beta=(a+sign(phi).*sqrt(a.^2+4*kappa...
	*Epsilon./(vS+Epsilon)))/2;

%weight calculation

w=vS./(s*mean(vS(:)));
w=1./(1+w.^2);
nor=imfilter(w,h,pad);

%final output

A=imfilter(Beta.*w,h,pad);

B=imfilter((mu-Beta.*nu).*w,h,pad);J=(G.*A+B)./nor;