close all; 
clear;
clc;

x=linspace(-200,20,1000);
y=linspace(0,10, 100);

[xx,yy]=meshgrid(x,y);
rr=sqrt(xx.^2+yy.^2);

Ve=20000;
Re=250e-6;
a=1;
b=1;
zz= Ve ./ ((1+a.*(rr-Re).^2).*(1+b.*(rr-Re).^2));

%f=@(r) Ve ./ ((1+a.(r-Re).^2).*((1+b.*(r-Re).^2)));
%rho=@(x,y) sqrt(x.^2+y.^2);

%rr=rho(xx,yy);
%zz=f(rr);

plot=surf(xx,yy,zz);
axis([-200,20,0,10]);


%%
close all; 
clear;
clc;

x=linspace(-2,0.20,1000);
y=linspace(0,0.10, 100);

[xx,yy]=meshgrid(x,y);
rr=sqrt(xx.^2+yy.^2);

Ve=20000;
Re=250e-6;
a=100000;
zz= Ve ./ (1+a^2.*(rr-Re).^2);

plot=surf(xx,yy,zz);
axis([-2,0.20,0,0.10]);

%% Grad exact
clc
x=256e-6;
y= 400e-6;

r = sqrt(x * x + y * y);

Ve = 20000;
Re = 250e-6;
a = 100;
dfdr = - Ve ./ ( 1 + ( (a.*(r - Re)).*(a*(r - Re)) ).*( (a.*(r - Re)).*(a*(r - Re)) ) ) .*  a.*a*2*(r-Re) .* x ./ r;
grad_x = dfdr .* x ./ r
grad_y= dfdr .* y ./ r


%%
close all; 
clear;
clc;

x=linspace(-2,0.20,1000);
y=linspace(0,0.10, 100);

[xx,yy]=meshgrid(x,y);
rr=sqrt(xx.^2+yy.^2);

Ve=20000;
Re=250e-6;

zz= Ve./ (2-rr./Re) .* (rr<2*Re);

plot=surf(xx,yy,zz);
axis([-2,0.20,0,0.10]);