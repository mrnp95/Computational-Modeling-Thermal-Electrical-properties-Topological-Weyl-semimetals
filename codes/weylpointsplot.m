clear
%clc
N=100;

c1=.5;
c2=-c1;

Q=1;
r=.5;
kz1=linspace(-r-Q,r-Q,N);
kx1=linspace(-r,r,N);

kz2=linspace(-r+Q,r+Q,N);
kx2=linspace(-r,r,N);

[kz1,kx1]=meshgrid(kz1,kx1);
[kz2,kx2]=meshgrid(kz2,kx2);

e11=c1*(kz1+Q) - (kx1.^2 + (kz1+Q).^2).^(1/2);
e12=c1*(kz1+Q) + (kx1.^2 + (kz1+Q).^2).^(1/2);

e21=c2*(kz2-Q) - (kx2.^2 + (kz2-Q).^2).^(1/2);
e22=c2*(kz2-Q) + (kx2.^2 + (kz2-Q).^2).^(1/2);

sur=zeros(N,N);

surf(kz1,kx1,e11)
hold on
surf(kz1,kx1,e12)
hold on
surf(kz1,kx1,sur)
hold on
surf(kz2,kx2,e21)
hold on
surf(kz2,kx2,e22)
hold on
surf(kz2,kx2,sur)

%%
clear
N=100;
kz=linspace(-3.14,3.14,N);
c=1;


e1=c*kz-kz;
e2=c*kz+kz;

plot(kz,e1)
hold on
plot(kz,e2);

