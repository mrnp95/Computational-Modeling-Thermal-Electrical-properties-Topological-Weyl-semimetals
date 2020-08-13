clear
%clc


syms gm ko m tx t;
syms kx ky kz;


so=[1 0;0 1];
sx=[0 1;1 0];
sy=[0 -1i;1i 0];
sz=[1 0;0 -1];

kzo=0;
kxo=0;
tx=t/2;
m=2*t;


H=gm*(cos(2*kz)-sin(kzo))*(cos(kx)-sin(kxo))*so-(m*(1-(cos(kx))^2-cos(ky))+2*tx*(cos(kz)-sin(kzo)))*sx-2*t*sin(ky)*sy-2*t*cos(kx)*sz;

[V,D]=eig(H);

dd=eig(H);


delHkx=diff(H,kx);
delHky=diff(H,ky);


