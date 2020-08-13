clear
%clc


syms gm ko m tx t;
syms kx ky kz;


so=[1 0;0 1];
sx=[0 1;1 0];
sy=[0 -1i;1i 0];
sz=[1 0;0 -1];

ko=0;
tx=t;
m=2*t;


H=gm*(cos(kz)-sin(ko))*so-(m*(2-cos(ky)-cos(-kx))+2*tx*(cos(kz)-sin(ko)))*sx-2*t*sin(ky)*sy-2*t*sin(-kx)*sz;

[V,D]=eig(H);

dd=eig(H);


delHkx=diff(H,kx);
delHky=diff(H,ky);


