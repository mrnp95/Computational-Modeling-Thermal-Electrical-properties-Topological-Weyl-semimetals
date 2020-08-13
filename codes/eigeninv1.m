clear
%clc


syms gm ko m tx t  bz;
syms kx ky kz;


so=[1 0;0 1];
sx=[0 1;1 0];
sy=[0 -1i;1i 0];
sz=[1 0;0 -1];

kzo=0;
kxo=0;
tx=t/2;
m=2*t;

%original hamiltonian
H=-(m*(1-(cos(kz))^2-cos(ky))+2*tx*(cos(kx)-sin(kxo)))*sx-2*t*sin(ky)*sy-2*t*cos(kz)*sz+bz*sz ;

%{
%H of watzman paper
H=-(m*(1-(cos(ky))^2-cos(kz))+2*tx*(cos(kx)-sin(kxo)))*sx-2*t*cos(ky)*sy-2*t*sin(kz)*sz;
%}

[V,D]=eig(H);

dd=eig(H);


delHkx=diff(H,kx);
delHky=diff(H,ky);


