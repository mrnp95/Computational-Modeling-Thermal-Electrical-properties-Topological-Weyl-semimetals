clear 
%clc

syms a b c d e ;
syms bz;
%syms H;
%syms zso xsz yso xsx xsy osz;
syms mo m1 m2 et bt gm ;
syms kx ky kz ;
%syms kz;
%kx=0;
%ky=0;

bz=0;              % Zeeman term


%a=mo -(m1).*(kz).^2-m2.*((kx).^2+(ky).^2);
%b=et.*(kx);
%c=-et.*(ky);
%d=(bt+gm).*(kz).*((ky).^2-(kx).^2);
%e=-2.*(bt-gm).*(kz).*(kx).*(ky);

%a=1;
%b=2;
%c=3;
%d=-1;
%e=-2;
%bz=4;

zso=[1 0 0 0;0 -1 0 0;0 0 1 0;0 0 0 -1];

xsz=[0 1 0 0 ;1 0 0 0;0 0 0 -1;0 0 -1 0];

yso=[0 -1i 0 0 ;1i 0 0 0 ;0 0 0 -1i;0 0 1i 0];

xsx=[0 0 0 1;0 0 1 0;0 1 0 0 ;1 0 0 0];

xsy=[0 0 0 -1i;0 0 -1i 0;0 1i 0 0;1i 0 0 0 ];

osz=[1 0 0 0 ;0 1 0 0;0 0 -1 0;0 0 0 -1];

H=a*zso+b*xsz+c*yso+d*xsx+e*xsy+bz*osz;
[V,D]=eig(H);
dd=eig(H);

delHkx=diff(H,kx);
delHky=diff(H,ky);

%E1=dd(1);
%berkx=diff(H,kx);
%berky=diff(H,ky);
%berry1=0;

%for i=2:4
    
%berry1=(((V(:,1)'*(berkx*V(:,i)))*(V(:,i)'*(berky*V(:,1)))-(V(:,1)'*(berky*V(:,i)))*(V(:,i)'*(berkx*V(:,1))))/(dd(1)-dd(i))^2)+berry1;

%end
%berry1=(((V(:,1).'*berkx*V(:,2))*(V(:,2).'*berky*V(:,1))-(V(:,1).'*berky*V(:,2))*(V(:,2).'*berkx*V(:,1)))/(dd(1)-dd(2))^2)+berry1;
%berry1= (diff(V(:,1),kx))'*diff(V(:,1),ky)-(diff(V(:,1),ky))'*diff(V(:,1),kx);


%berry1=berry1*i
%delekx=diff(E1,kx);
%del2ekx2=diff(delekx,kx);

%deleky=diff(E1,ky);
%del2eky2=diff(deleky,ky);

%del2ekykx=diff(deleky,kx);




