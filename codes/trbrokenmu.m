clear 
clc


%syms sigxx alpxx sigxy alpxy sxx axx sxy1 sxy2 axy1 axy2;
%syms sigz B h sk kb;
%syms kx ky kz;
%syms vx vy;
%syms e f delf mu tau E beta T;

%Dispersion relation with mag field
%syms bz;
%syms mo m1 m2 et bt gm ;
%syms a b c d e ;

et=50.*10.^-6;
bt=-et./5;
gm=et;
mo=-2*et;
m1=-4*et;
m2=-et/5;
%bz=30*et/5;
%bz=0;

%{
%tau=10.^-14;                %%%find it's value
tau=1;

temp=100;                %temp 100K
T=temp;                      %%%onsagar coefficient find it's value
%T=1;
mu=0*et/5;         %chemical potential
charge=1;
%h=6.626*10^-34/(2*pi);
h=4.135667662*10^(-15);      %plank constant in ev
%beta=1/(1.38*10^(-23)*temp); 
beta=1/(8.6173303*10^(-5)*temp); %bolzman constant in ev
kbolz=8.61733*10^-5;

%}
con=1;
ber=0;

charge=1;
tau=1;
temp=100;                %temp 100K
T=temp;                      %%%onsagar coefficient find it's value
h=1;
kbolz=8.61733*10^-5;
beta=1/(kbolz*T);
%mu=0.16*et/5;

NZ=20;                %no of values of kz
NX=NZ;
NY=NZ;

%uppr=3.16*10^0;                   %uppr limit of kz
%lowr=(-1)*uppr;                %lowr limit of kz
upprz=3.16;
lowrz=-3.16;
upprx=3.16;
lowrx=-3.16;
uppry=3.16;
lowry=-3.16;
deldiff=abs((upprz-lowrz)/(2*2*pi*NZ)); %difference between each value of kz
ddiff=abs((upprz-lowrz)/NZ);

%ky=0;
%kx=0;
kz=linspace(lowrz,upprz,NZ);
kx=linspace(lowrx,upprx,NX);
ky=linspace(lowry,uppry,NY);
%ky=0;
%kx=0;
[kz,kx,ky]=meshgrid(kz,kx,ky);
%[kz,kx]=meshgrid(kz,kx);

KZ=linspace(lowrz,upprz,NZ-2);
KX=linspace(lowrx,upprx,NX-2);
KY=linspace(lowry,uppry,NY-2);

%[Kz,Kx,Ky]=meshgrid(Kz,KX,Ky);
[KZ,KX,KY]=meshgrid(KZ,KX,KY);


a=mo -(m1).*(kz).^2-m2.*((kx).^2+(ky).^2);
b=et.*(kx);
c=-et.*(ky);
d=(bt+gm).*(kz).*((ky).^2-(kx).^2);
e=-2.*(bt-gm).*(kz).*(kx).*(ky);


%%%

maxB=100;                             %max value of B
minB=-maxB;                              % min value of B
noofB=10;                              %no of values of B 
zman=linspace(minB,maxB,noofB);

%nernst=zeros(1,noofB);


B=50;
%for  bfield =1:noofB;
    %B=zman(bfield);
    bz=5.8*10^-5*B;
    
noofmu=10;
mumax=.2*et/5;
muran=linspace(0,mumax,noofmu);

sigxyp=zeros(1,noofmu);
alpxyp=zeros(1,noofmu);

for munum=1:noofmu

    mu=muran(munum);
%Eigenvalues and eigenvector

H=zeros(4,4,NZ,NX,NY);

H(1,1,:,:,:)=a+bz;
H(2,1,:,:,:)=b+c*1i;
H(3,1,:,:,:)=0;
H(4,1,:,:,:)=d+e*1i;

H(1,2,:,:,:)=b-c*1i;
H(2,2,:,:,:)=bz-a;
H(3,2,:,:,:)=d+e*1i;
H(4,2,:,:,:)=0;

H(1,3,:,:,:)=0;
H(2,3,:,:,:)=d-e*1i;
H(3,3,:,:,:)=a-bz;
H(4,3,:,:,:)=-b+c*1i;

H(1,4,:,:,:)=d-e*1i;
H(2,4,:,:,:)=0;
H(3,4,:,:,:)=-b-c*1i;
H(4,4,:,:,:)=-a-bz;

%E1=(a.^2 + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2 - 2.*bz.*(a.^2 + b.^2 + c.^2).^(1/2)).^(1/2);
%E2=(a.^2 + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2 + 2.*bz.*(a.^2 + b.^2 + c.^2).^(1/2)).^(1/2);
%E3=-(a.^2 + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2 - 2.*bz.*(a.^2 + b.^2 + c.^2).^(1/2)).^(1/2);
%E4=-(a.^2 + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2 + 2.*bz.*(a.^2 + b.^2 + c.^2).^(1/2)).^(1/2);

%{
E1= (a.^2 + b.^2 + c.^2 + d.^2 + e.^2).^(1/2);
E2= (a.^2 + b.^2 + c.^2 + d.^2 + e.^2).^(1/2);
E3= -(a.^2 + b.^2 + c.^2 + d.^2 + e.^2).^(1/2);
E4= -(a.^2 + b.^2 + c.^2 + d.^2 + e.^2).^(1/2);
%}

E1=zeros(NZ,NX,NY);
E2=zeros(NZ,NX,NY);
E3=zeros(NZ,NX,NY);
E4=zeros(NZ,NX,NY);

efun1=zeros(4,NZ,NX,NY);
efun2=zeros(4,NZ,NX,NY);
efun3=zeros(4,NZ,NX,NY);
efun4=zeros(4,NZ,NX,NY);

for z=1:NZ
    for x=1:NX
        for y=1:NY
            %V=zeros(4,4);
            %D=zeros(4,4);
            [V,D]=eig(H(:,:,z,x,y));
            %dd=eig(H(:,:,z,x,y));
            E1(z,x,y)=D(1,1);
            E2(z,x,y)=D(2,2);
            E3(z,x,y)=D(3,3);
            E4(z,x,y)=D(4,4);
            
            efun1(:,z,x,y)=V(:,1);
            efun2(:,z,x,y)=V(:,2);
            efun3(:,z,x,y)=V(:,3);
            efun4(:,z,x,y)=V(:,4);
            
        end
    end
end

%%%
%Differentiation of Hamiltonian
delHkx=zeros(4,4,NZ,NX,NY);
delHky=zeros(4,4,NZ,NX,NY);

%delHkx
delHkx(1,1,:,:,:)=-2.*kx.*m2;
delHkx(2,1,:,:,:)=et;
delHkx(3,1,:,:,:)=0;
delHkx(4,1,:,:,:)=-ky.*kz.*(2*bt - 2*gm)*1i - 2.*kx.*kz.*(bt + gm);
            
delHkx(1,2,:,:,:)=et;
delHkx(2,2,:,:,:)=2.*kx.*m2;
delHkx(3,2,:,:,:)=-ky.*kz.*(2*bt - 2*gm)*1i - 2.*kx.*kz*(bt + gm);
delHkx(4,2,:,:,:)=0;
            
delHkx(1,3,:,:,:)=0;
delHkx(2,3,:,:,:)=ky.*kz.*(2*bt - 2*gm)*1i - 2*kx.*kz.*(bt + gm);
delHkx(3,3,:,:,:)=-2.*kx.*m2;
delHkx(4,3,:,:,:)=-et;
            
delHkx(1,4,:,:,:)=ky.*kz.*(2*bt - 2*gm)*1i - 2.*kx.*kz.*(bt + gm);
delHkx(2,4,:,:,:)=0;
delHkx(3,4,:,:,:)=-et;
delHkx(4,4,:,:,:)=2.*kx.*m2;
    

%delHky
delHky(1,1,:,:,:)=-2.*ky.*m2;
delHky(2,1,:,:,:)=-et*1i;
delHky(3,1,:,:,:)=0;
delHky(4,1,:,:,:)=2.*ky.*kz.*(bt + gm) - kx.*kz.*(2*bt - 2*gm)*1i;

delHky(1,2,:,:,:)=et*1i;
delHky(2,2,:,:,:)=2.*ky.*m2;
delHky(3,2,:,:,:)=2.*ky.*kz.*(bt + gm) - kx.*kz.*(2*bt - 2*gm)*1i;
delHky(4,2,:,:,:)=0;

delHky(1,3,:,:,:)=0;
delHky(2,3,:,:,:)=kx.*kz.*(2*bt - 2*gm)*1i + 2.*ky.*kz.*(bt + gm);
delHky(3,3,:,:,:)=-2.*ky.*m2;
delHky(4,3,:,:,:)=-et*1i;

delHky(1,4,:,:,:)=kx.*kz.*(2*bt - 2*gm)*1i + 2.*ky.*kz.*(bt + gm);
delHky(2,4,:,:,:)=0;
delHky(3,4,:,:,:)=et*1i;
delHky(4,4,:,:,:)=2.*ky.*m2;



%fermi distribution
f1=1./(1+exp((E1-mu)*beta));
mf1=exp(E1-mu)*beta./(1+exp((E1-mu)*beta));
df1dE1=-(beta*exp(beta*(E1 - mu)))./((exp(beta*(E1 - mu)) + 1).^2);    %actually delf/delE
%sf1=-(1./(1+exp((E1-mu)*beta))).*log(1./(1+exp((E1-mu)*beta)))-(exp(E1-mu)*beta./(1+exp((E1-mu)*beta))).*log(exp(E1-mu)*beta./(1+exp((E1-mu)*beta)));
sf1=-f1.*log(f1)-(1-f1).*log(1-f1);
%sf1=-f1.*log(f1)-mf1.*log(mf1);

f2=1./(1+exp((E2-mu)*beta));
mf2=exp(E2-mu)*beta./(1+exp((E2-mu)*beta));
df2dE2=-(beta*exp(beta*(E2 - mu)))./((exp(beta*(E2 - mu)) + 1).^2);    %actually delf/delE
%sf2=-(1./(1+exp((E2-mu)*beta))).*log(1./(1+exp((E2-mu)*beta)))-(exp(E2-mu)*beta./(1+exp((E2-mu)*beta))).*log(exp(E2-mu)*beta./(1+exp((E2-mu)*beta)));
sf2=-f2.*log(f2)-(1-f2).*log(1-f2);
%sf2=-f2.*log(f2)-mf2.*log(mf2);

f3=1./(1+exp((E3-mu)*beta));
mf3=exp(E3-mu)*beta./(1+exp((E3-mu)*beta));
df3dE3=-(beta*exp(beta*(E3 - mu)))./((exp(beta*(E3 - mu)) + 1).^2);    %actually delf/delE
%sf3=-(1./(1+exp((E3-mu)*beta))).*log(1./(1+exp((E3-mu)*beta)))-(exp(E3-mu)*beta./(1+exp((E3-mu)*beta))).*log(exp(E3-mu)*beta./(1+exp((E3-mu)*beta)));
sf3=-f3.*log(f3)-(1-f3).*log(1-f3);
%sf3=-f3.*log(f3)-mf3.*log(mf3);

f4=1./(1+exp((E4-mu)*beta));
mf4=exp(E4-mu)*beta./(1+exp((E4-mu)*beta));
df4dE4=-(beta*exp(beta*(E4 - mu)))./((exp(beta*(E4 - mu)) + 1).^2);    %actually delf/delE
%sf4=-(1./(1+exp((E4-mu)*beta))).*log(1./(1+exp((E4-mu)*beta)))-(exp(E4-mu)*beta./(1+exp((E4-mu)*beta))).*log(exp(E4-mu)*beta./(1+exp((E4-mu)*beta)));
sf4=-f4.*log(f4)-(1-f4).*log(1-f4);
%sf4=-f4.*log(f4)-mf4.*log(mf4);



%{

%V11=(a.^2 + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2 - 2.*bz.*(a.^2 + b.^2 + c.^2).^(1./2)).^(1./2)./(d + e.*1i) - (a.^2 + b.^2 - bz.^2 + c.^2 + d.^2 + e.^2)./(2.*bz.*(d + e.*1i)) + (a.^2 + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2 - 2.*bz.*(a.^2 + b.^2 + c.^2).^(1./2))./(2.*bz.*(d + e.*1i));
%V12=(a.^2 + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2 - 2.*bz.*(a.^2 + b.^2 + c.^2).^(1./2)).^(3./2)./(2.*(d + e.*1i).*(b.*bz - bz.*c.*1i)) + (a.^3 + a.^2.*bz + a.*b.^2 - a.*bz.^2 + a.*c.^2 + a.*d.^2 + a.*e.^2 + b.^2.*bz - bz.^3 + bz.*c.^2 - bz.*d.^2 - bz.*e.^2)./(2.*(d + e.*1i).*(b.*bz - bz.*c.*1i)) - ((a.^2 + 2.*a.*bz + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2).*(a.^2 + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2 - 2.*bz.*(a.^2 + b.^2 + c.^2).^(1./2)).^(1./2))./(2.*(d + e.*1i).*(b.*bz - bz.*c.*1i)) - ((a - bz).*(a.^2 + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2 - 2.*bz.*(a.^2 + b.^2 + c.^2).^(1./2)))./(2.*(d + e.*1i).*(b.*bz - bz.*c.*1i));
%V13=(a.^2 + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2 - 2.*bz.*(a.^2 + b.^2 + c.^2).^(1./2))./(2.*(b.*bz - bz.*c.*1i)) - (a.^2 + 2.*a.*bz + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2)./(2.*(b.*bz - bz.*c.*1i));
%V14=1;

V11=(b - c.*1i)./(d + e.*1i);
V12=- a./(d + e.*1i) + (a.^2 + b.^2 + c.^2 + d.^2 + e.^2).^(1/2)./(d + e.*1i);
V13=1;
V14=0;

efun1=zeros(4,NZ,NX,NY);

         efun1(1,:,:,:)=V11;   
         efun1(2,:,:,:)=V12;
         efun1(3,:,:,:)=V13;
         efun1(4,:,:,:)=V14;
            
        
%V21=(a.^2 + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2 + 2.*bz.*(a.^2 + b.^2 + c.^2).^(1./2)).^(1./2)./(d + e.*1i) - (a.^2 + b.^2 - bz.^2 + c.^2 + d.^2 + e.^2)./(2.*bz.*(d + e.*1i)) + (a.^2 + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2 + 2.*bz.*(a.^2 + b.^2 + c.^2).^(1./2))./(2.*bz.*(d + e.*1i));
%V22=(a.^2 + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2 + 2.*bz.*(a.^2 + b.^2 + c.^2).^(1./2)).^(3./2)./(2.*(d + e.*1i).*(b.*bz - bz.*c.*1i)) + (a.^3 + a.^2.*bz + a.*b.^2 - a.*bz.^2 + a.*c.^2 + a.*d.^2 + a.*e.^2 + b.^2.*bz - bz.^3 + bz.*c.^2 - bz.*d.^2 - bz.*e.^2)./(2.*(d + e.*1i).*(b.*bz - bz.*c.*1i)) - ((a.^2 + 2.*a.*bz + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2).*(a.^2 + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2 + 2.*bz.*(a.^2 + b.^2 + c.^2).^(1./2)).^(1./2))./(2.*(d + e.*1i).*(b.*bz - bz.*c.*1i)) - ((a - bz).*(a.^2 + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2 + 2.*bz.*(a.^2 + b.^2 + c.^2).^(1./2)))./(2.*(d + e.*1i).*(b.*bz - bz.*c.*1i));
%V23= (a.^2 + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2 + 2.*bz.*(a.^2 + b.^2 + c.^2).^(1./2))./(2.*(b.*bz - bz.*c.*1i)) - (a.^2 + 2.*a.*bz + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2)./(2.*(b.*bz - bz.*c.*1i));
%V24=1;

V21= a./(d + e.*1i) + (a.^2 + b.^2 + c.^2 + d.^2 + e.^2).^(1/2)./(d + e.*1i);
V22=(b + c.*1i)./(d + e.*1i);
V23=0;
V24=1;

efun2=zeros(4,NZ,NX,NY);

         efun2(1,:,:,:)=V21;   
         efun2(2,:,:,:)=V22;
         efun2(3,:,:,:)=V23;
         efun2(4,:,:,:)=V24;
            
         


%V31=- (a.^2 + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2 - 2.*bz.*(a.^2 + b.^2 + c.^2).^(1./2)).^(1./2)./(d + e.*1i) - (a.^2 + b.^2 - bz.^2 + c.^2 + d.^2 + e.^2)./(2.*bz.*(d + e.*1i)) + (a.^2 + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2 - 2.*bz.*(a.^2 + b.^2 + c.^2).^(1./2))./(2.*bz.*(d + e.*1i));
%V32=- (a.^2 + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2 - 2.*bz.*(a.^2 + b.^2 + c.^2).^(1./2)).^(3./2)./(2.*(d + e.*1i).*(b.*bz - bz.*c.*1i)) + (a.^3 + a.^2.*bz + a.*b.^2 - a.*bz.^2 + a.*c.^2 + a.*d.^2 + a.*e.^2 + b.^2.*bz - bz.^3 + bz.*c.^2 - bz.*d.^2 - bz.*e.^2)./(2.*(d + e.*1i).*(b.*bz - bz.*c.*1i)) + ((a.^2 + 2.*a.*bz + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2).*(a.^2 + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2 - 2.*bz.*(a.^2 + b.^2 + c.^2).^(1./2)).^(1./2))./(2.*(d + e.*1i).*(b.*bz - bz.*c.*1i)) - ((a - bz).*(a.^2 + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2 - 2.*bz.*(a.^2 + b.^2 + c.^2).^(1./2)))./(2.*(d + e.*1i).*(b.*bz - bz.*c.*1i));
%V33=(a.^2 + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2 - 2.*bz.*(a.^2 + b.^2 + c.^2).^(1./2))./(2.*(b.*bz - bz.*c.*1i)) - (a.^2 + 2.*a.*bz + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2)./(2.*(b.*bz - bz.*c.*1i));
%V34=1;

V31=(b - c.*1i)./(d + e.*1i);
V32=- a./(d + e.*1i) - (a.^2 + b.^2 + c.^2 + d.^2 + e.^2).^(1/2)./(d + e.*1i);
V33=1;
V34=0;

efun3=zeros(4,NZ,NX,NY);

         efun3(1,:,:,:)=V31;   
         efun3(2,:,:,:)=V32;
         efun3(3,:,:,:)=V33;
         efun3(4,:,:,:)=V34;
            
         

%V41=- (a.^2 + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2 + 2.*bz.*(a.^2 + b.^2 + c.^2).^(1./2)).^(1./2)./(d + e.*1i) - (a.^2 + b.^2 - bz.^2 + c.^2 + d.^2 + e.^2)./(2.*bz.*(d + e.*1i)) + (a.^2 + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2 + 2.*bz.*(a.^2 + b.^2 + c.^2).^(1./2))./(2.*bz.*(d + e.*1i));
%V42=- (a.^2 + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2 + 2.*bz.*(a.^2 + b.^2 + c.^2).^(1./2)).^(3./2)./(2.*(d + e.*1i).*(b.*bz - bz.*c.*1i)) + (a.^3 + a.^2.*bz + a.*b.^2 - a.*bz.^2 + a.*c.^2 + a.*d.^2 + a.*e.^2 + b.^2.*bz - bz.^3 + bz.*c.^2 - bz.*d.^2 - bz.*e.^2)./(2.*(d + e.*1i).*(b.*bz - bz.*c.*1i)) + ((a.^2 + 2.*a.*bz + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2).*(a.^2 + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2 + 2.*bz.*(a.^2 + b.^2 + c.^2).^(1./2)).^(1./2))./(2.*(d + e.*1i).*(b.*bz - bz.*c.*1i)) - ((a - bz).*(a.^2 + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2 + 2.*bz.*(a.^2 + b.^2 + c.^2).^(1./2)))./(2.*(d + e.*1i).*(b.*bz - bz.*c.*1i));
%V43=(a.^2 + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2 + 2.*bz.*(a.^2 + b.^2 + c.^2).^(1./2))./(2.*(b.*bz - bz.*c.*1i)) - (a.^2 + 2.*a.*bz + b.^2 + bz.^2 + c.^2 + d.^2 + e.^2)./(2.*(b.*bz - bz.*c.*1i));
%V44=1;

V41=a./(d + e.*1i) - (a.^2 + b.^2 + c.^2 + d.^2 + e.^2).^(1/2)./(d + e.*1i);
V42=(b + c.*1i)./(d + e.*1i);
V43=0;
V44=1;

efun4=zeros(4,NZ,NX,NY);

         efun4(1,:,:,:)=V41;   
         efun4(2,:,:,:)=V42;
         efun4(3,:,:,:)=V43;
         efun4(4,:,:,:)=V44;
 %}           
 



% new berry phase
berry12=zeros(NZ-2,NX-2,NY-2);
berry13=zeros(NZ-2,NX-2,NY-2);
berry14=zeros(NZ-2,NX-2,NY-2);

berry21=zeros(NZ-2,NX-2,NY-2);
berry23=zeros(NZ-2,NX-2,NY-2);
berry24=zeros(NZ-2,NX-2,NY-2);

berry31=zeros(NZ-2,NX-2,NY-2);
berry32=zeros(NZ-2,NX-2,NY-2);
berry34=zeros(NZ-2,NX-2,NY-2);

berry41=zeros(NZ-2,NX-2,NY-2);
berry42=zeros(NZ-2,NX-2,NY-2);
berry43=zeros(NZ-2,NX-2,NY-2);

for z=1:NZ-2
    for x=1:NX-2
        for y=1:NY-2
            berry12(z,x,y)=((efun1(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun2(:,z+1,x+1,y+1))*(efun2(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun1(:,z+1,x+1,y+1))-(efun1(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun2(:,z+1,x+1,y+1))*(efun2(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun1(:,z+1,x+1,y+1)))/(E1(z+1,x+1,y+1)-E2(z+1,x+1,y+1))^2;
            berry13(z,x,y)=((efun1(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun3(:,z+1,x+1,y+1))*(efun3(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun1(:,z+1,x+1,y+1))-(efun1(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun3(:,z+1,x+1,y+1))*(efun3(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun1(:,z+1,x+1,y+1)))/(E1(z+1,x+1,y+1)-E3(z+1,x+1,y+1))^2;
            berry14(z,x,y)=((efun1(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun4(:,z+1,x+1,y+1))*(efun4(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun1(:,z+1,x+1,y+1))-(efun1(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun4(:,z+1,x+1,y+1))*(efun4(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun1(:,z+1,x+1,y+1)))/(E1(z+1,x+1,y+1)-E4(z+1,x+1,y+1))^2;
            
            %berry13(z,x,y)=((efun1(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun3(:,z+1,x+1,y+1))*(efun3(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun1(:,z+1,x+1,y+1))-(efun1(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun3(:,z+1,x+1,y+1))*(efun3(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun1(:,z+1,x+1,y+1)));
            
            
            berry21(z,x,y)=((efun2(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun1(:,z+1,x+1,y+1))*(efun1(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun2(:,z+1,x+1,y+1))-(efun2(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun1(:,z+1,x+1,y+1))*(efun1(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun2(:,z+1,x+1,y+1)))/(E2(z+1,x+1,y+1)-E1(z+1,x+1,y+1))^2;
            berry23(z,x,y)=((efun2(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun3(:,z+1,x+1,y+1))*(efun3(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun2(:,z+1,x+1,y+1))-(efun2(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun3(:,z+1,x+1,y+1))*(efun3(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun2(:,z+1,x+1,y+1)))/(E2(z+1,x+1,y+1)-E3(z+1,x+1,y+1))^2;
            berry24(z,x,y)=((efun2(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun4(:,z+1,x+1,y+1))*(efun4(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun2(:,z+1,x+1,y+1))-(efun2(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun4(:,z+1,x+1,y+1))*(efun4(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun2(:,z+1,x+1,y+1)))/(E2(z+1,x+1,y+1)-E4(z+1,x+1,y+1))^2;
            
            berry31(z,x,y)=((efun3(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun1(:,z+1,x+1,y+1))*(efun1(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun3(:,z+1,x+1,y+1))-(efun3(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun1(:,z+1,x+1,y+1))*(efun1(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun3(:,z+1,x+1,y+1)))/(E3(z+1,x+1,y+1)-E1(z+1,x+1,y+1))^2;
            berry32(z,x,y)=((efun3(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun2(:,z+1,x+1,y+1))*(efun2(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun3(:,z+1,x+1,y+1))-(efun3(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun2(:,z+1,x+1,y+1))*(efun2(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun3(:,z+1,x+1,y+1)))/(E3(z+1,x+1,y+1)-E2(z+1,x+1,y+1))^2;
            berry34(z,x,y)=((efun3(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun4(:,z+1,x+1,y+1))*(efun4(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun3(:,z+1,x+1,y+1))-(efun3(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun4(:,z+1,x+1,y+1))*(efun4(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun3(:,z+1,x+1,y+1)))/(E3(z+1,x+1,y+1)-E4(z+1,x+1,y+1))^2;
            
            berry41(z,x,y)=((efun4(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun1(:,z+1,x+1,y+1))*(efun1(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun4(:,z+1,x+1,y+1))-(efun4(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun1(:,z+1,x+1,y+1))*(efun1(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun4(:,z+1,x+1,y+1)))/(E4(z+1,x+1,y+1)-E1(z+1,x+1,y+1))^2;
            berry42(z,x,y)=((efun4(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun2(:,z+1,x+1,y+1))*(efun2(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun4(:,z+1,x+1,y+1))-(efun4(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun2(:,z+1,x+1,y+1))*(efun2(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun4(:,z+1,x+1,y+1)))/(E4(z+1,x+1,y+1)-E2(z+1,x+1,y+1))^2;
            berry43(z,x,y)=((efun4(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun3(:,z+1,x+1,y+1))*(efun3(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun4(:,z+1,x+1,y+1))-(efun4(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun3(:,z+1,x+1,y+1))*(efun3(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun4(:,z+1,x+1,y+1)))/(E4(z+1,x+1,y+1)-E3(z+1,x+1,y+1))^2;
            
        end
    end
end


berry1=real(1i*(berry12+berry13+berry14));
berry2=real(1i*(berry21+berry23+berry24));
berry3=real(1i*(berry31+berry32+berry34));
berry4=real(1i*(berry41+berry42+berry43));

%{
berry1=real(1i*(berry13+berry14));
berry2=real(1i*(berry23+berry24));
berry3=real(1i*(berry31+berry32));
berry4=real(1i*(berry41+berry42));
%}

%berry1=real(1i*(berry13));




%{
%Berry phase 

%differentiation of wavefunction for berry phase calculation
dfun1dkx=zeros(4,NZ,NX,NY);
dfun2dkx=zeros(4,NZ,NX,NY);
dfun3dkx=zeros(4,NZ,NX,NY);
dfun4dkx=zeros(4,NZ,NX,NY);
dfun1dky=zeros(4,NZ,NX,NY);
dfun2dky=zeros(4,NZ,NX,NY);
dfun3dky=zeros(4,NZ,NX,NY);
dfun4dky=zeros(4,NZ,NX,NY);

for z=2:NZ-1
    for x=2:NX-1
        for y=2:NY-1
          dfun1dkx(:,z-1,x-1,y-1)=efun1(:,z+1,x,y)-efun1(:,z-1,x,y);
          dfun2dkx(:,z-1,x-1,y-1)=efun2(:,z+1,x,y)-efun2(:,z-1,x,y);
          dfun3dkx(:,z-1,x-1,y-1)=efun3(:,z+1,x,y)-efun3(:,z-1,x,y);
          dfun4dkx(:,z-1,x-1,y-1)=efun4(:,z+1,x,y)-efun4(:,z-1,x,y);
       end
    end
end

for z=2:NZ-1
    for y=2:NY-1
        for x=2:NX-1
          dfun1dky(:,z-1,x-1,y-1)=efun1(:,z,x,y+1)-efun1(:,z,x,y-1);
          dfun2dky(:,z-1,x-1,y-1)=efun2(:,z,x,y+1)-efun2(:,z,x,y-1);
          dfun3dky(:,z-1,x-1,y-1)=efun3(:,z,x,y+1)-efun3(:,z,x,y-1);
          dfun4dky(:,z-1,x-1,y-1)=efun4(:,z,x,y+1)-efun4(:,z,x,y-1);
         end
    end
end


dfun1dkx=dfun1dkx/(2*ddiff);
dfun1dky=dfun1dky/(2*ddiff);

dfun2dkx=dfun2dkx/(2*ddiff);
dfun2dky=dfun2dky/(2*ddiff);

dfun3dkx=dfun3dkx/(2*ddiff);
dfun3dky=dfun3dky/(2*ddiff);

dfun4dkx=dfun4dkx/(2*ddiff);
dfun4dky=dfun4dky/(2*ddiff);

%berry phase calculation
berry1=zeros(4,NZ-2,NX-2,NY-2);
for x=1:NX-2
    for y=1:NY-2
        for z=1:NZ-2
          berry1(z,x,y)=(dfun1dkx(:,z,x,y)'*dfun1dky(:,z,x,y))-(dfun1dky(:,z,x,y)'*dfun1dkx(:,z,x,y));
         end
    end
end
berry1=1i*berry1;

berry2=zeros(4,NZ-2,NX-2,NY-2);
for x=1:NX-2
    for y=1:NY-2
        for z=1:NZ-2
          berry2(z,x,y)=(dfun2dkx(:,z,x,y)'*dfun2dky(:,z,x,y))-(dfun2dky(:,z,x,y)'*dfun2dkx(:,z,x,y));
         end
    end
end
berry2=1i*berry2;

berry3=zeros(4,NZ-2,NX-2,NY-2);
for x=1:NX-2
    for y=1:NY-2
        for z=1:NZ-2
          berry3(z,x,y)=(dfun3dkx(:,z,x,y)'*dfun3dky(:,z,x,y))-(dfun3dky(:,z,x,y)'*dfun3dkx(:,z,x,y));
         end
    end
end
berry3=1i*berry3;

berry4=zeros(4,NZ-2,NX-2,NY-2);
for x=1:NX-2
    for y=1:NY-2
        for z=1:NZ-2
          berry4(z,x,y)=(dfun4dkx(:,z,x,y)'*dfun4dky(:,z,x,y))-(dfun4dky(:,z,x,y)'*dfun4dkx(:,z,x,y));
         end
    end
end
berry4=1i*berry4;


%}



%%

%velocity vx vy vxvy vyvy
vxo1=zeros(NZ-2,NX-2,NY-2);
vxo2=zeros(NZ-2,NX-2,NY-2);
vxo3=zeros(NZ-2,NX-2,NY-2);
vxo4=zeros(NZ-2,NX-2,NY-2);

for z=2:NZ-1
    for x=2:NX-1
        for y=2:NY-1
            vxo1(z-1,x-1,y-1)=E1(z+1,x,y)-E1(z-1,x,y);
            vxo2(z-1,x-1,y-1)=E2(z+1,x,y)-E2(z-1,x,y);
            vxo3(z-1,x-1,y-1)=E3(z+1,x,y)-E3(z-1,x,y);
            vxo4(z-1,x-1,y-1)=E4(z+1,x,y)-E4(z-1,x,y);
        end
    end
end
vxo1=vxo1/(2*ddiff);
vxo2=vxo2/(2*ddiff);
vxo3=vxo3/(2*ddiff);
vxo4=vxo4/(2*ddiff);


vyo1=zeros(NZ-2,NX-2,NY-2);
vyo2=zeros(NZ-2,NX-2,NY-2);
vyo3=zeros(NZ-2,NX-2,NY-2);
vyo4=zeros(NZ-2,NX-2,NY-2);
for z=2:NZ-1
    for y=2:NY-1
        for x=2:NX-1
            vyo1(z-1,x-1,y-1)=E1(z,x,y+1)-E1(z,x,y-1);
            vyo2(z-1,x-1,y-1)=E2(z,x,y+1)-E2(z,x,y-1);
            vyo3(z-1,x-1,y-1)=E3(z,x,y+1)-E3(z,x,y-1);
            vyo4(z-1,x-1,y-1)=E4(z,x,y+1)-E4(z,x,y-1);
        end
    end
end
vyo1=vyo1/(2*ddiff);
vyo2=vyo2/(2*ddiff);
vyo3=vyo3/(2*ddiff);
vyo4=vyo4/(2*ddiff);

vyovyo1=zeros(NZ-2,NX-2,NY-2);
vyovyo2=zeros(NZ-2,NX-2,NY-2);
vyovyo3=zeros(NZ-2,NX-2,NY-2);
vyovyo4=zeros(NZ-2,NX-2,NY-2);
for z=2:NZ-1
    for y=2:NY-1
        for x=2:NX-1
            vyovyo1(z-1,x-1,y-1)=E1(z,x,y+1)+E1(z,x,y-1)-2*E1(z,x,y);
            vyovyo2(z-1,x-1,y-1)=E2(z,x,y+1)+E2(z,x,y-1)-2*E2(z,x,y);
            vyovyo3(z-1,x-1,y-1)=E3(z,x,y+1)+E3(z,x,y-1)-2*E3(z,x,y);
            vyovyo4(z-1,x-1,y-1)=E4(z,x,y+1)+E4(z,x,y-1)-2*E4(z,x,y);
        end
    end
end
vyovyo1=vyovyo1/(ddiff^2);
vyovyo2=vyovyo2/(ddiff^2);
vyovyo3=vyovyo3/(ddiff^2);
vyovyo4=vyovyo4/(ddiff^2);

vxovyo1=zeros(NZ-2,NX-2,NY-2);
vxovyo2=zeros(NZ-2,NX-2,NY-2);
vxovyo3=zeros(NZ-2,NX-2,NY-2);
vxovyo4=zeros(NZ-2,NX-2,NY-2);
for z=2:NZ-1
    for x=2:NX-1
        for y=2:NY-1
           vxovyo1(z-1,x-1,y-1)=E1(z+1,x,y+1)+E1(z-1,x,y-1)-E1(z-1,x,y+1)-E1(z+1,x,y-1);
           vxovyo2(z-1,x-1,y-1)=E2(z+1,x,y+1)+E2(z-1,x,y-1)-E2(z-1,x,y+1)-E2(z+1,x,y-1);
           vxovyo3(z-1,x-1,y-1)=E3(z+1,x,y+1)+E3(z-1,x,y-1)-E3(z-1,x,y+1)-E3(z+1,x,y-1);
           vxovyo4(z-1,x-1,y-1)=E4(z+1,x,y+1)+E4(z-1,x,y-1)-E4(z-1,x,y+1)-E4(z+1,x,y-1);
        end
    end
end

vxovyo1=vxovyo1/(4*ddiff^2);
vxovyo2=vxovyo2/(4*ddiff^2);
vxovyo3=vxovyo3/(4*ddiff^2);
vxovyo4=vxovyo4/(4*ddiff^2);


vx1=vxo1/h;
vy1=vyo1/h;
vx2=vxo2/h;
vy2=vyo2/h;
vx3=vxo3/h;
vy3=vyo3/h;
vx4=vxo4/h;
vy4=vyo4/h;



Nx=NX-2;
Ny=Nx;
Nz=Nx;


%%

%sigxy calculation
sxy1=zeros;
sxy2=zeros;
sxy3=zeros;
sxy4=zeros;
sxy=zeros;
for z=1:Nz
   for x=1:Nx
    for y=1:Ny
      sxy1(z,x,y)=con*((charge*tau^2*B)*(-df1dE1(z+1,x+1,y+1))*(vx1(z,x,y)^2*vyovyo1(z,x,y)-vx1(z,x,y)*vy1(z,x,y)*vxovyo1(z,x,y)))+ber*berry1(z,x,y)*f1(z+1,x+1,y+1);
      %sxy1(z,x,y)=berry1(z,x,y)*f1(z+1,x+1,y+1);
      %sxy1(z,x,y)=((charge*tau^2*B)*(-df1dE1(z+1,x+1,y+1))*(vx1(z,x,y)^2*vyovyo1(z,x,y)-vx1(z,x,y)*vy1(z,x,y)*vxovyo1(z,x,y)));
                  
      sxy2(z,x,y)=con*((charge*tau^2*B)*(-df2dE2(z+1,x+1,y+1))*(vx2(z,x,y)^2*vyovyo2(z,x,y)-vx2(z,x,y)*vy2(z,x,y)*vxovyo2(z,x,y)))+ber*berry2(z,x,y)*f2(z+1,x+1,y+1);
      %sxy2(z,x,y)=berry2(z,x,y)*f2(z+1,x+1,y+1);
      %sxy2(z,x,y)=((charge*tau^2*B)*(-df2dE2(z+1,x+1,y+1))*(vx2(z,x,y)^2*vyovyo2(z,x,y)-vx2(z,x,y)*vy2(z,x,y)*vxovyo2(z,x,y)));
      
      sxy3(z,x,y)=con*((charge*tau^2*B)*(-df3dE3(z+1,x+1,y+1))*(vx3(z,x,y)^2*vyovyo3(z,x,y)-vx3(z,x,y)*vy3(z,x,y)*vxovyo3(z,x,y)))+ber*berry3(z,x,y)*f3(z+1,x+1,y+1);
      %sxy3(z,x,y)=berry3(z,x,y)*f3(z+1,x+1,y+1);
      %sxy3(z,x,y)=((charge*tau^2*B)*(-df3dE3(z+1,x+1,y+1))*(vx3(z,x,y)^2*vyovyo3(z,x,y)-vx3(z,x,y)*vy3(z,x,y)*vxovyo3(z,x,y)));
      
      sxy4(z,x,y)=con*((charge*tau^2*B)*(-df4dE4(z+1,x+1,y+1))*(vx4(z,x,y)^2*vyovyo4(z,x,y)-vx4(z,x,y)*vy4(z,x,y)*vxovyo4(z,x,y)))+ber*berry4(z,x,y)*f4(z+1,x+1,y+1);
      %sxy4(z,x,y)=berry4(z,x,y)*f4(z+1,x+1,y+1);
      %sxy4(z,x,y)=((charge*tau^2*B)*(-df4dE4(z+1,x+1,y+1))*(vx4(z,x,y)^2*vyovyo4(z,x,y)-vx4(z,x,y)*vy4(z,x,y)*vxovyo4(z,x,y)));
      
      sxy(z,x,y)=sxy1(z,x,y)+sxy2(z,x,y)+sxy3(z,x,y)+sxy4(z,x,y);
    end
  end
end


sigxy=trapzoidl(sxy,Nz,Nx,Ny,deldiff);
sigxy=(charge^2/h)*sigxy;


sigxy;


%%
%alpxy calculation
axy1=zeros;
axy2=zeros;
axy3=zeros;
axy4=zeros;
axy=zeros;
for z=1:Nz
   for x=1:Nx
    for y=1:Ny
      axy1(z,x,y)=con*((charge*tau^2*B/(T*kbolz))*(E1(z+1,x+1,y+1)-mu)*(-df1dE1(z+1,x+1,y+1))*(vx1(z,x,y)^2*vyovyo1(z,x,y)-vx1(z,x,y)*vy1(z,x,y)*vxovyo1(z,x,y)))+ber*berry1(z,x,y)*sf1(z+1,x+1,y+1);
      %axy1(z,x,y)=berry1(z,x,y)*sf1(z+1,x+1,y+1);
      %axy1(z,x,y)=((charge*tau^2*B/(T*kbolz))*(E1(z+1,x+1,y+1)-mu)*(-df1dE1(z+1,x+1,y+1))*(vx1(z,x,y)^2*vyovyo1(z,x,y)-vx1(z,x,y)*vy1(z,x,y)*vxovyo1(z,x,y)));
      
      axy2(z,x,y)=con*((charge*tau^2*B/(T*kbolz))*(E2(z+1,x+1,y+1)-mu)*(-df2dE2(z+1,x+1,y+1))*(vx2(z,x,y)^2*vyovyo2(z,x,y)-vx2(z,x,y)*vy2(z,x,y)*vxovyo2(z,x,y)))+ber*berry2(z,x,y)*sf2(z+1,x+1,y+1);
      %axy2(z,x,y)=berry2(z,x,y)*sf2(z+1,x+1,y+1);
      %axy2(z,x,y)=((charge*tau^2*B/(T*kbolz))*(E2(z+1,x+1,y+1)-mu)*(-df2dE2(z+1,x+1,y+1))*(vx2(z,x,y)^2*vyovyo2(z,x,y)-vx2(z,x,y)*vy2(z,x,y)*vxovyo2(z,x,y)));
      
      axy3(z,x,y)=con*((charge*tau^2*B/(T*kbolz))*(E3(z+1,x+1,y+1)-mu)*(-df3dE3(z+1,x+1,y+1))*(vx3(z,x,y)^2*vyovyo3(z,x,y)-vx3(z,x,y)*vy3(z,x,y)*vxovyo3(z,x,y)))+ber*berry3(z,x,y)*sf3(z+1,x+1,y+1);
      %axy3(z,x,y)=berry3(z,x,y)*sf3(z+1,x+1,y+1);
      %axy3(z,x,y)=((charge*tau^2*B/(T*kbolz))*(E3(z+1,x+1,y+1)-mu)*(-df3dE3(z+1,x+1,y+1))*(vx3(z,x,y)^2*vyovyo3(z,x,y)-vx3(z,x,y)*vy3(z,x,y)*vxovyo3(z,x,y)));
      
      axy4(z,x,y)=con*((charge*tau^2*B/(T*kbolz))*(E4(z+1,x+1,y+1)-mu)*(-df4dE4(z+1,x+1,y+1))*(vx4(z,x,y)^2*vyovyo4(z,x,y)-vx4(z,x,y)*vy4(z,x,y)*vxovyo4(z,x,y)))+ber*berry4(z,x,y)*sf4(z+1,x+1,y+1);
      %axy4(z,x,y)=berry4(z,x,y)*sf4(z+1,x+1,y+1);
      %axy4(z,x,y)=((charge*tau^2*B/(T*kbolz))*(E4(z+1,x+1,y+1)-mu)*(-df4dE4(z+1,x+1,y+1))*(vx4(z,x,y)^2*vyovyo4(z,x,y)-vx4(z,x,y)*vy4(z,x,y)*vxovyo4(z,x,y)));
      
      axy(z,x,y)=axy1(z,x,y)+axy2(z,x,y)+axy3(z,x,y)+axy4(z,x,y);
    end
  end
end


alpxy=trapzoidl(axy,Nz,Nx,Ny,deldiff);
%alpxy
alpxy=(charge*kbolz/h)*alpxy;

alpxy;

%%

%sigxx calculation
sxx1=zeros;
sxx2=zeros;
sxx3=zeros;
sxx4=zeros;
sxx=zeros;
for z=1:Nz
    for x=1:Nx
        for y=1:Ny
         sxx1(z,x,y)=vx1(z,x,y)^2*(-df1dE1(z+1,x+1,y+1));
         sxx2(z,x,y)=vx2(z,x,y)^2*(-df2dE2(z+1,x+1,y+1));
         sxx3(z,x,y)=vx3(z,x,y)^2*(-df3dE3(z+1,x+1,y+1));
         sxx4(z,x,y)=vx4(z,x,y)^2*(-df4dE4(z+1,x+1,y+1));
         sxx(z,x,y)=sxx1(z,x,y)+sxx2(z,x,y)+sxx3(z,x,y)+sxx4(z,x,y);
        end
    end
end


sigxx=trapzoidl(sxx,Nz,Nx,Ny,deldiff);
sigxx=(charge^2)*tau*sigxx;

sigxx;

%%
%alpxx calculation
axx1=zeros;
axx2=zeros;
axx3=zeros;
axx4=zeros;
axx=zeros;

for z=1:Nz
   for x=1:Nx
    for y=1:Ny
      axx1(z,x,y)=vx1(z,x,y)^2*(E1(z+1,x+1,y+1)-mu)*(-df1dE1(z+1,x+1,y+1));
      axx2(z,x,y)=vx2(z,x,y)^2*(E2(z+1,x+1,y+1)-mu)*(-df2dE2(z+1,x+1,y+1));
      axx3(z,x,y)=vx3(z,x,y)^2*(E3(z+1,x+1,y+1)-mu)*(-df3dE3(z+1,x+1,y+1));
      axx4(z,x,y)=vx4(z,x,y)^2*(E4(z+1,x+1,y+1)-mu)*(-df4dE4(z+1,x+1,y+1));
      axx(z,x,y)=axx1(z,x,y)+axx2(z,x,y)+axx3(z,x,y)+axx4(z,x,y);
    end
  end
end


alpxx=trapzoidl(axx,Nz,Nx,Ny,deldiff);
alpxx=-(charge*tau/T)*alpxx;

alpxx;

%%
thetaP=alpxy/alpxx;
thetaH=sigxy/sigxx;

%nernst(bfield)=(alpxx/sigxx)*(thetaP-thetaH);

sigxyp(munum)=sigxy;
alpxyp(munum)=alpxy;


%bfield*100/noofB
munum/noofmu*100

%nernst/T
end

muran1=linspace(0,mumax,noofmu-1);
plot(muran1,diff(sigxyp))
hold on
plot(muran,alpxyp)