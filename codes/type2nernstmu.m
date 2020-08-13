clear 
%clc


temp=200;                %temp 100K
T=temp;                      %%%onsagar coefficient find it's value
kbolz=8.61733*10^-5;
beta=1/(kbolz*T);
charge=1;
h=6.58*10^(-16);
tau=10^-14;


%mu=0.16*et/5;
mu=0;
noofmu=5;
diffofmu=.002;
mumax=mu+(noofmu)*diffofmu;
mu=mu+diffofmu;


Q=2*10^8;
v=1*10^6;
%C1=0;
%C2=0;





NZ=30;                        %no of values of kz
NX=NZ;
NY=NZ;

%uppr=3.14*10^0;                   %uppr limit of kz
%lowr=(-1)*uppr;                %lowr limit of kz
order=10^8;
upprz1=3.14*order;
lowrz1=.86*order;
upprx1=3.14*order;
lowrx1=-3.14*order;
uppry1=3.14*order;
lowry1=-3.14*order;

upprz2=-lowrz1;
lowrz2=-upprz1;
upprx2=3.14*order;
lowrx2=-3.14*order;
uppry2=3.14*order;
lowry2=-3.14*order;


deldiff=abs((upprz1-lowrz1)/(2*2*pi*NZ)); %difference between each value of kz
ddiff=abs((upprz1-lowrz1)/NZ);

%ky=0;
%kx=0;
kz1=linspace(lowrz1,upprz1,NZ);
kx1=linspace(lowrx1,upprx1,NX);
ky1=linspace(lowry1,uppry1,NY);
%ky=0;
%kx=0;
[kz1,kx1,ky1]=meshgrid(kz1,kx1,ky1);
%[kz,kx]=meshgrid(kz,kx);


kz2=linspace(lowrz2,upprz2,NZ);
kx2=linspace(lowrx2,upprx2,NX);
ky2=linspace(lowry2,uppry2,NY);
%ky=0;
%kx=0;
[kz2,kx2,ky2]=meshgrid(kz2,kx2,ky2);
%[kz,kx]=meshgrid(kz,kx);


KZ1=linspace(lowrz1,upprz1,NZ-2);
KX1=linspace(lowrx1,upprx1,NX-2);
KY1=linspace(lowry1,uppry1,NY-2);

%[Kz,Kx,Ky]=meshgrid(Kz,KX,Ky);
[KZ1,KX1,KY1]=meshgrid(KZ1,KX1,KY1);

KZ2=linspace(lowrz2,upprz2,NZ-2);
KX2=linspace(lowrx2,upprx2,NX-2);
KY2=linspace(lowry2,uppry2,NY-2);

%[Kz,Kx,Ky]=meshgrid(Kz,KX,Ky);
[KZ2,KX2,KY2]=meshgrid(KZ2,KX2,KY2);

s=0;
n=20;
%while mu<=mumax
%sigp=zeros(noofmu,2*n);
cp=zeros(1,2*n);
c=linspace(-4,4,n);


sigt=zeros(1,n);
alpt=zeros(1,n);
sigl=zeros(1,n);
alpl=zeros(1,n);
thetaH=zeros(1,n);
thetaP=zeros(1,n);
nernst=zeros(1,n);


%while s<noofmu
    


%mu=0.02;
MU=linspace(.01,.07,n);
for i=1:n
    
    C1=2*v;
    C2=-C1;
%C1=-2*v;
%C2=-C1;
  mu=MU(i);  
    
H1=zeros(2,2,NZ,NX,NY);
H2=zeros(2,2,NZ,NX,NY);

H1(1,1,:,:,:)=h*Q*v - h*C1*(Q - kz1) - h*kz1*v;
H1(1,2,:,:,:)=- h*kx1*v + h*ky1*v*1i;
H1(2,1,:,:,:)=- h*kx1*v - h*ky1*v*1i;
H1(2,2,:,:,:)=h*kz1*v - h*Q*v - h*C1*(Q - kz1);

H2(1,1,:,:,:)=h*Q*v + h*kz2*v + h*C2*(Q + kz2);
H2(1,2,:,:,:)=h*kx2*v - h*ky2*v*1i;
H2(2,1,:,:,:)=h*kx2*v + h*ky2*v*1i;
H2(2,2,:,:,:)=h*C2*(Q + kz2) - h*kz2*v - h*Q*v;

E11=zeros(NZ,NX,NY);
E12=zeros(NZ,NX,NY);
E21=zeros(NZ,NX,NY);
E22=zeros(NZ,NX,NY);


efun11=zeros(2,NZ,NX,NY);
efun12=zeros(2,NZ,NX,NY);
efun21=zeros(2,NZ,NX,NY);
efun22=zeros(2,NZ,NX,NY);

for z=1:NZ
    for x=1:NX
        for y=1:NY
            [V1,D1]=eig(H1(:,:,z,x,y));
            [V2,D2]=eig(H2(:,:,z,x,y));
            %dd1=eig(H1(:,:,z,x,y));
            %dd2=eig(H2(:,:,z,x,y));
            E11(z,x,y)=D1(1,1);
            E12(z,x,y)=D1(2,2);
            E21(z,x,y)=D2(1,1);
            E22(z,x,y)=D2(2,2);
            
            
            efun11(:,z,x,y)=V1(:,1);
            efun12(:,z,x,y)=V1(:,2);
            efun21(:,z,x,y)=V2(:,1);
            efun22(:,z,x,y)=V2(:,2);
            
        end
    end
end


%%%
%Differentiation of Hamiltonian
delHkx1=zeros(2,2,NZ,NX,NY);
delHkx2=zeros(2,2,NZ,NX,NY);

delHky1=zeros(2,2,NZ,NX,NY);
delHky2=zeros(2,2,NZ,NX,NY);

%delHkx
delHkx1(1,1,:,:,:)=0;
delHkx1(2,1,:,:,:)=-h*v;
delHkx1(1,2,:,:,:)=-h*v;
delHkx1(2,2,:,:,:)=0;

delHkx2(1,1,:,:,:)=0;
delHkx2(2,1,:,:,:)=h*v;
delHkx2(1,2,:,:,:)=h*v;
delHkx2(2,2,:,:,:)=0;

delHky1(1,1,:,:,:)=0;
delHky1(2,1,:,:,:)=-h*v*1i;
delHky1(1,2,:,:,:)=h*v*1i;
delHky1(2,2,:,:,:)=0;

delHky2(1,1,:,:,:)=0;
delHky2(2,1,:,:,:)=h*v*1i;
delHky2(1,2,:,:,:)=-h*v*1i;
delHky2(2,2,:,:,:)=0;

 
%%
%fermi distribution
f11=1./(1+exp((E11-mu)*beta));
dfdE11=-(beta*exp(beta*(E11 - mu)))./((exp(beta*(E11 - mu)) + 1).^2);    %actually delf/delE
sf11=-f11.*log(f11)-(1-f11).*log(1-f11);


f12=1./(1+exp((E12-mu)*beta));
dfdE12=-(beta*exp(beta*(E12 - mu)))./((exp(beta*(E12 - mu)) + 1).^2);    %actually delf/delE
sf12=-f12.*log(f12)-(1-f12).*log(1-f12);


f21=1./(1+exp((E21-mu)*beta));
dfdE21=-(beta*exp(beta*(E21 - mu)))./((exp(beta*(E21 - mu)) + 1).^2);    %actually delf/delE
sf21=-f21.*log(f21)-(1-f21).*log(1-f21);


f22=1./(1+exp((E22-mu)*beta));
dfdE22=-(beta*exp(beta*(E22 - mu)))./((exp(beta*(E22 - mu)) + 1).^2);    %actually delf/delE
sf22=-f22.*log(f22)-(1-f22).*log(1-f22);

%%

% new berry phase
berry1112=zeros(NZ-2,NX-2,NY-2);
%berry13=zeros(NZ-2,NX-2,NY-2);
%berry14=zeros(NZ-2,NX-2,NY-2);

berry1211=zeros(NZ-2,NX-2,NY-2);
%berry23=zeros(NZ-2,NX-2,NY-2);
%berry24=zeros(NZ-2,NX-2,NY-2);

berry2122=zeros(NZ-2,NX-2,NY-2);
%berry32=zeros(NZ-2,NX-2,NY-2);
%berry34=zeros(NZ-2,NX-2,NY-2);

berry2221=zeros(NZ-2,NX-2,NY-2);
%berry42=zeros(NZ-2,NX-2,NY-2);
%berry43=zeros(NZ-2,NX-2,NY-2);

for z=1:NZ-2
    for x=1:NX-2
        for y=1:NY-2
            berry1112(z,x,y)=((efun11(:,z+1,x+1,y+1)'*delHkx1(:,:,z+1,x+1,y+1)*efun12(:,z+1,x+1,y+1))*(efun12(:,z+1,x+1,y+1)'*delHky1(:,:,z+1,x+1,y+1)*efun11(:,z+1,x+1,y+1))-(efun11(:,z+1,x+1,y+1)'*delHky1(:,:,z+1,x+1,y+1)*efun12(:,z+1,x+1,y+1))*(efun12(:,z+1,x+1,y+1)'*delHkx1(:,:,z+1,x+1,y+1)*efun11(:,z+1,x+1,y+1)))/(E11(z+1,x+1,y+1)-E12(z+1,x+1,y+1))^2;
            %berry13(z,x,y)=((efun1(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun3(:,z+1,x+1,y+1))*(efun3(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun1(:,z+1,x+1,y+1))-(efun1(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun3(:,z+1,x+1,y+1))*(efun3(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun1(:,z+1,x+1,y+1)))/(E1(z+1,x+1,y+1)-E3(z+1,x+1,y+1))^2;
            %berry14(z,x,y)=((efun1(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun4(:,z+1,x+1,y+1))*(efun4(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun1(:,z+1,x+1,y+1))-(efun1(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun4(:,z+1,x+1,y+1))*(efun4(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun1(:,z+1,x+1,y+1)))/(E1(z+1,x+1,y+1)-E4(z+1,x+1,y+1))^2;
            
            berry1211(z,x,y)=((efun12(:,z+1,x+1,y+1)'*delHkx1(:,:,z+1,x+1,y+1)*efun11(:,z+1,x+1,y+1))*(efun11(:,z+1,x+1,y+1)'*delHky1(:,:,z+1,x+1,y+1)*efun12(:,z+1,x+1,y+1))-(efun12(:,z+1,x+1,y+1)'*delHky1(:,:,z+1,x+1,y+1)*efun11(:,z+1,x+1,y+1))*(efun11(:,z+1,x+1,y+1)'*delHkx1(:,:,z+1,x+1,y+1)*efun12(:,z+1,x+1,y+1)))/(E12(z+1,x+1,y+1)-E11(z+1,x+1,y+1))^2;
            %berry23(z,x,y)=((efun2(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun3(:,z+1,x+1,y+1))*(efun3(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun2(:,z+1,x+1,y+1))-(efun2(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun3(:,z+1,x+1,y+1))*(efun3(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun2(:,z+1,x+1,y+1)))/(E2(z+1,x+1,y+1)-E3(z+1,x+1,y+1))^2;
            %berry24(z,x,y)=((efun2(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun4(:,z+1,x+1,y+1))*(efun4(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun2(:,z+1,x+1,y+1))-(efun2(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun4(:,z+1,x+1,y+1))*(efun4(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun2(:,z+1,x+1,y+1)))/(E2(z+1,x+1,y+1)-E4(z+1,x+1,y+1))^2;
            
            berry2122(z,x,y)=((efun21(:,z+1,x+1,y+1)'*delHkx2(:,:,z+1,x+1,y+1)*efun22(:,z+1,x+1,y+1))*(efun22(:,z+1,x+1,y+1)'*delHky2(:,:,z+1,x+1,y+1)*efun21(:,z+1,x+1,y+1))-(efun21(:,z+1,x+1,y+1)'*delHky2(:,:,z+1,x+1,y+1)*efun22(:,z+1,x+1,y+1))*(efun22(:,z+1,x+1,y+1)'*delHkx2(:,:,z+1,x+1,y+1)*efun21(:,z+1,x+1,y+1)))/(E21(z+1,x+1,y+1)-E22(z+1,x+1,y+1))^2;
            %berry32(z,x,y)=((efun3(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun2(:,z+1,x+1,y+1))*(efun2(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun3(:,z+1,x+1,y+1))-(efun3(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun2(:,z+1,x+1,y+1))*(efun2(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun3(:,z+1,x+1,y+1)))/(E3(z+1,x+1,y+1)-E2(z+1,x+1,y+1))^2;
            %berry34(z,x,y)=((efun3(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun4(:,z+1,x+1,y+1))*(efun4(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun3(:,z+1,x+1,y+1))-(efun3(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun4(:,z+1,x+1,y+1))*(efun4(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun3(:,z+1,x+1,y+1)))/(E3(z+1,x+1,y+1)-E4(z+1,x+1,y+1))^2;
            
            berry2221(z,x,y)=((efun22(:,z+1,x+1,y+1)'*delHkx2(:,:,z+1,x+1,y+1)*efun21(:,z+1,x+1,y+1))*(efun21(:,z+1,x+1,y+1)'*delHky2(:,:,z+1,x+1,y+1)*efun22(:,z+1,x+1,y+1))-(efun22(:,z+1,x+1,y+1)'*delHky2(:,:,z+1,x+1,y+1)*efun21(:,z+1,x+1,y+1))*(efun21(:,z+1,x+1,y+1)'*delHkx2(:,:,z+1,x+1,y+1)*efun22(:,z+1,x+1,y+1)))/(E22(z+1,x+1,y+1)-E21(z+1,x+1,y+1))^2;
            %berry42(z,x,y)=((efun4(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun2(:,z+1,x+1,y+1))*(efun2(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun4(:,z+1,x+1,y+1))-(efun4(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun2(:,z+1,x+1,y+1))*(efun2(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun4(:,z+1,x+1,y+1)))/(E4(z+1,x+1,y+1)-E2(z+1,x+1,y+1))^2;
            %berry43(z,x,y)=((efun4(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun3(:,z+1,x+1,y+1))*(efun3(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun4(:,z+1,x+1,y+1))-(efun4(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun3(:,z+1,x+1,y+1))*(efun3(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun4(:,z+1,x+1,y+1)))/(E4(z+1,x+1,y+1)-E3(z+1,x+1,y+1))^2;
            
        end
    end
end


berry11=real(1i*(berry1112));
berry12=real(1i*(berry1211));
berry21=real(1i*(berry2122));
berry22=real(1i*(berry2221));

%{
berry1=real(1i*(berry13+berry14));
berry2=real(1i*(berry23+berry24));
berry3=real(1i*(berry31+berry32));
berry4=real(1i*(berry41+berry42));
%}

%%
%velocity vx vy vxvy vyvy
vxo11=zeros(NZ-2,NX-2,NY-2);
vxo12=zeros(NZ-2,NX-2,NY-2);
vxo21=zeros(NZ-2,NX-2,NY-2);
vxo22=zeros(NZ-2,NX-2,NY-2);

for z=2:NZ-1
    for x=2:NX-1
        for y=2:NY-1
            vxo11(z-1,x-1,y-1)=E11(z+1,x,y)-E11(z-1,x,y);
            vxo12(z-1,x-1,y-1)=E12(z+1,x,y)-E12(z-1,x,y);
            vxo21(z-1,x-1,y-1)=E21(z+1,x,y)-E21(z-1,x,y);
            vxo22(z-1,x-1,y-1)=E22(z+1,x,y)-E22(z-1,x,y);
        end
    end
end
vxo11=vxo11/(2*ddiff);
vxo12=vxo12/(2*ddiff);
vxo21=vxo21/(2*ddiff);
vxo22=vxo22/(2*ddiff);


vyo11=zeros(NZ-2,NX-2,NY-2);
vyo12=zeros(NZ-2,NX-2,NY-2);
vyo21=zeros(NZ-2,NX-2,NY-2);
vyo22=zeros(NZ-2,NX-2,NY-2);
for z=2:NZ-1
    for y=2:NY-1
        for x=2:NX-1
            vyo11(z-1,x-1,y-1)=E11(z,x,y+1)-E11(z,x,y-1);
            vyo12(z-1,x-1,y-1)=E12(z,x,y+1)-E12(z,x,y-1);
            vyo21(z-1,x-1,y-1)=E21(z,x,y+1)-E21(z,x,y-1);
            vyo22(z-1,x-1,y-1)=E22(z,x,y+1)-E22(z,x,y-1);
        end
    end
end
vyo11=vyo11/(2*ddiff);
vyo12=vyo12/(2*ddiff);
vyo21=vyo21/(2*ddiff);
vyo22=vyo22/(2*ddiff);

vyovyo11=zeros(NZ-2,NX-2,NY-2);
vyovyo12=zeros(NZ-2,NX-2,NY-2);
vyovyo21=zeros(NZ-2,NX-2,NY-2);
vyovyo22=zeros(NZ-2,NX-2,NY-2);
for z=2:NZ-1
    for y=2:NY-1
        for x=2:NX-1
            vyovyo11(z-1,x-1,y-1)=E11(z,x,y+1)+E11(z,x,y-1)-2*E11(z,x,y);
            vyovyo12(z-1,x-1,y-1)=E12(z,x,y+1)+E12(z,x,y-1)-2*E12(z,x,y);
            vyovyo21(z-1,x-1,y-1)=E21(z,x,y+1)+E21(z,x,y-1)-2*E21(z,x,y);
            vyovyo22(z-1,x-1,y-1)=E22(z,x,y+1)+E22(z,x,y-1)-2*E22(z,x,y);
        end
    end
end
vyovyo11=vyovyo11/(ddiff^2);
vyovyo12=vyovyo12/(ddiff^2);
vyovyo21=vyovyo21/(ddiff^2);
vyovyo22=vyovyo22/(ddiff^2);

vxovyo11=zeros(NZ-2,NX-2,NY-2);
vxovyo12=zeros(NZ-2,NX-2,NY-2);
vxovyo21=zeros(NZ-2,NX-2,NY-2);
vxovyo22=zeros(NZ-2,NX-2,NY-2);
for z=2:NZ-1
    for x=2:NX-1
        for y=2:NY-1
           vxovyo11(z-1,x-1,y-1)=E11(z+1,x,y+1)+E11(z-1,x,y-1)-E11(z-1,x,y+1)-E11(z+1,x,y-1);
           vxovyo12(z-1,x-1,y-1)=E12(z+1,x,y+1)+E12(z-1,x,y-1)-E12(z-1,x,y+1)-E12(z+1,x,y-1);
           vxovyo21(z-1,x-1,y-1)=E21(z+1,x,y+1)+E21(z-1,x,y-1)-E21(z-1,x,y+1)-E21(z+1,x,y-1);
           vxovyo22(z-1,x-1,y-1)=E22(z+1,x,y+1)+E22(z-1,x,y-1)-E22(z-1,x,y+1)-E22(z+1,x,y-1);
        end
    end
end

vxovyo11=vxovyo11/(4*ddiff^2);
vxovyo12=vxovyo12/(4*ddiff^2);
vxovyo21=vxovyo21/(4*ddiff^2);
vxovyo22=vxovyo22/(4*ddiff^2);


vx11=vxo11/h;
vy11=vyo11/h;
vx12=vxo12/h;
vy12=vyo12/h;
vx21=vxo21/h;
vy21=vyo21/h;
vx22=vxo22/h;
vy22=vyo22/h;



Nz=NZ-2;
Nx=Nz;
Ny=Nz;

%%

%sigxy calculation
sxy11=zeros(Nz,Nx,Ny);
sxy12=zeros(Nz,Nx,Ny);
sxy21=zeros(Nz,Nx,Ny);
sxy22=zeros(Nz,Nx,Ny);
%sxy=zeros;

for z=1:Nz
   for x=1:Nx
    for y=1:Ny
      sxy11(z,x,y)=berry11(z,x,y)*f11(z+1,x+1,y+1);
      
                  
      sxy12(z,x,y)=berry12(z,x,y)*f12(z+1,x+1,y+1);
      
      
      sxy21(z,x,y)=berry21(z,x,y)*f21(z+1,x+1,y+1);
     
      
      sxy22(z,x,y)=berry22(z,x,y)*f22(z+1,x+1,y+1);
      
      
      
    end
  end
end

sxy1=sxy11+sxy12;
sxy2=sxy21+sxy22;

sigxy1=trapzoidl(sxy1,Nz,Nx,Ny,deldiff);
sigxy2=trapzoidl(sxy2,Nz,Nx,Ny,deldiff);
sigxy=sigxy1+sigxy2;
sigxy=(charge^2/h)*sigxy;

%sigxy;
sigt(i)=sigxy;

%%
%alpxy calculation
axy11=zeros(Nz,Nx,Ny);
axy12=zeros(Nz,Nx,Ny);
axy21=zeros(Nz,Nx,Ny);
axy22=zeros(Nz,Nx,Ny);

for z=1:Nz
   for x=1:Nx
    for y=1:Ny
      axy11(z,x,y)=berry11(z,x,y)*sf11(z+1,x+1,y+1);
      
      
      axy12(z,x,y)=berry12(z,x,y)*sf12(z+1,x+1,y+1);
      
      
      axy21(z,x,y)=berry21(z,x,y)*sf21(z+1,x+1,y+1);
      
      
      axy22(z,x,y)=berry22(z,x,y)*sf22(z+1,x+1,y+1);
      
      
      
      
    end
  end
end
axy1=axy11+axy12;
axy2=axy21+axy22;

alpxy1=trapzoidl(axy1,Nz,Nx,Ny,deldiff);
alpxy2=trapzoidl(axy2,Nz,Nx,Ny,deldiff);
alpxy=alpxy1+alpxy2;
alpxy=(kbolz*charge/h)*alpxy;


%alpxy;
alpt(i)=alpxy;


%%
%sigxx calculation
sxx11=zeros(Nz,Nx,Ny);
sxx12=zeros(Nz,Nx,Ny);
sxx21=zeros(Nz,Nx,Ny);
sxx22=zeros(Nz,Nx,Ny);

for z=1:Nz
    for x=1:Nx
        for y=1:Ny
         sxx11(z,x,y)=vx11(z,x,y)^2*(-dfdE11(z+1,x+1,y+1));
         sxx12(z,x,y)=vx12(z,x,y)^2*(-dfdE12(z+1,x+1,y+1));
         sxx21(z,x,y)=vx21(z,x,y)^2*(-dfdE21(z+1,x+1,y+1));
         sxx22(z,x,y)=vx22(z,x,y)^2*(-dfdE22(z+1,x+1,y+1));
        end
    end
end

sxx1=sxx11+sxx12;
sxx2=sxx21+sxx22;

sigxx1=trapzoidl(sxx1,Nz,Nx,Ny,deldiff);
sigxx2=trapzoidl(sxx2,Nz,Nx,Ny,deldiff);
sigxx=sigxx1+sigxx2;
sigxx=(charge^2)*tau*sigxx;

%sigxx;
sigl(i)=sigxx;
%%
%alpxx calculation
axx11=zeros(Nz,Nx,Ny);
axx12=zeros(Nz,Nx,Ny);
axx21=zeros(Nz,Nx,Ny);
axx22=zeros(Nz,Nx,Ny);

for z=1:Nz
   for x=1:Nx
    for y=1:Ny
      axx11(z,x,y)=vx11(z,x,y)^2*(E11(z+1,x+1,y+1)-mu)*(-dfdE11(z+1,x+1,y+1));
      axx12(z,x,y)=vx12(z,x,y)^2*(E12(z+1,x+1,y+1)-mu)*(-dfdE12(z+1,x+1,y+1));
      axx21(z,x,y)=vx21(z,x,y)^2*(E21(z+1,x+1,y+1)-mu)*(-dfdE21(z+1,x+1,y+1));
      axx22(z,x,y)=vx22(z,x,y)^2*(E22(z+1,x+1,y+1)-mu)*(-dfdE22(z+1,x+1,y+1));
      
    end
  end
end

axx1=axx11+axx12;
axx2=axx21+axx22;

alpxx1=trapzoidl(axx1,Nz,Nx,Ny,deldiff);
alpxx2=trapzoidl(axx2,Nz,Nx,Ny,deldiff);
alpxx=alpxx1+alpxx2;
alpxx=-(charge*tau/T)*alpxx;

%alpxx;
alpl(i)=alpxx;

%%

thetaH(i)=sigt(i)/sigl(i);
thetaP(i)=alpt(i)/alpl(i);
nernst(i)=(thetaP(i)-thetaH(i))*alpl(i)/sigl(i);

i/n*100
end
%{
thetaH=sigxy/sigxx;
thetaP=alpxy/alpxx;
nernst=(thetaP-thetaH)*alpxy/sigxy

%}

%{
mu=mu+diffofmu;
%plot(c,sig)
for j=1:n
    sigp(s+1,j)=-sig(n+1-j)*10^12;
end
for j=n+1:2*n
    sigp(s+1,j)=sig(j-n)*10^12;
end
%hold on
s=s+1;

%end
for j=1:n
    cp(j)=-c(n+1-j);
end
for j=n+1:2*n
    cp(j)=c(j-n);
end
%}
    
%plot(c,sig)
plot(MU,nernst)
hold on