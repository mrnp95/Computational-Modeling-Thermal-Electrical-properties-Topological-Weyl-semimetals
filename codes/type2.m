clear 
%clc


temp=40;                %temp 100K
T=temp;                      %%%onsagar coefficient find it's value
kbolz=8.61733*10^-5;
beta=1/(kbolz*T);
%mu=0.16*et/5;
mu=0;
noofmu=5;
diffofmu=.002;

mumax=mu+(noofmu)*diffofmu;
mu=mu+diffofmu;
charge=1;
h=6.58*10^(-16);

Q=2;
v=1;
%C1=0;
%C2=0;





NZ=20;                        %no of values of kz
NX=NZ;
NY=NZ;

%uppr=3.14*10^0;                   %uppr limit of kz
%lowr=(-1)*uppr;                %lowr limit of kz
upprz1=3.14;
lowrz1=.86;
upprx1=3.14;
lowrx1=-3.14;
uppry1=3.14;
lowry1=-3.14;

upprz2=-lowrz1;
lowrz2=-upprz1;
upprx2=3.14;
lowrx2=-3.14;
uppry2=3.14;
lowry2=-3.14;


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
n=30;
%while mu<=mumax
sigp=zeros(noofmu,2*n);
cp=zeros(1,2*n);
c=linspace(-1,5,n);
sig=zeros(1,n);

while s<noofmu
    



for i=1:n
    
    C1=c(i);
    C2=-C1;
%C1=-1;
%C2=-C1;
    
    
H1=zeros(2,2,NZ,NX,NY);
H2=zeros(2,2,NZ,NX,NY);

H1(1,1,:,:,:)=Q*v - C1*(Q - kz1) - kz1*v;
H1(1,2,:,:,:)=- kx1*v + ky1*v*1i;
H1(2,1,:,:,:)=- kx1*v - ky1*v*1i;
H1(2,2,:,:,:)=kz1*v - Q*v - C1*(Q - kz1);

H2(1,1,:,:,:)=Q*v + kz2*v + C2*(Q + kz2);
H2(1,2,:,:,:)=kx2*v - ky2*v*1i;
H2(2,1,:,:,:)=kx2*v + ky2*v*1i;
H2(2,2,:,:,:)=C2*(Q + kz2) - kz2*v - Q*v;

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

break
%%%
%Differentiation of Hamiltonian
delHkx1=zeros(2,2,NZ,NX,NY);
delHkx2=zeros(2,2,NZ,NX,NY);

delHky1=zeros(2,2,NZ,NX,NY);
delHky2=zeros(2,2,NZ,NX,NY);

%delHkx
delHkx1(1,1,:,:,:)=0;
delHkx1(2,1,:,:,:)=-v;
delHkx1(1,2,:,:,:)=-v;
delHkx1(2,2,:,:,:)=0;

delHkx2(1,1,:,:,:)=0;
delHkx2(2,1,:,:,:)=v;
delHkx2(1,2,:,:,:)=v;
delHkx2(2,2,:,:,:)=0;

delHky1(1,1,:,:,:)=0;
delHky1(2,1,:,:,:)=-v*1i;
delHky1(1,2,:,:,:)=v*1i;
delHky1(2,2,:,:,:)=0;

delHky2(1,1,:,:,:)=0;
delHky2(2,1,:,:,:)=v*1i;
delHky2(1,2,:,:,:)=-v*1i;
delHky2(2,2,:,:,:)=0;

 
%%
%fermi distribution
f11=1./(1+exp((E11-mu)*beta*10^14*h));
%mf11=exp(E11-mu)*beta./(1+exp((E11-mu)*beta));
%dfdE11=-(beta*exp(beta*(E11 - mu)))./((exp(beta*(E11 - mu)) + 1).^2);    %actually delf/delE
sf11=-f11.*log(f11)-(1-f11).*log(1-f11);


f12=1./(1+exp((E12-mu)*beta*10^14*h));
%mf12=exp(E12-mu)*beta./(1+exp((E12-mu)*beta));
%dfdE12=-(beta*exp(beta*(E12 - mu)))./((exp(beta*(E12 - mu)) + 1).^2);    %actually delf/delE
sf12=-f12.*log(f12)-(1-f12).*log(1-f12);


f21=1./(1+exp((E21-mu)*beta*10^14*h));
%mf21=exp(E21-mu)*beta./(1+exp((E21-mu)*beta));
%dfdE21=-(beta*exp(beta*(E21 - mu)))./((exp(beta*(E21 - mu)) + 1).^2);    %actually delf/delE
sf21=-f21.*log(f21)-(1-f21).*log(1-f21);


f22=1./(1+exp((E22-mu)*beta*10^14*h));
%mf22=exp(E22-mu)*beta./(1+exp((E22-mu)*beta));
%dfdE22=-(beta*exp(beta*(E22 - mu)))./((exp(beta*(E22 - mu)) + 1).^2);    %actually delf/delE
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

Nx=NX-2;
Ny=Nx;
Nz=Nx;

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

sigxy;
constant=charge^2*Q*10^8/(2*3.14*h);
sig(i)=sigxy/constant;

%i/n*100
s*100/noofmu+(i/n*100)/noofmu

end
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

end
for j=1:n
    cp(j)=-c(n+1-j);
end
for j=n+1:2*n
    cp(j)=c(j-n);
end

    
plot(cp,sigp)