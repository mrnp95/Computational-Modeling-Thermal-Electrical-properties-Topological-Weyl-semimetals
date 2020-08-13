clear 
%clc



charge=1;
%tau=7.5*10^-19;
tau=10^-14;
temp=5;                %temp 100K
T=temp;                      %%%onsagar coefficient find it's value
noofT=3;
Tmax=17;
h=(4.135*10^-15)/(2*3.14);
%h=6.626*10^-34/(2*3.14);
kbolz=8.61733*10^-5;
%kbolz=1.38*10^-23;
beta=1/(kbolz*T);
t=.0005*charge;                     %lattice parameter
%bz=10^-5;
%m=9.1*10^-31;
%c=3*10^8;

mu=0.001*charge;                    %chemical potential

consig=1;
conalp=1;
bersig=0;
beralp=0;

NZ=20;                        %no of values of kz
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
sur=zeros(NZ,NX,NY);

KZ=linspace(lowrz,upprz,NZ-2);
KX=linspace(lowrx,upprx,NX-2);
KY=linspace(lowry,uppry,NY-2);

%[Kz,Kx,Ky]=meshgrid(Kz,KX,Ky);
[KZ,KX,KY]=meshgrid(KZ,KX,KY);


%s=0;
%while T<Tmax+1
%Eigenvalues and eigenvector

maxB=200;                             %max value of B
minB=-maxB;                              % min value of B
n=25;                              %no of values of B 
zman=linspace(minB,maxB,n);
%zman=linspace(10,80,n);

%nernstp=zeros(n,1);
sigl=zeros(n,1);
alpl=zeros(n,1);
sigt=zeros(n,1);
alpt=zeros(n,1);
B=zeros(n,1);

%B=5;
%bz=5.8025*10^-5*B;
for i=1:n
    
   B=zman(i);
   %B=3;
   %T=zman(i);
   bz=5.8025*10^-5*B*charge;
   %omega=charge*B/m;
   %omega=1000;
H=zeros(2,2,NZ,NX,NY);


%oroginal hamiltonian
H(1,1,:,:,:)= bz- 2*t*cos(kz);
H(2,1,:,:,:)=2*t*(cos(kz).^2 + cos(ky) - 1) - t*cos(kx) - t*sin(ky)*2i;
H(1,2,:,:,:)=2*t*(cos(kz).^2 + cos(ky) - 1) - t*cos(kx) + t*sin(ky)*2i;
H(2,2,:,:,:)=2*t*cos(kz)-bz;


%{
%H of watzman paper
H(1,1,:,:,:)=-2*t*sin(kz);
H(2,1,:,:,:)= 2*t*(cos(ky).^2 + cos(kz) - 1) - t*cos(kx) - t*cos(ky)*2i;
H(1,2,:,:,:)= 2*t*(cos(ky).^2 + cos(kz) - 1) - t*cos(kx) + t*cos(ky)*2i;
H(2,2,:,:,:)= 2*t*sin(kz);
%}


E1=zeros(NZ,NX,NY);
E2=zeros(NZ,NX,NY);

efun1=zeros(2,NZ,NX,NY);
efun2=zeros(2,NZ,NX,NY);

for z=1:NZ
    for x=1:NX
        for y=1:NY
            [V,D]=eig(H(:,:,z,x,y));
            %dd=eig(H(:,:,z,x,y));
            E1(z,x,y)=D(1,1);
            E2(z,x,y)=D(2,2);
            
            
            efun1(:,z,x,y)=V(:,1);
            efun2(:,z,x,y)=V(:,2);
           
        end
    end
end
%kx=kx*2;
%break
%%
%Differentiation of Hamiltonian
delHkx=zeros(2,2,NZ,NX,NY);
delHky=zeros(2,2,NZ,NX,NY);


%
%original hamiltonian
%delHkx
delHkx(1,1,:,:,:)=0;
delHkx(2,1,:,:,:)=t*sin(kx);
delHkx(1,2,:,:,:)=t*sin(kx);
delHkx(2,2,:,:,:)=0;
            
%delHky
delHky(1,1,:,:,:)=0;
delHky(2,1,:,:,:)=- t*cos(ky)*2i - 2*t*sin(ky);
delHky(1,2,:,:,:)=t*cos(ky)*2i - 2*t*sin(ky);
delHky(2,2,:,:,:)=0;
%

%{
%H of watzman paper
delHkx(1,1,:,:,:)=0;
delHkx(2,1,:,:,:)=t*sin(kx);
delHkx(1,2,:,:,:)=t*sin(kx);
delHkx(2,2,:,:,:)=0;

delHky(1,1,:,:,:)=0;
delHky(2,1,:,:,:)=t*sin(ky)*2i - 4*t*cos(ky).*sin(ky);
delHky(1,2,:,:,:)=- t*sin(ky)*2i - 4*t*cos(ky).*sin(ky);
delHky(2,2,:,:,:)=0;
%}

%fermi distribution
f1=1./(1+exp((E1-mu)*beta));
df1dE1=-(beta*exp(beta*(E1 - mu)))./((exp(beta*(E1 - mu)) + 1).^2);    %actually delf/delE
sf1=-f1.*log(f1)-(1-f1).*log(1-f1);


f2=1./(1+exp((E2-mu)*beta));
df2dE2=-(beta*exp(beta*(E2 - mu)))./((exp(beta*(E2 - mu)) + 1).^2);    %actually delf/delE
sf2=-f2.*log(f2)-(1-f2).*log(1-f2);





% new berry phase
berry12=zeros(NZ-2,NX-2,NY-2);
%berry13=zeros(NZ-2,NX-2,NY-2);
%berry14=zeros(NZ-2,NX-2,NY-2);

berry21=zeros(NZ-2,NX-2,NY-2);
%berry23=zeros(NZ-2,NX-2,NY-2);
%berry24=zeros(NZ-2,NX-2,NY-2);


for z=1:NZ-2
    for x=1:NX-2
        for y=1:NY-2
            berry12(z,x,y)=((efun1(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun2(:,z+1,x+1,y+1))*(efun2(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun1(:,z+1,x+1,y+1))-(efun1(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun2(:,z+1,x+1,y+1))*(efun2(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun1(:,z+1,x+1,y+1)))/(E1(z+1,x+1,y+1)-E2(z+1,x+1,y+1))^2;
            %berry13(z,x,y)=((efun1(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun3(:,z+1,x+1,y+1))*(efun3(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun1(:,z+1,x+1,y+1))-(efun1(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun3(:,z+1,x+1,y+1))*(efun3(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun1(:,z+1,x+1,y+1)))/(E1(z+1,x+1,y+1)-E3(z+1,x+1,y+1))^2;
            %berry14(z,x,y)=((efun1(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun4(:,z+1,x+1,y+1))*(efun4(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun1(:,z+1,x+1,y+1))-(efun1(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun4(:,z+1,x+1,y+1))*(efun4(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun1(:,z+1,x+1,y+1)))/(E1(z+1,x+1,y+1)-E4(z+1,x+1,y+1))^2;
            
            berry21(z,x,y)=((efun2(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun1(:,z+1,x+1,y+1))*(efun1(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun2(:,z+1,x+1,y+1))-(efun2(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun1(:,z+1,x+1,y+1))*(efun1(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun2(:,z+1,x+1,y+1)))/(E2(z+1,x+1,y+1)-E1(z+1,x+1,y+1))^2;
            %berry23(z,x,y)=((efun2(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun3(:,z+1,x+1,y+1))*(efun3(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun2(:,z+1,x+1,y+1))-(efun2(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun3(:,z+1,x+1,y+1))*(efun3(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun2(:,z+1,x+1,y+1)))/(E2(z+1,x+1,y+1)-E3(z+1,x+1,y+1))^2;
            %berry24(z,x,y)=((efun2(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun4(:,z+1,x+1,y+1))*(efun4(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun2(:,z+1,x+1,y+1))-(efun2(:,z+1,x+1,y+1)'*delHky(:,:,z+1,x+1,y+1)*efun4(:,z+1,x+1,y+1))*(efun4(:,z+1,x+1,y+1)'*delHkx(:,:,z+1,x+1,y+1)*efun2(:,z+1,x+1,y+1)))/(E2(z+1,x+1,y+1)-E4(z+1,x+1,y+1))^2;
            
            
        end
    end
end


berry1=real(1i*berry12);
berry2=real(1i*berry21);
%{
berfac1=zeros(NZ-2,NX-2,NY-2);
berfac2=zeros(NZ-2,NX-2,NY-2);

for z=1:NZ-2
    for x=1:NX-2
        for y=1:NY-2
          berfac1(z,x,y)=(1+(charge/(h*c))*B*berry1(z,x,y));
          berfac2(z,x,y)=(1+(charge/(h*c))*B*berry2(z,x,y));
               
        end
    end
end
%}

%break
%%

%velocity vx vy vxvy vyvy
vxo1=zeros(NZ-2,NX-2,NY-2);
vxo2=zeros(NZ-2,NX-2,NY-2);

for z=2:NZ-1
    for x=2:NX-1
        for y=2:NY-1
            vxo1(z-1,x-1,y-1)=E1(z+1,x,y)-E1(z-1,x,y);
            vxo2(z-1,x-1,y-1)=E2(z+1,x,y)-E2(z-1,x,y);
        end
    end
end
vxo1=vxo1/(2*ddiff);
vxo2=vxo2/(2*ddiff);



vyo1=zeros(NZ-2,NX-2,NY-2);
vyo2=zeros(NZ-2,NX-2,NY-2);

for z=2:NZ-1
    for y=2:NY-1
        for x=2:NX-1
            vyo1(z-1,x-1,y-1)=E1(z,x,y+1)-E1(z,x,y-1);
            vyo2(z-1,x-1,y-1)=E2(z,x,y+1)-E2(z,x,y-1);
            
        end
    end
end
vyo1=vyo1/(2*ddiff);
vyo2=vyo2/(2*ddiff);


vyovyo1=zeros(NZ-2,NX-2,NY-2);
vyovyo2=zeros(NZ-2,NX-2,NY-2);

for z=2:NZ-1
    for y=2:NY-1
        for x=2:NX-1
            vyovyo1(z-1,x-1,y-1)=E1(z,x,y+1)+E1(z,x,y-1)-2*E1(z,x,y);
            vyovyo2(z-1,x-1,y-1)=E2(z,x,y+1)+E2(z,x,y-1)-2*E2(z,x,y);
            
        end
    end
end
vyovyo1=vyovyo1/(ddiff^2);
vyovyo2=vyovyo2/(ddiff^2);


vxovyo1=zeros(NZ-2,NX-2,NY-2);
vxovyo2=zeros(NZ-2,NX-2,NY-2);

for z=2:NZ-1
    for x=2:NX-1
        for y=2:NY-1
           vxovyo1(z-1,x-1,y-1)=E1(z+1,x,y+1)+E1(z-1,x,y-1)-E1(z-1,x,y+1)-E1(z+1,x,y-1);
           vxovyo2(z-1,x-1,y-1)=E2(z+1,x,y+1)+E2(z-1,x,y-1)-E2(z-1,x,y+1)-E2(z+1,x,y-1);
           
        end
    end
end

vxovyo1=vxovyo1/(4*ddiff^2);
vxovyo2=vxovyo2/(4*ddiff^2);



vx1=vxo1/h;
vy1=vyo1/h;
vx2=vxo2/h;
vy2=vyo2/h;




Nx=NX-2;
Ny=Nx;
Nz=Nx;


%%

%sigxy calculation
sxy1=zeros(Nz,Nx,Ny);
sxy2=zeros(Nz,Nx,Ny);
%sxy=zeros;
for z=1:Nz
   for x=1:Nx
    for y=1:Ny
      
      sxy1(z,x,y)=consig*((charge*tau^2*B)*(-df1dE1(z+1,x+1,y+1))*(vx1(z,x,y)^2*vyovyo1(z,x,y)-vx1(z,x,y)*vy1(z,x,y)*vxovyo1(z,x,y)))+bersig*berry1(z,x,y)*f1(z+1,x+1,y+1);
      %sxy1(z,x,y)=berry1(z,x,y)*f1(z+1,x+1,y+1);
      %sxy1(z,x,y)=((charge*tau^2*B)*(-df1dE1(z+1,x+1,y+1))*(vx1(z,x,y)^2*vyovyo1(z,x,y)-vx1(z,x,y)*vy1(z,x,y)*vxovyo1(z,x,y)));
                  
      sxy2(z,x,y)=consig*((charge*tau^2*B)*(-df2dE2(z+1,x+1,y+1))*(vx2(z,x,y)^2*vyovyo2(z,x,y)-vx2(z,x,y)*vy2(z,x,y)*vxovyo2(z,x,y)))+bersig*berry2(z,x,y)*f2(z+1,x+1,y+1);
      %sxy2(z,x,y)=berry2(z,x,y)*f2(z+1,x+1,y+1);
      %sxy2(z,x,y)=((charge*tau^2*B)*(-df2dE2(z+1,x+1,y+1))*(vx2(z,x,y)^2*vyovyo2(z,x,y)-vx2(z,x,y)*vy2(z,x,y)*vxovyo2(z,x,y)));
      
      %{
      %sxy1(z,x,y)=(berfac1(z,x,y))^-1*(-df1dE1(z+1,x+1,y+1))*(vy1(z,x,y))^2*(berfac1(z,x,y)^2+omega^2*tau^2)^-1;
      %sxy2(z,x,y)=(berfac2(z,x,y))^-1*(-df2dE2(z+1,x+1,y+1))*(vy2(z,x,y))^2*(berfac2(z,x,y)^2+omega^2*tau^2)^-1;
      %}
      
    end
  end
end

sxy=sxy1+sxy2;
sigxy=trapzoidl(sxy,Nz,Nx,Ny,deldiff);
%sigxy=(charge^2*tau^2*omega)*sigxy;
sigxy=(charge^2/h)*sigxy;

sigt(i)=sigxy;


%%
%alpxy calculation
axy1=zeros(Nz,Nx,Ny);
axy2=zeros(Nz,Nx,Ny);
%axy=zeros;
for z=1:Nz
   for x=1:Nx
    for y=1:Ny
       
      axy1(z,x,y)=conalp*((charge*tau^2*B/(T*kbolz))*(E1(z+1,x+1,y+1)-mu)*(-df1dE1(z+1,x+1,y+1))*(vx1(z,x,y)^2*vyovyo1(z,x,y)-vx1(z,x,y)*vy1(z,x,y)*vxovyo1(z,x,y)))+beralp*berry1(z,x,y)*sf1(z+1,x+1,y+1);
      %axy1(z,x,y)=berry1(z,x,y)*sf1(z+1,x+1,y+1);
      %axy1(z,x,y)=((charge*tau^2*B/(T*kbolz))*(E1(z+1,x+1,y+1)-mu)*(-df1dE1(z+1,x+1,y+1))*(vx1(z,x,y)^2*vyovyo1(z,x,y)-vx1(z,x,y)*vy1(z,x,y)*vxovyo1(z,x,y)));
      
      axy2(z,x,y)=conalp*((charge*tau^2*B/(T*kbolz))*(E2(z+1,x+1,y+1)-mu)*(-df2dE2(z+1,x+1,y+1))*(vx2(z,x,y)^2*vyovyo2(z,x,y)-vx2(z,x,y)*vy2(z,x,y)*vxovyo2(z,x,y)))+beralp*berry2(z,x,y)*sf2(z+1,x+1,y+1);
      %axy2(z,x,y)=berry2(z,x,y)*sf2(z+1,x+1,y+1);
      %axy2(z,x,y)=((charge*tau^2*B/(T*kbolz))*(E2(z+1,x+1,y+1)-mu)*(-df2dE2(z+1,x+1,y+1))*(vx2(z,x,y)^2*vyovyo2(z,x,y)-vx2(z,x,y)*vy2(z,x,y)*vxovyo2(z,x,y)));
      
      %{  
      %axy1(z,x,y)=(berfac1(z,x,y))^-1*(E1(z+1,x+1,y+1)-mu)*(-df1dE1(z+1,x+1,y+1))*(vy1(z,x,y))^2*(berfac1(z,x,y)^2+omega^2*tau^2)^-1; 
      %axy1(z,x,y)=(berfac1(z,x,y))^-1*(E1(z+1,x+1,y+1)-mu)*(-df1dE1(z+1,x+1,y+1))*(vy1(z,x,y))^2;
      %axy2(z,x,y)=(berfac2(z,x,y))^-1*(E2(z+1,x+1,y+1)-mu)*(-df2dE2(z+1,x+1,y+1))*(vy2(z,x,y))^2*(berfac2(z,x,y)^2+omega^2*tau^2)^-1; 
      %axy2(z,x,y)=(berfac2(z,x,y))^-1*(E2(z+1,x+1,y+1)-mu)*(-df2dE2(z+1,x+1,y+1))*(vy2(z,x,y))^2;
      %}
    end
  end
end
axy=axy1+axy2;

alpxy=trapzoidl(axy,Nz,Nx,Ny,deldiff);
%alpxy=-(charge*tau^2*omega/T)*alpxy;
alpxy=(charge*kbolz/h)*alpxy;

alpt(i)=alpxy;
%break
%%

%sigxx calculation
sxx1=zeros(Nz,Nx,Ny);
sxx2=zeros(Nz,Nx,Ny);
%sxx=zeros;
for z=1:Nz
    for x=1:Nx
        for y=1:Ny
         
         sxx1(z,x,y)=vx1(z,x,y)^2*(-df1dE1(z+1,x+1,y+1));
         sxx2(z,x,y)=vx2(z,x,y)^2*(-df2dE2(z+1,x+1,y+1));
         %{
         sxx1(z,x,y)=(berfac1(z,x,y))*(-df1dE1(z+1,x+1,y+1))*(vx1(z,x,y))^2*((1+berfac1(z,x,y))^2+omega^2*tau^2)^-1;
         sxx2(z,x,y)=(berfac2(z,x,y))*(-df2dE2(z+1,x+1,y+1))*(vx2(z,x,y))^2*((1+berfac2(z,x,y))^2+omega^2*tau^2)^-1; 
         %}
        end
    end
end

sxx=sxx1+sxx2;
sigxx=trapzoidl(sxx,Nz,Nx,Ny,deldiff);
sigxx=(charge^2)*tau*sigxx;

sigl(i)=sigxx;

%%
%alpxx calculation
axx1=zeros(Nz,Nx,Ny);
axx2=zeros(Nz,Nx,Ny);
%axx=zeros;

for z=1:Nz
   for x=1:Nx
    for y=1:Ny
      
      axx1(z,x,y)=vx1(z,x,y)^2*(E1(z+1,x+1,y+1)-mu)*(-df1dE1(z+1,x+1,y+1));
      axx2(z,x,y)=vx2(z,x,y)^2*(E2(z+1,x+1,y+1)-mu)*(-df2dE2(z+1,x+1,y+1));
      %{
      axx1(z,x,y)=(berfac1(z,x,y))*(E1(z+1,x+1,y+1)-mu)*(-df1dE1(z+1,x+1,y+1))*(vx1(z,x,y))^2*((1+berfac1(z,x,y))^2+omega^2*tau^2)^-1;   
      axx2(z,x,y)=(berfac2(z,x,y))*(E2(z+1,x+1,y+1)-mu)*(-df2dE2(z+1,x+1,y+1))*(vx2(z,x,y))^2*((1+berfac2(z,x,y))^2+omega^2*tau^2)^-1;   
      %}
      
    end
  end
end

axx=axx1+axx2;
alpxx=trapzoidl(axx,Nz,Nx,Ny,deldiff);
alpxx=-(charge*tau/T)*alpxx;

alpl(i)=alpxx;
%break
%%
%{
thetaP=alpxy/alpxx;
thetaH=sigxy/sigxx;


nernst=(alpxx/sigxx)*(thetaP-thetaH);

%Tnernst=nernst/T
nernstp(i)=nernst;
%s*100/noofT+(magfld/noofB*100)/noofT
%}
i/n*100


end

thetaP=alpt./alpl;
thetaH=sigt./sigl;
nernstp=(alpl./sigl).*(thetaP-thetaH);

%{
Bp=zeros(2*noofB,1);
Tnernstp=zeros(2*noofB,1);
for r=1:noofB
    Bp(r)=-zman(noofB+1-r)/10;
    Tnernstp(r)=-Tnernst(noofB+1-r)*10^-12;
end
for r=noofB+1:2*noofB
    Bp(r)=zman(r-noofB)/10;
    Tnernstp(r)=Tnernst(r-noofB)*10^-12;
end 


%{
Bp=zeros(noofB,1);
Tnernstp=zeros(noofB,1);
for r=1:noofB
    Bp(r)=zman(r);
    Tnernstp(r)=Tnernst(r);
end
%}
plot(Bp,Tnernstp)
hold on
s=s+1;


T=T+3;
%end
%}
zmann=zman/20;
nernstn=nernstp/T;
plot(zmann,nernstp,'p')
grid on
hold on
