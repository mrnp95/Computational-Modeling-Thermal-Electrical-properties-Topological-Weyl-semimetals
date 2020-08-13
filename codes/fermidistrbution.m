clear 
clc
%{
temp=100
et=50*10^-3;     %chemical potential
mu=5*et/5;         %.01 ev
%beta=1/(8.617*10^-5*temp);        %temp =100K
beta=1/(10^-4*temp);

k=linspace(-3.16,3.16,10000);
E=et*k;
%E=linspace(-1,1,10000);
%}

%f=1./(1+exp((E-mu)*beta));
N=100000;
maxi=.0000000001;
mini=1-maxi;

f1=linspace(0,maxi,N);
f2=linspace(mini,1,N);
%delf=(beta*exp(beta*(E - mu)))./((exp(beta*(E - mu)) + 1).^2);
%sf=-(1./(1+exp((E-mu)*beta))).*log(1./(1+exp((E-mu)*beta)))-(exp(E-mu)*beta./(1+exp((E-mu)*beta))).*log(exp(E-mu)*beta./(1+exp((E-mu)*beta)));
sf1=-f1.*log(f1)-(1-f1).*log(1-f1);
sf2=-f2.*log(f2)-(1-f2).*log(1-f2);

%x=log((1-f)./f);
%y=log(1+exp(x))-x.*exp(x);

%plot(E,f)
%hold on
%plot(E,delf,'r')
%hold on
%plot(E,sf,'g')
%plot(f,y)

%plot(f1,sf1)
%hold on
plot(f2,sf2)
