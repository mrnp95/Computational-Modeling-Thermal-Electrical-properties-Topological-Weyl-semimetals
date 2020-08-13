
%plot of kz vs E
nz=NZ;
plot(kz(nz/2,:,nz/2),E1(nz/2,:,nz/2))
hold on
plot(kz(nz/2,:,nz/2),E2(nz/2,:,nz/2),'r')
%hold on
%plot(kz(nz/2,:,nz/2),E3(nz/2,:,nz/2),'g')
%hold on
%plot(kz(nz/2,:,nz/2),E4(nz/2,:,nz/2),'y')

%%
%surf of kz,kx vs E  
nz=NZ;
surf(kz(:,:,nz/2),kx(:,:,nz/2),(E1(:,:,nz/2)/t)+.105)
hold on
surf(kz(:,:,nz/2),kx(:,:,nz/2),(E2(:,:,nz/2)/t)-.105)
%hold on
%surf(kz(:,:,nz/2),kx(:,:,nz/2),E3(:,:,nz/2))
%hold on
%surf(kz(:,:,nz/2),kx(:,:,nz/2),E4(:,:,nz/2))
%hold on
%surf(kz(:,:,nz/2),kx(:,:,nz/2),sur(:,:,nz/2))
%hold on
%surf(kz(:,:,nz/2),kx(:,:,nz/2),sur(:,:,nz/2))
%%
%plot of kx vs E
nz=NZ;
plot(kx(:,nz/2+nz/4,nz/2),E1(:,nz/2+nz/4,nz/2))
hold on
plot(kx(:,nz/2+nz/4,nz/2),E2(:,nz/2+nz/4,nz/2),'r')
%hold on
%plot(kx(:,nz/2,nz/2),E3(:,nz/2,nz/2),'g')
%hold on
%plot(kx(:,nz/2,nz/2),E4(:,nz/2,nz/2),'y')
hold on
plot(kx(:,nz/2+nz/4,nz/2),sur(:,nz/2+nz/4,nz/2),'m')
%%
%plot of kx vs vxo
plot(KX(:,nz/2,nz/2),vxo1(:,nz/2,nz/2))
hold on
plot(KX(:,nz/2,nz/2),vxo2(:,nz/2,nz/2),'r')
hold on
plot(KX(:,nz/2,nz/2),vxo3(:,nz/2,nz/2),'g')
hold on
plot(KX(:,nz/2,nz/2),vxo4(:,nz/2,nz/2),'y')
%%
%plot of kz vs E(ky=0) 
plot(kz(nz/2,:),E1(nz/2,:))
hold on
plot(kz(nz/2,:),E2(nz/2,:),'r')
hold on
plot(kz(nz/2,:),E3(nz/2,:),'g')
hold on
plot(kz(nz/2,:),E4(nz/2,:),'y')

%%
%surf of kz,kx vs E(ky=0)
surf(kz(:,:),kx(:,:),E1(:,:))
hold on
surf(kz(:,:),kx(:,:),E2(:,:))
hold on
surf(kz(:,:),kx(:,:),E3(:,:))
hold on
surf(kz(:,:),kx(:,:),E4(:,:))

%%
%plot of kz vs E(ky=kx=0)
plot(kz,E1)
hold on
plot(kz,E2,'r')
hold on
plot(kz,E3,'g')
hold on
plot(kz,E4,'y')

%%
%plot of E vs f
nz=NZ;
plot(E1(nz/2,:,nz/2),f1(nz/2,:,nz/2))
hold on
plot(E2(nz/2,:,nz/2),f2(nz/2,:,nz/2),'r')
hold on
plot(E3(nz/2,:,nz/2),f3(nz/2,:,nz/2),'g')
hold on
plot(E4(nz/2,:,nz/2),f4(nz/2,:,nz/2),'y')
 
%%
%plot of E vs sf
nz=NZ;
plot(E1(nz/2,:,nz/2),sf1(nz/2,:,nz/2))
hold on
plot(E2(nz/2,:,nz/2),sf2(nz/2,:,nz/2),'r')
hold on
plot(E3(nz/2,:,nz/2),sf3(nz/2,:,nz/2),'g')
hold on
plot(E4(nz/2,:,nz/2),sf4(nz/2,:,nz/2),'y')

%%
%plot of kz vs f
nz=NZ;
plot(kz(nz/2,:,nz/2),f1(nz/2,:,nz/2))
hold on
plot(kz(nz/2,:,nz/2),f2(nz/2,:,nz/2),'r')
hold on
plot(kz(nz/2,:,nz/2),f3(nz/2,:,nz/2),'g')
hold on
plot(kz(nz/2,:,nz/2),f4(nz/2,:,nz/2),'y')

%%
%plot of kz vs sf
nz=NZ;
plot(kz(nz/2,:,nz/2),sf1(nz/2,:,nz/2))
hold on
plot(kz(nz/2,:,nz/2),sf2(nz/2,:,nz/2),'r')
hold on
plot(kz(nz/2,:,nz/2),sf3(nz/2,:,nz/2),'g')
hold on
plot(kz(nz/2,:,nz/2),sf4(nz/2,:,nz/2),'y')

%%
%plot of E vs dfdE
nz=NZ;
plot(E1(nz/2,:,nz/2),df1dE1(nz/2,:,nz/2))
hold on
plot(E2(nz/2,:,nz/2),df2dE2(nz/2,:,nz/2),'r')
hold on
plot(E3(nz/2,:,nz/2),df3dE3(nz/2,:,nz/2),'g')
hold on
plot(E4(nz/2,:,nz/2),df4dE4(nz/2,:,nz/2),'y')
 
%%
%plot of kz vs berry
nz=NZ;
plot(KZ(nz/2,:,nz/2),berry1(nz/2,:,nz/2))
hold on
plot(KZ(nz/2,:,nz/2),berry2(nz/2,:,nz/2),'r')
hold on
plot(KZ(nz/2,:,nz/2),berry3(nz/2,:,nz/2),'g')
hold on
plot(KZ(nz/2,:,nz/2),berry4(nz/2,:,nz/2),'y')

%%
%plot of kx vs berry
nz=NZ;
plot(KX(:,nz/2,nz/2),berry1(:,nz/2,nz/2))
hold on
plot(KX(:,nz/2,nz/2),berry2(:,nz/2,nz/2),'r')
hold on
plot(KX(:,nz/2,nz/2),berry3(:,nz/2,nz/2),'g')
hold on
plot(KX(:,nz/2,nz/2),berry4(:,nz/2,nz/2),'y')

 
%%
%surf of kz,kx vs berry
nz=NZ-2;
surf(KZ(:,:,nz/2),KX(:,:,nz/2),berry1(:,:,nz/2))
hold on
surf(KZ(:,:,nz/2),KX(:,:,nz/2),berry2(:,:,nz/2))
hold on
surf(KZ(:,:,nz/2),KX(:,:,nz/2),berry3(:,:,nz/2))
hold on
surf(KZ(:,:,nz/2),KX(:,:,nz/2),berry4(:,:,nz/2))
%%
%surf of kz,kx vs berfac
nz=NZ-2;
surf(KZ(:,:,nz/2),KX(:,:,nz/2),berfac1(:,:,nz/2))
hold on
surf(KZ(:,:,nz/2),KX(:,:,nz/2),berfac2(:,:,nz/2))
%hold on
%surf(KZ(:,:,nz/2),KX(:,:,nz/2),berry3(:,:,nz/2))
%hold on
%surf(KZ(:,:,nz/2),KX(:,:,nz/2),berry4(:,:,nz/2))



%%
%plot of KZ vs sxy
nz=NZ;
plot(KZ(nz/2,:,nz/2),sxy1(nz/2,:,nz/2))
hold on
plot(KZ(nz/2,:,nz/2),sxy2(nz/2,:,nz/2),'r')
%hold on
%plot(KZ(nz/2,:,nz/2),sxy3(nz/2,:,nz/2),'g')
%hold on
%plot(KZ(nz/2,:,nz/2),sxy4(nz/2,:,nz/2),'y')

%%
%surf of KZ,kx vs sxy
nz=NZ;
surf(KZ(:,:,nz/2),KX(:,:,nz/2),sxy1(:,:,nz/2))
hold on
surf(KZ(:,:,nz/2),KX(:,:,nz/2),sxy2(:,:,nz/2))
%hold on
%plot(KZ(nz/2,:,nz/2),sxy3(nz/2,:,nz/2),'g')
%hold on
%plot(KZ(nz/2,:,nz/2),sxy4(nz/2,:,nz/2),'y')
%%
%plot of KZ vs sxy
nz=NZ;
plot(KZ(nz/2,:,nz/2),axy1(nz/2,:,nz/2))
hold on
plot(KZ(nz/2,:,nz/2),axy2(nz/2,:,nz/2),'r')
%hold on
%plot(KZ(nz/2,:,nz/2),sxy3(nz/2,:,nz/2),'g')
%hold on
%plot(KZ(nz/2,:,nz/2),sxy4(nz/2,:,nz/2),'y')

%%
%surf of KZ,kx vs sxy
nz=NZ;
surf(KZ(:,:,nz/2),KX(:,:,nz/2),axy1(:,:,nz/2))
hold on
surf(KZ(:,:,nz/2),KX(:,:,nz/2),axy2(:,:,nz/2))
%hold on
%plot(KZ(nz/2,:,nz/2),sxy3(nz/2,:,nz/2),'g')
%hold on
%plot(KZ(nz/2,:,nz/2),sxy4(nz/2,:,nz/2),'y')

%%
%plot of KZ vs sxy
nz=NZ;
plot(KZ(nz/2,:,nz/2),sxx1(nz/2,:,nz/2))
hold on
plot(KZ(nz/2,:,nz/2),sxx2(nz/2,:,nz/2),'r')
%hold on
%plot(KZ(nz/2,:,nz/2),sxy3(nz/2,:,nz/2),'g')
%hold on
%plot(KZ(nz/2,:,nz/2),sxy4(nz/2,:,nz/2),'y')

%%
%surf of KZ,kx vs sxy
nz=NZ;
surf(KZ(:,:,nz/2),KX(:,:,nz/2),sxx1(:,:,nz/2))
hold on
surf(KZ(:,:,nz/2),KX(:,:,nz/2),sxx2(:,:,nz/2))
%hold on
%plot(KZ(nz/2,:,nz/2),sxy3(nz/2,:,nz/2),'g')
%hold on
%plot(KZ(nz/2,:,nz/2),sxy4(nz/2,:,nz/2),'y')
%%
%plot of KZ vs sxy
nz=NZ;
plot(KZ(nz/2,:,nz/2),axx1(nz/2,:,nz/2))
hold on
plot(KZ(nz/2,:,nz/2),axx2(nz/2,:,nz/2),'r')
%hold on
%plot(KZ(nz/2,:,nz/2),sxy3(nz/2,:,nz/2),'g')
%hold on
%plot(KZ(nz/2,:,nz/2),sxy4(nz/2,:,nz/2),'y')

%%
%surf of KZ,kx vs sxy
nz=NZ;
surf(KZ(:,:,nz/2),KX(:,:,nz/2),axx1(:,:,nz/2))
hold on
surf(KZ(:,:,nz/2),KX(:,:,nz/2),axx2(:,:,nz/2))
%hold on
%plot(KZ(nz/2,:,nz/2),sxy3(nz/2,:,nz/2),'g')
%hold on
%plot(KZ(nz/2,:,nz/2),sxy4(nz/2,:,nz/2),'y')
%%
%surf of kz,kx vs sf
nz=NZ;
surf(kz(:,:,nz/2),kx(:,:,nz/2),sf1(:,:,nz/2))
hold on
surf(kz(:,:,nz/2),kx(:,:,nz/2),sf2(:,:,nz/2))
%hold on
%surf(kz2(:,:,nz/2),kx2(:,:,nz/2),sf21(:,:,nz/2))
%hold on
%surf(kz2(:,:,nz/2),kx2(:,:,nz/2),sf22(:,:,nz/2))

%%
%surf of kz,kx vs axy
nz=NZ;
surf(KZ(:,:,nz/2),KX(:,:,nz/2),axy1(:,:,nz/2))
hold on
surf(KZ(:,:,nz/2),KX(:,:,nz/2),axy2(:,:,nz/2))
%hold on
%surf(kz2(:,:,nz/2),kx2(:,:,nz/2),sf21(:,:,nz/2))
%hold on
%surf(kz2(:,:,nz/2),kx2(:,:,nz/2),sf22(:,:,nz/2))




%%
%plot of kz vs axy
nz=NZ;
plot(KZ(nz/2,:,nz/2),axy1(nz/2,:,nz/2))
hold on
plot(KZ(nz/2,:,nz/2),axy2(nz/2,:,nz/2),'r')
hold on
plot(KZ(nz/2,:,nz/2),axy3(nz/2,:,nz/2),'g')
hold on
plot(KZ(nz/2,:,nz/2),axy4(nz/2,:,nz/2),'y')
%%
%plot of kz vs vxo
nz=NZ;
plot(KZ(nz/2,:,nz/2),vxo1(nz/2,:,nz/2))
hold on
plot(KZ(nz/2,:,nz/2),vxo2(nz/2,:,nz/2),'r')
hold on
plot(KZ(nz/2,:,nz/2),vxo3(nz/2,:,nz/2),'g')
hold on
plot(KZ(nz/2,:,nz/2),vxo4(nz/2,:,nz/2),'y')

%%
%surf of kz,kx vs vxo
nz=NZ;
surf(KZ(:,:,nz/2),KX(:,:,nz/2),vxo1(:,:,nz/2))
hold on
surf(KZ(:,:,nz/2),KX(:,:,nz/2),vxo2(:,:,nz/2))
hold on
surf(KZ(:,:,nz/2),KX(:,:,nz/2),vxo2(:,:,nz/2))
hold on
surf(KZ(:,:,nz/2),KX(:,:,nz/2),vxo2(:,:,nz/2))

%%
%plot of kz vs vyo
plot(KZ(nz/2,:),vyo1(nz/2,:,nz/2))
hold on
plot(KZ(nz/2,:),vyo2(nz/2,:,nz/2),'r')
hold on
plot(KZ(nz/2,:),vyo3(nz/2,:,nz/2),'g')
hold on
plot(KZ(nz/2,:),vyo4(nz/2,:,nz/2),'y')

%%
%plot of kz vs vyovyo
plot(KZ(nz/2,:),vyovyo1(nz/2,:,nz/2))
hold on
plot(KZ(nz/2,:),vyovyo2(nz/2,:,nz/2),'r')
hold on
plot(KZ(nz/2,:),vyovyo3(nz/2,:,nz/2),'g')
hold on
plot(KZ(nz/2,:),vyovyo4(nz/2,:,nz/2),'y')
%%
%plot of kz vs vxovyo
plot(KZ(nz/2,:),vxovyo1(nz/2,:,nz/2))
hold on
plot(KZ(nz/2,:),vxovyo2(nz/2,:,nz/2),'r')
hold on
plot(KZ(nz/2,:),vxovyo3(nz/2,:,nz/2),'g')
hold on
plot(KZ(nz/2,:),vxovyo4(nz/2,:,nz/2),'y')

%%
%plot of kz vs vxo,vyo,vyovyo,vxovyo
plot(KZ(nz/2,:),vxo1(nz/2,:,nz/2))
hold on
plot(KZ(nz/2,:),vyo1(nz/2,:,nz/2),'r')
hold on
plot(KZ(nz/2,:),vxovyo1(nz/2,:,nz/2),'g')
hold on
plot(KZ(nz/2,:),vxovyo1(nz/2,:,nz/2),'y')


%%
%plot of kz vs V1
plot(kz(nz/2,:,nz/2),V11(nz/2,:,nz/2))
hold on
plot(kz(nz/2,:,nz/2),V12(nz/2,:,nz/2),'r')
hold on
plot(kz(nz/2,:,nz/2),V13(nz/2,:,nz/2),'g')

%%
%plot of kz vs V2
plot(kz(nz/2,:,nz/2),V21(nz/2,:,nz/2))
hold on
plot(kz(nz/2,:,nz/2),V22(nz/2,:,nz/2),'r')
hold on
plot(kz(nz/2,:,nz/2),V23(nz/2,:,nz/2),'g')

%%
%plot of kz vs V3
plot(kz(nz/2,:,nz/2),V31(nz/2,:,nz/2))
hold on
plot(kz(nz/2,:,nz/2),V32(nz/2,:,nz/2),'r')
hold on
plot(kz(nz/2,:,nz/2),V33(nz/2,:,nz/2),'g')

%%
%plot of kz vs V4
plot(kz(nz/2,:,nz/2),V41(nz/2,:,nz/2))
hold on
plot(kz(nz/2,:,nz/2),V42(nz/2,:,nz/2),'r')
hold on
plot(kz(nz/2,:,nz/2),V43(nz/2,:,nz/2),'g')

%%
%surf of kz,kx vs V1
surf(kz(:,:,nz/2),kx(:,:,nz/2),real(V11(:,:,nz/2)))
hold on
surf(kz(:,:,nz/2),kx(:,:,nz/2),real(V12(:,:,nz/2)))
hold on
surf(kz(:,:,nz/2),kx(:,:,nz/2),real(V13(:,:,nz/2)))

%%
%surf of kz,kx vs V2
surf(kz(:,:,nz/2),kx(:,:,nz/2),real(V21(:,:,nz/2)))
hold on
surf(kz(:,:,nz/2),kx(:,:,nz/2),real(V22(:,:,nz/2)))
hold on
surf(kz(:,:,nz/2),kx(:,:,nz/2),real(V23(:,:,nz/2)))
%%
%surf of kz,kx vs V3
surf(kz(:,:,nz/2),kx(:,:,nz/2),real(V31(:,:,nz/2)))
hold on
surf(kz(:,:,nz/2),kx(:,:,nz/2),real(V32(:,:,nz/2)))
hold on
surf(kz(:,:,nz/2),kx(:,:,nz/2),real(V33(:,:,nz/2)))
%%
%surf of kz,kx vs V4
surf(kz(:,:,nz/2),kx(:,:,nz/2),real(V41(:,:,nz/2)))
hold on
surf(kz(:,:,nz/2),kx(:,:,nz/2),real(V42(:,:,nz/2)))
hold on
surf(kz(:,:,nz/2),kx(:,:,nz/2),real(V43(:,:,nz/2)))

%%
%plot of kz vs V_1
plot(kz(nz/2,:,nz/2),V11(nz/2,:,nz/2))
hold on
plot(kz(nz/2,:,nz/2),V21(nz/2,:,nz/2),'r')
hold on
plot(kz(nz/2,:,nz/2),V31(nz/2,:,nz/2),'g')
hold on
plot(kz(nz/2,:,nz/2),V41(nz/2,:,nz/2),'y')

%%
%plot of kz vs V_2
plot(kz(nz/2,:,nz/2),V12(nz/2,:,nz/2))
hold on
plot(kz(nz/2,:,nz/2),V22(nz/2,:,nz/2),'r')
hold on
plot(kz(nz/2,:,nz/2),V32(nz/2,:,nz/2),'g')
hold on
plot(kz(nz/2,:,nz/2),V42(nz/2,:,nz/2),'y')

%%
%plot of kz vs V_3
plot(kz(nz/2,:,nz/2),V13(nz/2,:,nz/2))
hold on
plot(kz(nz/2,:,nz/2),V23(nz/2,:,nz/2),'r')
hold on
plot(kz(nz/2,:,nz/2),V33(nz/2,:,nz/2),'g')
hold on
plot(kz(nz/2,:,nz/2),V43(nz/2,:,nz/2),'y')

%%
%surf of kz,kx vs V_1
surf(kz(:,:,nz/2),kx(:,:,nz/2),real(V11(:,:,nz/2)))
hold on
surf(kz(:,:,nz/2),kx(:,:,nz/2),real(V21(:,:,nz/2)))
hold on
surf(kz(:,:,nz/2),kx(:,:,nz/2),real(V31(:,:,nz/2)))
hold on
surf(kz(:,:,nz/2),kx(:,:,nz/2),real(V41(:,:,nz/2)))
%%
%surf of kz,kx vs V_2
surf(kz(:,:,nz/2),kx(:,:,nz/2),real(V12(:,:,nz/2)))
hold on
surf(kz(:,:,nz/2),kx(:,:,nz/2),real(V22(:,:,nz/2)))
hold on
surf(kz(:,:,nz/2),kx(:,:,nz/2),real(V32(:,:,nz/2)))
hold on
surf(kz(:,:,nz/2),kx(:,:,nz/2),real(V42(:,:,nz/2)))

%%
%surf of kz,kx vs V_3
surf(kz(:,:,nz/2),kx(:,:,nz/2),real(V13(:,:,nz/2)))
hold on
surf(kz(:,:,nz/2),kx(:,:,nz/2),real(V23(:,:,nz/2)))
hold on
surf(kz(:,:,nz/2),kx(:,:,nz/2),real(V33(:,:,nz/2)))
hold on
surf(kz(:,:,nz/2),kx(:,:,nz/2),real(V43(:,:,nz/2)))
%%
%plot of KZ vs sxx
plot(KZ(nz/2,:,nz/2),sxx1(nz/2,:,nz/2))
hold on
plot(KZ(nz/2,:,nz/2),sxx2(nz/2,:,nz/2),'r')
hold on
plot(KZ(nz/2,:,nz/2),sxx3(nz/2,:,nz/2),'g')
hold on
plot(KZ(nz/2,:,nz/2),sxx4(nz/2,:,nz/2),'y')
%%
%plot of KZ vs sxy
plot(KZ(nz/2,:,nz/2),sxy1(nz/2,:,nz/2))
hold on
plot(KZ(nz/2,:,nz/2),sxy2(nz/2,:,nz/2),'r')
hold on
plot(KZ(nz/2,:,nz/2),sxy3(nz/2,:,nz/2),'g')
hold on
plot(KZ(nz/2,:,nz/2),sxy4(nz/2,:,nz/2),'y')
%%
%plot of KX vs sxy
plot(KX(:,nz/2,nz/2),sxy1(:,nz/2,nz/2))
hold on
plot(KX(:,nz/2,nz/2),sxy2(:,nz/2,nz/2),'r')
hold on
plot(KX(:,nz/2,nz/2),sxy3(:,nz/2,nz/2),'g')
hold on
plot(KX(:,nz/2,nz/2),sxy4(:,nz/2,nz/2),'y')
%%
%plot of KZ vs axx
nz=NZ-2;
plot(KZ(nz/2,:,nz/2),axx1(nz/2,:,nz/2))
hold on
plot(KZ(nz/2,:,nz/2),axx2(nz/2,:,nz/2),'r')
hold on
plot(KZ(nz/2,:,nz/2),axx3(nz/2,:,nz/2),'g')
hold on
plot(KZ(nz/2,:,nz/2),axx4(nz/2,:,nz/2),'y')
%%
%plot of KX vs axy
nz=NZ-2;
plot(KX(:,nz/2,nz/2),axy1(:,nz/2,nz/2))
hold on
plot(KX(:,nz/2,nz/2),axy2(:,nz/2,nz/2),'r')
hold on
plot(KX(:,nz/2,nz/2),axy3(:,nz/2,nz/2),'g')
hold on
plot(KX(:,nz/2,nz/2),axy4(:,nz/2,nz/2),'y')
%%
%surf of kz,kx vs axy
nz=NZ;
surf(KZ(:,:,nz/2),KX(:,:,nz/2),axx1(:,:,nz/2))
hold on
surf(KZ(:,:,nz/2),KX(:,:,nz/2),axx2(:,:,nz/2))
hold on
surf(KZ(:,:,nz/2),KX(:,:,nz/2),axx3(:,:,nz/2))
hold on
surf(KZ(:,:,nz/2),KX(:,:,nz/2),axx4(:,:,nz/2))
%surf(KZ(:,:),KX(:,:),axx4(:,:,nz/2))
%%
%plot of kx vs axy
nz=NZ;
plot(KX(:,nz/2),axy1(:,nz/2,nz/2))
hold on
plot(KX(:,nz/2),axy2(:,nz/2,nz/2),'r')
hold on
plot(KX(:,nz/2),axy3(:,nz/2,nz/2),'g')
hold on
plot(KX(:,nz/2),axy4(:,nz/2,nz/2),'y')




%%
%plot of k_y vs E_y
plot(k_y,E1_y)
hold on
plot(k_y,E2_y,'r')
hold on
plot(k_y,E3_y,'g')
hold on
plot(k_y,E4_y,'y')

%%
%plot of K_Y vs v_yo
plot(K_Y,v_yo1)
hold on
plot(K_Y,v_yo2,'r')
hold on
plot(K_Y,v_yo3,'g')
hold on
plot(K_Y,v_yo4,'y')
%%
%plot of K_Y vs v_xo
plot(K_Y,v_xo1)
hold on
plot(K_Y,v_xo2,'r')
hold on
plot(K_Y,v_xo3,'g')
hold on
plot(K_Y,v_xo4,'y')
%%
%plot of K_Y vs v_yov_yo
plot(K_Y,v_yov_yo1);
hold on
plot(K_Y,v_yov_yo2,'r');
hold on
plot(K_Y,v_yov_yo3,'g');
hold on
plot(K_Y,v_yov_yo4,'y');
%%
%plot of K_Y vs vxov_yo
plot(K_Y,vxov_yo1);
hold on
plot(K_Y,vxov_yo2,'r');
hold on
plot(K_Y,vxov_yo3,'g');
hold on
plot(K_Y,vxov_yo4,'y');

%%
%surf of K_X,K_Y vs v_xo
surf(K_X,K_Y,v_xo1)
hold on
surf(K_X,K_Y,v_xo2)
hold on
surf(K_X,K_Y,v_xo3)
hold on
surf(K_X,K_Y,v_xo4)
%%
%surf of K_X,K_Y vs v_yo
surf(K_X,K_Y,v_yo1)
hold on
surf(K_X,K_Y,v_yo2)
hold on
surf(K_X,K_Y,v_yo3)
hold on
surf(K_X,K_Y,v_yo4)

%%
%plot of k_y vs E
plot(k_y(:,nz/2),E1_xy(:,nz/2))
hold on
plot(k_y(:,nz/2),E2_xy(:,nz/2),'r')
hold on
plot(k_y(:,nz/2),E3_xy(:,nz/2),'g')
hold on
plot(k_y(:,nz/2),E4_xy(:,nz/2),'y')

%%
%plot of K_Y vs v_xov_yo
plot(K_Y(:,nz/2),v_xov_yo1(:,nz/2));
hold on
plot(K_Y(:,nz/2),v_xov_yo2(:,nz/2),'r');
hold on
plot(K_Y(:,nz/2),v_xov_yo3(:,nz/2),'g');
hold on
plot(K_Y(:,nz/2),v_xov_yo4(:,nz/2),'y');
%%
%plot of K_Y vs v_xo
plot(K_Y(:,nz/2),v_xo1(:,nz/2));
hold on
plot(K_Y(:,nz/2),v_xo2(:,nz/2),'r');
hold on
plot(K_Y(:,nz/2),v_xo3(:,nz/2),'g');
hold on
plot(K_Y(:,nz/2),v_xo4(:,nz/2),'y');
%%
%plot of K_Y vs v_yo
plot(K_Y(:,nz/2),v_yo1(:,nz/2));
hold on
plot(K_Y(:,nz/2),v_yo2(:,nz/2),'r');
hold on
plot(K_Y(:,nz/2),v_yo3(:,nz/2),'g');
hold on
plot(K_Y(:,nz/2),v_yo4(:,nz/2),'y');
%%
%plot of K_Y vs v_yov_yo
plot(K_Y(:,nz/2),v_yov_yo1(:,nz/2));
hold on
plot(K_Y(:,nz/2),v_yov_yo2(:,nz/2),'r');
hold on
plot(K_Y(:,nz/2),v_yov_yo3(:,nz/2),'g');
hold on
plot(K_Y(:,nz/2),v_yov_yo4(:,nz/2),'y');
%%
%plot of K_X vs v_xov_yo
plot(K_X(nz/2,:),v_xo1(nz/2,:));
hold on
plot(K_X(nz/2,:),v_xo2(nz/2,:),'r');
hold on
plot(K_X(nz/2,:),v_xo3(nz/2,:),'g');
hold on
plot(K_X(nz/2,:),v_xo4(nz/2,:),'y');
%%
%plot of K_X vs v_xo
plot(K_X(nz/2,:),v_xov_yo1(nz/2,:));
hold on
plot(K_X(nz/2,:),v_xov_yo2(nz/2,:),'r');
hold on
plot(K_X(nz/2,:),v_xov_yo3(nz/2,:),'g');
hold on
plot(K_X(nz/2,:),v_xov_yo4(nz/2,:),'y');


