function [gamk,Bk,Bkf,Bkzeta,BkS,BkSe,Bkff,Bkfzeta,BkfS,BkfSe,Bkzetazeta,BkzetaS,BkzetaSe,BkSS,BkSSe,BkSeSe] = Bkinit(x,nfk,Dk,Dkf,Dkzeta,Dkinvff,Dkinvfzeta,Dkinvzetazeta)

%x(1) = S
%x(2) = Se

gamk = ( (x(1))/(x(2)) ).*Dk ;
Bk =  1./( 1 + 1./gamk )  ;


Bkinvf = -(1./gamk).*(1./Dk).*Dkf ;
Bkinvzeta = -(1./gamk).*(1./Dk).*Dkzeta ;
BkinvS = -(1./gamk)*( 1/x(1) ) ;
BkinvSe = (1./gamk)*( 1/x(2)) ;

Bkf = -(Bk.^2).*(Bkinvf) ;
Bkzeta = -(Bk.^2).*(Bkinvzeta) ;
BkS= -(Bk.^2).*(BkinvS) ;
BkSe = -(Bk.^2).*(BkinvSe) ;

Bkinvff = (1./gamk).*Dk.*Dkinvff ;
Bkinvfzeta = (1./gamk).*Dk.*Dkinvfzeta ; %symmetric with Bkinvzetaf
Bkinvzetazeta = (1./gamk).*Dk.*Dkinvzetazeta ;
BkinvfS = ( 1/x(1) )*( 1./gamk ).* (1./Dk).*(Dkf) ;
BkinvfSe = -(1./gamk)*(1/x(2)).*(1./Dk).*(Dkf) ;
BkinvzetaS = ( 1/x(1) )*( 1./gamk ).* (1./Dk).*(Dkzeta) ;
BkinvzetaSe = -(1./gamk)*(1/x(2)).*(1./Dk).*(Dkzeta) ;
BkinvSS = 2*(1./gamk)*(1/x(1))^2 ;
BkinvSSe = -(1./gamk)*(1/x(2))*(1/x(1)) ;
BkinvSeSe = zeros(1,nfk) ;

Bkff = 2*(Bk.^3).*(Bkinvf).*(Bkinvf) - (Bk.^2).*(Bkinvff) ;
Bkfzeta = 2*(Bk.^3).*(Bkinvf).*(Bkinvzeta) - (Bk.^2).*(Bkinvfzeta) ;
Bkzetazeta = 2*(Bk.^3).*(Bkinvzeta).*(Bkinvzeta) - (Bk.^2).*(Bkinvzetazeta) ;
BkfS = 2*(Bk.^3).*(Bkinvf).*(BkinvS) - (Bk.^2).*(BkinvfS) ;
BkfSe = 2*(Bk.^3).*(Bkinvf).*(BkinvSe) - (Bk.^2).*(BkinvfSe) ;
BkzetaS = 2*(Bk.^3).*(Bkinvzeta).*(BkinvS) - (Bk.^2).*(BkinvzetaS) ;
BkzetaSe = 2*(Bk.^3).*(Bkinvzeta).*(BkinvSe) - (Bk.^2).*(BkinvzetaSe) ;
BkSS = 2*(Bk.^3).*(BkinvS).*(BkinvS) - (Bk.^2).*(BkinvSS) ;
BkSSe = 2*(Bk.^3).*(BkinvS).*(BkinvSe) - (Bk.^2).*(BkinvSSe) ;
BkSeSe = 2*(Bk.^3).*(BkinvSe).*(BkinvSe) - (Bk.^2).*(BkinvSeSe) ;







end