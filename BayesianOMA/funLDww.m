function LDww = funLDww(x,mFk,nfk,Dk,Dkf,Dkzeta,Dkff,Dkzetazeta,Dkfzeta,gamk,Bk)

%x(1) = S
%x(2) = Se

LDff = sum( Bk.*(1./Dk).*Dkff - Bk.^2.*((1./Dk).^2).*Dkf.*Dkf ) ;  
LDfzeta = sum( Bk.*(1./Dk).*Dkfzeta - Bk.^2.*((1./Dk).^2).*Dkf.*Dkzeta ) ;  
LDzetazeta = sum( Bk.*(1./Dk).*Dkzetazeta - Bk.^2.*((1./Dk).^2).*Dkzeta.*Dkzeta ) ;

LDfS = (1/x(1))*sum( (1./gamk).*(Bk.^2).*(1./Dk).*Dkf ) ;                    
LDfSe = -(1/x(2))*sum( (1./gamk).*(Bk.^2).*(1./Dk).*Dkf ) ;
LDzetaS = (1/x(1))*sum( (1./gamk).*(Bk.^2).*(1./Dk).*Dkzeta ) ;                    
LDzetaSe = -(1/x(2))*sum( (1./gamk).*(Bk.^2).*(1./Dk).*Dkzeta ) ;

LDSS =  -(1/x(1))^2*sum(Bk.^2) ;
LDSSe = -(1/x(1))*(1/x(2))*sum( (1./gamk).*Bk.^2 ) ;

LDSeSe = -(1/x(2))^2*nfk*(mFk - 1) - (1/x(2))^2*sum(  ((1./gamk).^2).*Bk.^2 ) ;


LDww = [LDff LDfzeta LDfS LDfSe; LDfzeta LDzetazeta LDzetaS LDzetaSe; ...
    LDfS LDzetaS LDSS LDSSe; LDfSe LDzetaSe LDSSe LDSeSe] ;






end