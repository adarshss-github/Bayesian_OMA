function [Dk,Dkf,Dkzeta,Dkinvff,Dkinvfzeta,Dkinvzetazeta,Dkff,Dkfzeta,Dkzetazeta] = Dkinit(x,fk,nfk,q)

%x(1) = fnat x(2) = zeta

betak = x(1)./fk ;

Dk = ( (2*pi*fk).^(-2*q) )./ ( ( 1-(betak).^2 ).^2 + ( 2*x(2)*( betak ) ).^2  ) ;

Dkinvf = 4*(1./fk).*(betak.^3 -betak*(1 - 2*x(2)^2 ) ).*(2*pi*fk).^(2*q) ;
Dkinvzeta = 8*x(2)*( (betak).^2 ).*((2*pi*fk).^(2*q)) ;

Dkf = -(Dk.^2).*(Dkinvf) ;
Dkzeta = -(Dk.^2).*(Dkinvzeta) ;

Dkinvff = 4*(fk.^(-2)).*( 3*betak.^2 - 1 + 2*x(2)^2 ).*(2*pi*fk).^(2*q) ;
Dkinvfzeta = 16*(1./fk)*x(2).*betak.*(2*pi*fk).^(2*q) ; %symmetrical with Dkinvzetaf
Dkinvzetazeta = 8*(betak).^2.*((2*pi*fk).^(2*q)) ;

Dkff =  2*(Dk.^3).*(Dkinvf).*(Dkinvf) - (Dk.^2).*(Dkinvff) ;
Dkfzeta = 2*(Dk.^3).*(Dkinvf).*(Dkinvzeta) - (Dk.^2).*(Dkinvfzeta) ; %symmetrical with Dkzetaf
Dkzetazeta = 2*(Dk.^3).*(Dkinvzeta).*(Dkinvzeta) - (Dk.^2).*(Dkinvzetazeta) ;

clear Dkinvf Dkinvzet


end