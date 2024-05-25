function [S0,Se0,d] = initValBayoma(Fk,fk,f0,zeta0,q)

[mFk,~] = size(Fk) ;
nfk = length(fk) ;

d = 0 ;
A0 = zeros(mFk,mFk) ;
Dk = ( (2*pi*fk).^(-2*q) )./ ( ( 1 - ( f0./fk).^2 ).^2 + ( 2*zeta0*( f0./fk ) ).^2  ) ;
dkDk = 0 ;


for i = 1:1:nfk
    
 d = d + Fk(:,i)'*Fk(:,i) ;
 Dcapk = real(Fk(:,i)*Fk(:,i)') ;
 A0 = A0 + (Dcapk) ;
 
    
end

clear Dcapk

[Phi,Lamda] = eig(A0) ;
lamda = diag(Lamda);
[lamda0,I] = max(lamda) ;
Phi0 = Phi(:,I) ;
[Phi0] = norm1mat(Phi0) ;
clear Phi Lamda lamda I

for i = 1:1:nfk
    
 Dcapk = real(Fk(:,i)*Fk(:,i)') ;
 dcapk = Phi0.'*Dcapk*Phi0 ;
 dkDk =  dkDk + (dcapk)/Dk(i) ;    

end

Se0 = (d - lamda0)/ ( (mFk-1)*nfk ) ;
S0 = (1/nfk)*dkDk ;

clear Fk fk f0 zeta0 q nfk A0 Dk dkDk dcapk lamda0 Phi0

end
