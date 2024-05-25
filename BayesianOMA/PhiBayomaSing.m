function [Phi,A,lamda0] = PhiBayomaSing(x,Fk,fk,q)

[mFk,~] = size(Fk) ;
 nfk = length(fk) ;
 Dcapk = zeros(mFk,mFk) ;
A = zeros(mFk,mFk) ;
Dk = ( (2*pi*fk).^(-2*q) )./ ( ( 1-( x(1)./fk).^2 ).^2 + ( 2*x(2)*( x(1)./fk ) ).^2  ) ;

for i = 1:1:nfk
    
    Dcapk = real(Fk(:,i)*Fk(:,i)') ;
    A = A + ( 1/( 1 + (x(4))/( x(3)*Dk(i) ) ) )*Dcapk ;
    
end

clear Dcapk Dk

[Phi,Lamda] = eig(A) ;
lamda = diag(Lamda);
[lamda0,I] = max(lamda) ;
Phi= Phi(:,I) ;
[Phi] = norm1mat(Phi) ;


end