function L = nllfsdof(x,Fk,fk,mFk,nfk,q,d)

N1 = mFk*nfk*log(pi) ;
N2 = 0 ;
N3 = (mFk - 1)*nfk*log(x(4)) ;
Dcapk = zeros(mFk,mFk) ;
A = zeros(mFk,mFk) ;

Dk = ( (2*pi*fk).^(-2*q) )./ ( ( 1-( x(1)./fk).^2 ).^2 + ( 2*x(2)*( x(1)./fk ) ).^2  ) ;

for i = 1:1:nfk
    
    N2 = N2 + log( x(3)*Dk(i) + x(4) ) ;
    Dcapk = real(Fk(:,i)*Fk(:,i)') ;
    A = A + ( 1/( 1 + (x(4))/( x(3)*Dk(i) ) ) )*Dcapk ;
    
end

lamda = max(eig(A)) ;

N4 = (d - lamda)/x(4) ;

L = N1 + N2 + N3 + N4 ;

end