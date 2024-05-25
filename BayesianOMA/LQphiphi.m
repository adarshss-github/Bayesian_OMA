function [Lphiphi,qk,Dcapkst] = LQphiphi(x,Fk,fk,mFk,nfk,q,Phi,Dk,A)

%x(1) = Se


Dcapk = zeros(mFk,mFk) ;
Dcapkst = zeros(mFk,mFk,nfk) ;
In = eye(mFk) ;
Lphiphi = zeros(mFk,mFk) ;
qk = zeros(1,nfk) ;

for i = 1:1:nfk
    
    Dcapk = real(Fk(:,i)*Fk(:,i)') ;
    Dcapkst(:,:,i) = Dcapk ; 
    qk(i) = Phi'*Dcapk*Phi ;
    
end

Lphiphi = -2*( 1/x(1) )*( A -2*(Phi*Phi')*A + 4*(Phi'*A*Phi)*(Phi*Phi') -2*(A*(Phi*Phi')) - (Phi'*A*Phi)*In ) ;        



end