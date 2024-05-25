function Lwphi = funLwphi(x,Phi,mFk,nfk,Dcapkst,Bk,Bkf,Bkzeta,BkS,BkSe,A)

%x(1) = Se

dumf = zeros(mFk,mFk) ;
dumzeta = zeros(mFk,mFk) ;
dumS = zeros(mFk,mFk) ;
dumSe = zeros(mFk,mFk) ;

for i = 1:1:nfk
    
    dumf = Bkf(i)*reshape(Dcapkst(:,:,i),mFk,mFk) + dumf ;
    dumzeta = Bkzeta(i)*reshape(Dcapkst(:,:,i),mFk,mFk) + dumzeta ;
    dumS = BkS(i)*reshape(Dcapkst(:,:,i),mFk,mFk) + dumS ;
    dumSe = BkSe(i)*reshape(Dcapkst(:,:,i),mFk,mFk) + dumSe ;
    
end

Lfphi = -2*Phi'*(1/x(1))*dumf ;
Lzetaphi = -2*Phi'*(1/x(1))*dumzeta ;
LSphi = -2*Phi'*(1/x(1))*dumS ;

LSephi = -2*Phi'*(1/x(1))*dumSe + 2*(1/(x(1)^2))*Phi'*A ;

Lwphi = [Lfphi;Lzetaphi;LSphi;LSephi] ;

end