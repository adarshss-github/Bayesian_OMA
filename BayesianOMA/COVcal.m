function [fnatCOV,zetaCOV,SCOV,SeCOV,PhiCOV,Covmat] = COVcal(HL,fnat,zeta,S,Se,Phi,Ampv,Lamda)


[mfnat,~] = size(fnat) ;
[mPhi,~] = size(Phi) ;
ntheta = 4 + mPhi ;
Covmat = zeros(4+mPhi,4+mPhi,mfnat) ;
fnatCOV = zeros(mfnat,1) ;
zetaCOV = zeros(mfnat,1) ;
SCOV = zeros(mfnat,1) ;
SeCOV = zeros(mfnat,1) ;
PhiCOV = zeros(mfnat,1) ;
HLdum = zeros(4+mPhi,4+mPhi) ;
Covmatdum = zeros(4+mPhi,4+mPhi) ;
AdumCov = zeros(mPhi,mPhi) ;
dvc  = zeros(4+mPhi,4+mPhi) ;
d2G = zeros(4+mPhi,4+mPhi) ;
d2Lc = zeros(4+mPhi,4+mPhi) ;


for i = 1:1:mfnat
    
    HLdum = HL(:,:,i) ;
    AdumCov = Ampv(:,:,i);
    dvc = [eye(4,4) zeros(4,mPhi); zeros(mPhi,4) eye(mPhi,mPhi)-Phi(:,i)*Phi(:,i)'] ;
   
    d2G = [zeros(4,4) zeros(4,mPhi); zeros(mPhi,4) -2*eye(mPhi,mPhi)] ;
    d2Lc = dvc'*( HLdum + (-1/ Se(i,1) )*Lamda(i)*d2G)*dvc ;
   [~,Sing,~] = svd(d2Lc) ;
    Covmatdum = dvc*pinv(d2Lc,Sing(4+mPhi-1,4+mPhi-1))*dvc' ;
    clear Sing
    Covmat(:,:,i) = Covmatdum ;
    fnatCOV(i,1) = ( sqrt(Covmatdum(1,1) )/fnat(i,1))*100 ;
    zetaCOV(i,1) = ( sqrt(Covmatdum(2,2) )/zeta(i,1))*100 ;
    SCOV(i,1) = ( sqrt(Covmatdum(3,3))/ S(i,1) )*100 ;
    SeCOV(i,1) = ( sqrt(Covmatdum(4,4))/Se(i,1) )*100 ;
    D = eig(Covmatdum(5:end,5:end)) ;
   PhiCOV(i,1) = sqrt(sum(D))*100 ;


end

clear D fnat HL HLdum i mfnat mPhi S Se zeta Covmatdum Ampv AdumCov Lamda dvc d2G d2Lc 


end