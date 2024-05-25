function H = HessL(theta,Phi,Fk,fk,mFk,nfk,q,d,A)

   [Dk,Dkf,Dkzeta,Dkinvff,Dkinvfzeta,Dkinvzetazeta,Dkff,Dkfzeta,Dkzetazeta] = Dkinit([theta(1),theta(2)],fk,nfk,q) ;
   [Lphiphi,qk,Dcapkst] = LQphiphi(theta(4),Fk,fk,mFk,nfk,q,Phi,Dk,A) ;
   [gamk,Bk,Bkf,Bkzeta,BkS,BkSe,Bkff,Bkfzeta,BkfS,BkfSe,Bkzetazeta,BkzetaS,BkzetaSe,BkSS,BkSSe,BkSeSe] = Bkinit([theta(3),theta(4)],nfk,Dk,Dkf,Dkzeta,Dkinvff,Dkinvfzeta,Dkinvzetazeta);
   LDww = funLDww([theta(3),theta(4)],mFk,nfk,Dk,Dkf,Dkzeta,Dkff,Dkzetazeta,Dkfzeta,gamk,Bk) ;
   LQww = funLQww(theta(4),Bk,Bkf,Bkzeta,BkS,BkSe,Bkff,Bkfzeta,BkfS,BkfSe,Bkzetazeta,BkzetaS,BkzetaSe,BkSS,BkSSe,BkSeSe,qk,d) ;
   Lwphi = funLwphi(theta(4),Phi,mFk,nfk,Dcapkst,Bk,Bkf,Bkzeta,BkS,BkSe,A) ;
   Lww = LQww + LDww ;
   H = [Lww Lwphi; Lwphi' Lphiphi] ;
   
    clear Dk Dkf Dkzeta Dkinvff Dkinvfzeta Dkinvzetazeta Dkff Dkfzeta Dkzetazeta
   clear Lphiphi qk Dcapkst LDww LQww Lwphi 
   clear gamk Bk Bkf Bkzeta BkS BkSe Bkff Bkfzeta BkfS BkfSe Bkzetazeta BkzetaS BkzetaSe BkSS BkSSe BkSeSe

end