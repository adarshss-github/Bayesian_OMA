function LQww = funLQww(x,Bk,Bkf,Bkzeta,BkS,BkSe,Bkff,Bkfzeta,BkfS,BkfSe,Bkzetazeta,BkzetaS,BkzetaSe,BkSS,BkSSe,BkSeSe,qk,d)

LQff = -(1/x(1))*sum( Bkff.*qk ) ;
LQfzeta = -(1/x(1))*sum( Bkfzeta.*qk ) ;
LQfS = -(1/x(1))*sum( BkfS.*qk ) ;
LQzetazeta = -(1/x(1))*sum( Bkzetazeta.*qk ) ;
LQzetaS = -(1/x(1))*sum( BkzetaS.*qk ) ;
LQSS = -(1/x(1))*sum( BkSS.*qk ) ;

LQfSe = (1/x(1)^2)*sum( Bkf.*qk) - (1/x(1))*sum(BkfSe.*qk ) ;
LQzetaSe = (1/x(1)^2)*sum( Bkzeta.*qk) - (1/x(1))*sum(BkzetaSe.*qk ) ;
LQSSe = (1/x(1)^2)*sum( BkS.*qk) - (1/x(1))*sum(BkSSe.*qk ) ;
LQSeSe = 2*(1/x(1)^3)*(d - sum(Bk.*qk)) + 2*(1/x(1)^2)*sum(BkSe.*qk) -(1/x(1))*sum( BkSeSe.*qk) ;

LQww = [LQff LQfzeta LQfS LQfSe; LQfzeta LQzetazeta LQzetaS LQzetaSe; ...
    LQfS LQzetaS LQSS LQSSe; LQfSe LQzetaSe LQSSe LQSeSe] ;


















end