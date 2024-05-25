function [fnat,zeta,S,Se,Phi,snr,fnatCOV,zetaCOV,SCOV,SeCOV,PhiCOV,HL,Covmat] = BayomaSing(Y,N,fs,nsing,opt,q) 

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
%[fnat,zeta,Phi,S,Se] = BayomaSing(Y,N,fs,nsing,opt,q) 
%By Adarsh S, Ph.D. Candidate IIT Kanpur
%
%Description
%------------
%This functions calculates the modal parameters by using Bayesian Operational Modal Analysis (Single mode algorithm)
%
%Input arguments
%----------------
%1) Y is the matrix containing the recored data, arranged in rows
% N,fs,nsing and opt are used for plotting the singular value spectrum
%---------------------------------------------------------------------
%2)If opt = 1, N is the number of points which are used for Fourier Transform, while calculating spectral densities (should be a power of 2),  
%and if opt = 2, N is the number of lags required in the calculation of the cross-correlations (should be even)
%3)fs is the sampling frequency
% nsing is the number of singular values to be plotted (1 to 4)
%4)If opt = 1, the spectral density matrix is calculated by FFT , and if opt = 2, the spectral density matrix is calculated from correlation matrix after applying an exponential window
%5)q = 2 if displacement data is used, q = 1 if velocity data is used and q = 0
%if acceleration data is used 
%
%Output arguments
%----------------
%1)fnat is the MPVs (most probable values) of the identified natural frequencies (Hz)
%2)zeta is the MPVs of the identified damping ratios
%3)S is the MPVs of the identified modal force PSD value (micro*g)^2/Hz
%4)Se is the MPVs of the identified prediction error PSD value (micro*g)^2/Hz (mainly accounts for measurement noise)
%5)Phi is the MPVs of the identified  mode shapes arranged in columns (Normalized to have unit length)
%6)snr is the signal-to-noise ratio of each mode
%7)fnatCOV is the coefficient of variance of natural frequencies in %
%8)zetaCOV is the coefficient of variance of damping ratios in %
%9)SCOV is the coefficient of variance of modal force PSD value in %
%10)SeCOV is the coefficient of variance of sensor noise PSD values in %
%11)PhiCOV is the coefficient of variance of modeshape hyper angles values in %
%12)HL contains the Hessian matrices of the Negative Log Likelihood function
%13)Covmat contains the Covariance matrices 
%
%
%Ex:[fnat,zeta,S,Se,Phi,snr,fnatCOV,zetaCOV,SCOV,SeCOV,PhiCOV,HL,Covmat] = BayomaSing(Y,4*1024,100,2,1,0) ;
%
%Notes:
%------
%1)The algorithm used is the Fast Bayesian FFT algorithm for well separated modes, so the modes need to be well seperated inorder for this
%algorithm to give reliable results
%2) The algorithm uses the scaled FFT valus of the data to formulate the likelihood function
%3) Uniform prior is assumed
%4) Calculation of MPVs is equivalent to Maximum Likelihood estimation
%5) The basic assumptions here are:
%   a) Damping is small and thus the mode shapes will be real
%   b) Scaled FFTs have complex circular symmetric Gaussian distribution and are independent
%   c) Each measurement cahnnel has measurements noises with same PSD 
%   d) Modes are well separated
%6) The covariance matrix is calculated by taking the inverse of the constrained Hessian of the negative log-likelihood function;
%  
%Warning: The Hessian is often ill-conditioned, so Monte Carlo methods
%might have to be use to estimate covariance matrix
%
%Reference: Au, S.K., 2011. Fast Bayesian FFT method for ambient modal identification with separated modes. Journal of Engineering Mechanics, 137(3), pp.214-226.
%
%                                     -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; May, 2020 || *********
%                                     -----------------------------------------------
%
%
%***********************************************************************************************************************************************************
%===========================================================================================================================================================

[Y] = removmean(Y) ; %Makes the data zero mean

[Mk,fini,mY] = filtfftsc(Y,N,fs,nsing,opt) ; % Plots the svd spectrum, and filters out the sclaed fft values of the modal frequency bands based on user input

lfini = length(fini) ;

modpara = zeros(lfini,4) ;
Phi = zeros(mY,lfini) ;
HL = zeros(4+mY,4+mY,lfini) ;
Ampv = zeros(mY,mY,lfini) ;
Lamda = zeros(lfini,1) ;

for i = 1:1:lfini
    
   Fk = Mk{i,1}{1,1} ;
   fk  = Mk{i,1}{2,1} ;
   [mFk,~] = size(Fk) ;
   nfk = length(fk) ;
   f0 = fini(i) ; %Initial value of natural frequency for optimization algorithm
   zeta0 = 0.01 ; %Initial value of damping ratio for optimization algorithm (Change if necessary)
   [S0,Se0,d] = initValBayoma(Fk,fk,f0,zeta0,q) ; %Computes the initial values of S and Se for optimization algorithm 
   f = @(x) nllfsdof(x,Fk,fk,mFk,nfk,q,d) ; %The negative log-likelihood function
   theta = fminsearch(f,[f0,zeta0,S0,Se0]) ; %Optimization step
   modpara(i,:) = theta ;
   [Phi(:,i),Ampv(:,:,i),Lamda(i)] = PhiBayomaSing(theta,Fk,fk,q) ; %Computes the modeshape at MPV

  HL(:,:,i) = HessL(theta,Phi(:,i),Fk,fk,mFk,nfk,q,d,Ampv(:,:,i)) ; %Computes the unconstrained Hessian matrix
 
     clear Fk fk mFk nfk f0 zeta0 S0 Se0 d f theta 

end

[fnatCOV,zetaCOV,SCOV,SeCOV,PhiCOV,Covmat] = COVcal(HL,modpara(:,1), modpara(:,2) ,modpara(:,3),modpara(:,4),Phi,Ampv,Lamda) ; %Computes the C.O.V values of modal parameters and the Covariance matrix

modpara(:,3:4) = modpara(:,3:4)/(10^-6*9.81)^2 ;

fnat = modpara(:,1) ;
zeta = modpara(:,2) ;
S = modpara(:,3) ;
Se = modpara(:,4) ;
snr = S./( (4*zeta.^2).*Se ) ;



clear modpara Mk Ampv Lamada;


end