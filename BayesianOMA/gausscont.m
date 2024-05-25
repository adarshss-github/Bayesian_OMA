function [theta] = gausscont(thetastar,covmat,p)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
%[theta] = gausscont(thetastar,covmat,p)
%
%By Adarsh S, Ph.D. Candidate IIT Kanpur
%
%Description
%------------
%This function calculates the values of the contours of a multivariable 2-D
%Gaussian distribution, for a given probability p
%
%Input arguments
%----------------
%1)thetastar is the mean valuse of the varables (column vector)
%2)covmat is the 2x2 covariance matrix
%3) p is the probability for which the contour is to be ploted
%
%Output arguments
%----------------
%1)theta consists of the x and y co-ordinates of the contour, arranged in
%rows, each column is an x-y co-ordinate
%
%Ex: [cont] = gausscont(thetastar,covmat,90/100) ;
%
%                                      -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; May, 2020 || *********
%                                      -----------------------------------------------
%
%
%***********************************************************************************************************************************************************
%===========================================================================================================================================================

L = [ sqrt(covmat(1,1)) 0 ; covmat(2,1)/sqrt(covmat(1,1)) sqrt( covmat(2,2) - ( covmat(2,1)^2/covmat(1,1) ) ) ] ;
beta = 0:0.01:2*pi ;
alpha = sqrt(-2*log(1-p)) ;
y = [alpha*cos(beta);alpha*sin(beta)] ;
theta = thetastar + L*y ;

clear L beta alpha y thetastar covmat
 
end