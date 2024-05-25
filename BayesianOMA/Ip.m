%%
load('MeasuredAccwithSensornoise.mat', 'Xasn')
[Y] = grabData([1,5],Xasn) ;
[modpara,Phi] = BayomaSing(Y,4*1024,100,2,1,0) ;