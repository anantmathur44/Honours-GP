clear all
n= 1000; a=0; b=1;
a = 1/(n+1);
x=linspace(a,b,n+1);

mu = 0;
sigma_sq = 1;

global sigma_sq_q;
global mu_q;
sigma_sq_q  =  1;