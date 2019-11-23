clear all
n= 1000 ; a=0; b=1;
a = 1/(n+1);
x=linspace(a,b,n+1);
mu = 0;
sigma_sq = 1;

global sigma_sq_q;
global mu_q;
sigma_sq_q  =  1;

%% Genereating covariance from kernal 
par = [0.1 1];
%lambda%alpha

cov = generate_cov_row(x,par);

%%Creating circulant row
circ = [cov; cov(end-1:-1:2)];

%SIMULATE KNOWN DATA
Y = circ_simulate(x,circ,mu,sigma_sq);
Y = Y(1:n+1)';


%%%%%%%%%%%%%%%% MLE Alpha
options = optimoptions('fminunc','Algorithm','quasi-newton','Display','iter','MaxFunctionEvaluations', 100000);

tic 
mle_ld = fminunc(@(i)(-likelihood_ld(x,Y,[i, par(2)],mu,n+1)) ,  0.1, options);
toc


options = optimoptions('fminunc','Algorithm','quasi-newton','Display','iter','MaxFunctionEvaluations', 200000,'MaxIterations',10000, 'SpecifyObjectiveGradient',true,'CheckGradients',false);

tic 
%mle_2 = fminunc(@(i)(mm_bound(x,Y,i,n)), [-ones(2*n-2,1);ones(n-1,1); 0.06; 1]' ,options);
toc
%mle_2(3*n-2:end)
tic 
mle_3 = fminunc(@(i)(mm_bound3(x,Y,[i par(2)],n)), [-0.5*ones(2*n-2,1)' 0.1 1] ,options);
%mle_2 = fminunc(@(i)(mm_bound2(x,Y,[i par(2)],n)), [-0.5*ones(2*n,1)' 0.1 1] ,options);
toc
mle_ld
mle_3(2*n)





