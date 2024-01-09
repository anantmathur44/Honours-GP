clear all
n=20000; a=0; b=1;
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


%MLE
options = optimoptions('fminunc','Algorithm','quasi-newton','Display','none');
%MLE Alpha
% tic 
% mle_ld = fminunc(@(i)(-likelihood_ld(x,Y,[par(1), i],mu,n+1)) , 0.1,options);
% toc

% tic 
% mle_tr = fminunc(@(i)(-likelihood_tr(x,Y,[par(1), i],n+1,10)) , 0.1,options);
% toc
% mle_tr

% 
% h_hat = par;
% m = 50;
% %[Z, s2] = mc_score_hutchison(x,Y,n,h_hat,m);
% 
% figure(1)
% linear_par = linspace(mle_ld-0.1,mle_ld+0.1,30);
% score_array = arrayfun(@(i) mc_score(x,Y,n,[par(1), i],m),linear_par);
% p_1 = plot(linear_par, score_array,'linewidth', 1);
% tic
% fzero(@(i) mc_score(x,Y,n,[par(1), i],m),2)
% toc
% mle_ld
% %%%%%%%%%%%%%%%% MLE Alpha
% 
% tic 
% mle_ld = fminunc(@(i)(-likelihood_ld(x,Y,[i, par(2)],mu,n+1)) , 0.1,options);
% toc




h_hat = par;
m = 1;
%[Z, s2] = mc_score_hutchison(x,Y,n,h_hat,m);


%linear_par = linspace(mle_ld-0.05,mle_ld+0.05,30);
%score_array = arrayfun(@(i) mc_score_p1(x,Y,n,[i, par(2)],m),linear_par);
%p_1 = plot(linear_par, score_array,'linewidth', 1);
tic
x_1 = fzero(@(i) mc_score_var(x,Y,n,[i, par(2)],m),0.1)
toc
%mle_ld
%mle_ld


tic
x_2 = fzero(@(i) mc_score_var2(x,Y,n,[par(1), i],1),1.3)
toc

% tic 
% mle_ld2 = fminunc(@(i)(-likelihood_ld(x,Y,[par(1),i],mu,n+1)) , 1.3 ,options)
% toc

mle_ld2-x_2
mle_ld-x_1

