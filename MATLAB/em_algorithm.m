clear all
n= 20000; a=0; b=1;
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

plot(x,Y)
%MLE

mu_hat = mean(Y);
mu_q = mu_hat;
options = optimoptions('fminunc','Algorithm','quasi-newton','Display','none');

fprintf('\nCholesky time:\n')
% chol(toeplitz(cov));

% tic
% mle_chol = fminunc(@(par)(-likelihood(x,Y,par,mu,n+1)) ,[0.01, 0.01], options);
% toc


fprintf('\nLevinson Durban time:\n')

% 
% tic 
% mle_ld = fminunc(@(par)(-likelihood_ld(x,Y,par,mu,n+1)) , [0.1 0.1], options);
% toc


%EM Algorithm
fprintf('\nMC-EM time:\n')
tic
%Starting search par


start_h = [0.02 0.1];
%%%%%%%%%% par,mu,sigma^2 %%%%%%%%%%%
h_hat = start_h;
h_hat_values_h = [start_h(1)]; 
h_hat_values_a = [start_h(2)]; 
mu_values = [];

%Sample size for monte carlo
m = 10;
%Iterations
iterations = 50;

k = floor(rand(50,1)*100);
% for j=[1:iterations]
%     h_hat = em_step(x,Y,n,h_hat,m,k,j);
%     h_hat_values_h = [h_hat_values_h h_hat(1)]; %x:x-values, Y:Observed data, h_hat:starting par, m = monte carlo size
%     h_hat_values_a = [h_hat_values_a h_hat(2)];
%     h_hat = [h_hat mu_q sigma_sq_q] ;
%     var_mu = mu_q;
%  
%     
% end



start_h = [0.1 0.1];
%%%%%%%%%% par,mu,sigma^2 %%%%%%%%%%%
h_hat = start_h;
h_hat_values_h_1 = [start_h(1)]; 
h_hat_values_a_1 = [start_h(2)]; 

for j=[1:iterations]
    h_hat = em_step_var(x,Y,n,h_hat,m,k,j);
    h_hat_values_h_1 = [h_hat_values_h_1 h_hat(1)]; %x:x-values, Y:Observed data, h_hat:starting par, m = monte carlo size
    h_hat_values_a_1 = [h_hat_values_a_1 h_hat(2)];
    h_hat = [h_hat mu_hat sigma_sq_q] ;
    exact_mu = mu_q;
    
end

toc

% % fprintf('\nCholesky MLE lambda: %d\n',mle_chol(1))
%fprintf('Levinson-Durban MLE lambda msd: %d\n',h_hat_values_h_1(end)-mle_ld(1))
%fprintf('MC-EM MLE lambda mse: %d\n',mle_ld(1)-0.1)
% 
% fprintf('\nCholesky MLE alpha: %d\n',mle_chol(2))
% fprintf('Levinson-Durban MLE alpha: %d\n',mle_ld(2))
 %fprintf('MC-EM MLE alphamsd: %d\n',h_hat_values_a_1(end)-mle_ld(2))
%fprintf('MC-EM MLE alphamse: %d\n',mle_ld(2)-1)
 

% [l,s_chol] = likelihood(x,Y,mle_ld,mu_hat,n+1);
% [l,s_ld] = likelihood_ld(x,Y,mle_chol,mu_hat,n+1);

% fprintf('\nCholesky sigma: %d\n',s_chol)
% fprintf('Levinson-Durban sigma: %d\n',s_ld)
% fprintf('MC-EM sigma: %d\n\n',sigma_sq_q)
% 
% fprintf('MLE mu: %d\n',mu_hat)

%figure(1)
%plot(h_hat_values_h(20:end))



%line([0 iterations-20],[mle_ld(1) mle_ld(1)],'Color',[1 0 0])
%title("Lambda MLE")

%figure(2)
%plot(h_hat_values_a)



line([0 iterations],[mle_ld(2) mle_ld(2)],'Color',[1 0 0])
title("Alpha MLE")


figure(3)
plot(h_hat_values_h_1(20:end))

% 
line([0 iterations-20],[mle_ld(1) mle_ld(1)],'Color',[1 0 0])
title("Lambda(2) MLE")

figure(4)
plot(h_hat_values_a_1)


mean(h_hat_values_a_1(10:end))
line([0 iterations],[mean(h_hat_values_a_1(10:end)) mean(h_hat_values_a_1(10:end))],'Color',[1 0 0])
title("Alpha(2) MLE")

figure(5)
plot(x,Y,'o')

var_h = sum((h_hat_values_h(10:end)-mle_ld(1)).^2)/(iterations-10)
% mean_h = mean(h_hat_values_h(50:end));
var_h_2 = sum((h_hat_values_h_1(10:end)-mle_ld(1)).^2)/(iterations-10)
var_a = sum((h_hat_values_a(10:end)-mle_ld(2)).^2)/(iterations-10)
% mean_a = mean(h_hat_values_a(50:end));
var_a_2 = sum((h_hat_values_a_1(10:end)-mle_ld(2)).^2)/(iterations-10)

