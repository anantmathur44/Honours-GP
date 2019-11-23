
n= 100; a=0; b=1;
a = 1/(n+1);
x=linspace(a,b,n+1);

mu = 0;
sigma_sq = 1;

global sigma_sq_q;
global mu_q;
sigma_sq_q  =  1;


%% Genereating covariance from kernal 
par = [0.04 2];
%lambda%alpha

cov = generate_cov_row(x,par);

%%Creating circulant row
circ = [cov; cov(end-1:-1:2)];

%SIMULATE KNOWN DATA
Y = circ_simulate(x,circ,mu,sigma_sq);
Y = Y(1:n+1)';

% plot(x,Y)
%MLE


options = optimoptions('fminunc','Algorithm','quasi-newton','Display','none');

fprintf('\nlevinson ime:\n')
tic 
mle_ld = fminunc(@(i)(-likelihood_ld(x,Y,[i par(2)],mu,n+1)) , 0.1, options);
toc



h_hat = [mle_ld par(2)]

mu_hat = 0;
sigma_sq_hat = 1;

m = 2
cov = generate_cov_row(x,h_hat);%Generate covariance row for h_hat
circ = [cov; cov(end-1:-1:2)];  %Generate circulant row (to simulate)
%Generating Z_i(conditional)
Z = [];

%coo_row = circ(1:n+1,1);

lambda = fft(circ);


sol = conjgrad_circ(lambda/(2*n),Y-mu_hat, n ,100);    
padded_x = [sol ; zeros(n-1,1)];
b = matrix_vector_multi(padded_x,lambda);
mode = b(n+2:end)+mu_hat;


for i=[1:m]
    Z_sim = circ_simulate(x,circ,mu_hat,sigma_sq_hat);
    %Conditional computations 
    %Z_cond = Z_sim(n+2:end)' + c_uo*((c_oo)\(Y- Z_sim(1:n+1)')); 
    %Z_cond = Z_sim(n+2:end)' + c_uo*pcg(c_oo,(Y- Z_sim(1:n+1)'));
    %Z_cond = Z_sim(n+2:end)' + c_uo*(levinson_durbin(coo_row,(Y- Z_sim(1:n+1)')));
    sol = conjgrad_circ(lambda/(2*n),Y- Z_sim(1:n+1)', n ,100);
    %sol = levinson_durbin(coo_row,Y- Z_sim(1:n+1)');
    %sol = toeplitz(coo_row)\(Y- Z_sim(1:n+1)');
    padded_x = [sol ; zeros(n-1,1)];
    %padded_x = [levinson_durbin(coo_row,Y- Z_sim(1:n+1)') ; zeros(n-1,1)];
    %padded_x = [(toeplitz(coo_row)\(Y- Z_sim(1:n+1)')) ; zeros(n-1,1)];    
    b = matrix_vector_multi(padded_x,lambda);
    Z_cond = Z_sim(n+2:end)' + b(n+2:end);
    Z = [Z Z_cond];       
end
figure(1)

hax=axes;
hold on;
linear_par = linspace(mle_ld-0.005,mle_ld+0.005,1000);
q3_array = arrayfun(@(i) q_step_var2(x,Y,Z,[i par(2)],m,mu_hat,n,mode,h_hat),linear_par);
p_1 = plot(linear_par, q3_array,'linewidth', 1);

q2_array = arrayfun(@(i) q_step_var(x,Y,Z,[i par(2)],m,mu_hat,n,mode,h_hat),linear_par);
p_2 = plot(linear_par, q2_array,'linewidth', 1);

q1_array = arrayfun(@(i) q_step(x,Y,Z,[i par(2)],m,mu_hat,n),linear_par);
p_3 = plot(linear_par, q1_array,'linewidth', 1);

exact_array = arrayfun(@(i) q_step_exact(x,Y,mode,[i par(2)],m,h_hat),linear_par);

p_4 = plot(linear_par, exact_array,'linewidth', 1);


exact_max = fminunc(@(i) -q_step_exact(x,Y,mode,[i par(2)],m,h_hat) , mle_ld(1), options)

exact_q3 = fminunc(@(i)  -q_step_var2(x,Y,Z,[i par(2)],m,mu_hat,n,mode,h_hat), mle_ld(1),options)


exact_q2 = fminunc(@(i)  -q_step_var(x,Y,Z,[i par(2)],m,mu_hat,n,mode,h_hat), mle_ld(1),options)


exact_q1 = fminunc(@(i)  - q_step(x,Y,Z,[i par(2)],m,mu_hat,n), mle_ld(1),options)

yt = yticks;
yt(1);
 
line([exact_max exact_max], [yt(1) q_step_exact(x,Y,mode,[exact_max par(2)],m,h_hat)],'Color', [0.4940 0.1840 0.5560])
 
line([exact_q3 exact_q3], [yt(1) q_step_var2(x,Y,Z,[exact_q3 par(2)],m,mu_hat,n,mode,h_hat)],'Color', [0 0.4470 0.7410])

line([exact_q2 exact_q2], [yt(1) q_step_var(x,Y,Z,[exact_q2 par(2)],m,mu_hat,n,mode,h_hat)],'Color', [0.8500 0.3250 0.0980])


line([exact_q1 exact_q1], [yt(1)  q_step(x,Y,Z,[exact_q1 par(2)],m,mu_hat,n)],'Color', [0.9290 0.6940 0.1250])

legend('q3','q2','q1', 'Exact')

hold off;





