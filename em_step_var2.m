function h_hat = em_step_var2(x,Y,n,h_hat,m)
%x:x-values, Y:Observed data, h_hat:starting par, m = monte carlo size
%E-Step
%Generating parameter_t

mu_hat = 0;
sigma_sq_hat = 1;

h_hat = [h_hat(1) h_hat(2)];
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


%M-step
%Fminsearch options

options = optimoptions('fminunc','Algorithm','quasi-newton','Display','none');


h_hat_q3= fminunc(@(i) -q_step_var2(x,Y,Z,i,m,mu_hat,n,mode,h_hat) , [h_hat(1), h_hat(2)], options);
h_hat_exact = fminunc(@(i) -q_step_exact(x,Y,mode,i,m,h_hat) , [h_hat(1), h_hat(2)], options)
%h_hat_original = fminunc(@(i) -q_step(x,Y,Z,i,m,mu_hat,n) , [h_hat(1),h_hat(2)], options);

max_likelihood_q3_par_using_q3 = q_step_var2(x,Y,Z,h_hat_q3,m,mu_hat,n,mode,h_hat);
max_likelihood_exact_par_using_q3 = q_step_var2(x,Y,Z,h_hat_exact,m,mu_hat,n,mode,h_hat);



max_likelihood_q3_par_using_q3 = max_likelihood_q3_par_using_q3
max_likelihood_q3_par_using_exact = q_step_exact(x,Y,mode,h_hat_q3,m,h_hat)

max_likelihood_exact_par_using_q3 = max_likelihood_exact_par_using_q3
max_likelihood_exact_par_using_exact = q_step_exact(x,Y,mode,h_hat_exact,m,h_hat)

h_hat = h_hat_q3;
% vr_max_likelihood = q_step_var2(x,Y,Z,h_hat_vr,m,mu_hat,n,mode,h_hat)
% original_max_likelihood = q_step(x,Y,Z,h_hat_original,m,mu_hat,n)


% em_update_vr = h_hat_vr
% em_update_original = h_hat_original
%em_update_vr = h_hat_new


