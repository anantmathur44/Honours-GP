function h_hat = em_step(x,Y,n,h_hat,m,k,j)
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

% coo_row = circ(1:n+1,1);

lambda = fft(circ);


seed = k(j);
for i=[1:m]
    Z_sim = circ_simulate_seed(x,circ,mu_hat,sigma_sq_hat,seed+i);
    %Conditional computations 
    %Z_cond = Z_sim(n+2:end)' + c_uo*((c_oo)\(Y- Z_sim(1:n+1)')); 
    %Z_cond = Z_sim(n+2:end)' + c_uo*pcg(c_oo,(Y- Z_sim(1:n+1)'));
    %Z_cond = Z_sim(n+2:end)' + c_uo*(levinson_durbin(coo_row,(Y- Z_sim(1:n+1)')));
    
    sol = conjgrad_circ(lambda/(2*n),Y- Z_sim(1:n+1)', n ,100);
    %sol = levinson_durbin(coo_row,Y- Z_sim(1:n+1)');
    
    %sol = toeplitz(coo_row)\(Y- Z_sim(1:n+1)');
    
    padded_x = [sol ; zeros(n-1,1)];
    %padded_x = [levinson_durbin(coo_row, Y- Z_sim(1:n+1)') ; zeros(n-1,1)];
    %padded_x = [(toeplitz(coo_row)\(Y- Z_sim(1:n+1)')) ; zeros(n-1,1)];
    
    b = matrix_vector_multi(padded_x,lambda);
    Z_cond = Z_sim(n+2:end)' + b(n+2:end);
    Z = [Z Z_cond];
       
end


%   M-step
%Fminsearch options

options = optimoptions('fminunc','Algorithm','quasi-newton','Display','none');


h_hat = fminunc(@(i) -q_step(x,Y,Z,i,m,mu_hat,n) , [h_hat(1),h_hat(2)], options);



