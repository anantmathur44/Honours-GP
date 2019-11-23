function score= mc_score_p1(x,Y,n,par,m)
%x:x-values, Y:Observed data, h_hat:starting par, m = monte carlo size

%Generating parameter_t


cov = generate_cov_row(x,par);%Generate covariance row for h_hat
circ = [cov; cov(end-1:-1:2)];  %Generate circulant row (to simulate)
lambda1 = fft(circ);

Z = [];

for i=[1:m]
    Z_sim = circ_simulate_seed(x,circ,0,1,i);   
    Z_sim = conjgrad_circ(lambda1/(2*n), Z_sim(1:n+1)', n ,100);
    Z = [Z Z_sim];       
end


%generating circ row
col = generate_cov_row_d1(x,par);
circ = [col; col(end-1:-1:2)];
prods = [];

lambda = ifft(circ);
lambda = lambda*2*n;



for i = Z
    comp_z = [i; zeros(n-1,1)];
    %comp_z'*(toeplitz(circ))*comp_z
    quad = fft((lambda).*ifft(comp_z));
    quad = comp_z'*quad;
    prods = [prods real(quad)];
    %prods = [prods quadratic_product(lambda,comp_z,1)/n];
end


%trace(toeplitz(col)*inv(toeplitz(cov)))





estimate = mean(prods);

length(Y);
length(lambda1);

ky = conjgrad_circ(lambda1/(2*n), Y, n ,100);
ky = [ky; zeros(n-1,1)];
quad = fft((lambda).*ifft(ky));
quad = ky'*quad;


score = -0.5*estimate + 0.5*real(quad);


