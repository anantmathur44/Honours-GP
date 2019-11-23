function [s1,s2]= mc_score_hutchison(x,Y,n,par,m)
%x:x-values, Y:Observed data, h_hat:starting par, m = monte carlo size

%Generating parameter_t



cov = generate_cov_row(x,par);%Generate covariance row for h_hat
circ = [cov; cov(end-1:-1:2)];  %Generate circulant row (to simulate)
lambda = fft(circ);

Z = [];
H = [];

for i=[1:m]
    h = (rand(1,n+1)<0.5).*2.-1;  
    Z_sim = conjgrad_circ(lambda/(2*n), h', n ,100);
    Z = [Z Z_sim]; 
    H = [H h'];
end


%generating circ row
col = generate_cov_row_d1(x,par);
circ = [col; col(end-1:-1:2)];
prods = [];

lambda = ifft(circ);
lambda = lambda*2*n;



for i = 1:size(H,2)
    r = [Z(:,i); zeros(n-1,1)];
    comp_r = [H(:,i); zeros(n-1,1)];
    %comp_z'*(toeplitz(circ))*comp_z
    quad = fft((lambda).*ifft(r));
    quad = comp_r'*quad;
    prods = [prods real(quad)];
    %prods = [prods quadratic_product(lambda,comp_z,1)/n];
end


%trace(toeplitz(col)*inv(toeplitz(cov)))



%var_prods = var(prods)

estimate = mean(prods)












s1 = Z;
s2 = 2;