function l  = likelihood_tr(x,Y,p,n,m)

cov = generate_cov_row(x,p);%Generate covariance row for h_hat
circ = [cov; cov(end-1:-1:2)]; %Generate circulant row (to simulate)
lambda = fft(circ);
circ1 = log(circ);

Z = [];

for i=[1:m]
    Z_sim = circ_simulate_seed(x,circ1,0,1,i);      
    Z = [Z Z_sim];       
end

length(Y)
length(lambda)
Z_c = conjgrad_circ(lambda/(2*n), Y, n ,100);
qp = Y*Z_c;


for i = Z
    quad = i'*i;
    prods = [prods quad];
    %prods = [prods quadratic_product(lambda,comp_z,1)/n];
end

l = -0.5*sum(prods)+0.5*qp;