function score= mc_score_exact(x,Y,n,par,m)
%x:x-values, Y:Observed data, h_hat:starting par, m = monte carlo size

%Generating parameter_t


cov = generate_cov_row(x,par);%Generate covariance row for h_hat
circ = [cov; cov(end-1:-1:2)];  %Generate circulant row (to simulate)
lambda1 = fft(circ);


%generating circ row
col = generate_cov_row_d2(x,par);
circ = [col; col(end-1:-1:2)];


lambda = ifft(circ);
lambda = lambda*2*n;





ky = conjgrad_circ(lambda1/(2*n), Y, n ,100);
ky = [ky; zeros(n-1,1)];
quad = fft((lambda).*ifft(ky));
quad = ky'*quad;


score = -0.5*trace(toeplitz(col)*inv(toeplitz(cov))) + 0.5*real(quad);


