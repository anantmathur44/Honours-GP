n= 10; a=0; b=1;
x=linspace(a,b,n+1);
%% Genereating covariance from kernal 
h =  [0.5 1];
cov = generate_cov_row(x,h);
Z = [];
circ = [cov; cov(end-1:-1:2)];

for i=[1:5000]    
    Z_sim = circ_simulate(x,circ,5,1);
    Z = [Z Z_sim'];       
end

prods = [];
prods2 = [];



for i=Z 
    prods = [prods i'*M*i] ;
    prods2 = [prods2 (i-5)'*M*(i-5)] ;
    
end
var(prods)
var(prods2)

T = (M + M')/2;

CC = T*toeplitz(circ);

variance = [var(prods2) var(prods)]
diff_act = var(prods) - var(prods2)
diff = 4*5*ones(2*n,1)'*CC*T*5*ones(2*n,1)
trace_term = 2*trace(CC*CC)
diff+trace_term


%  
% b = Z_sim(1:n+1);
% 
% 
% 
% norm(levinson_durbin(cov,b) - (toeplitz(cov))\b)
% 
% norm(conjgrad_circ(ifft(circ),b,n,100) - (toeplitz(cov))\b)
% 
% 
% % 
% % d = b - toeplitz(cov)*x_;
% % 
% % c_ = inv(toeplitz(circ));
% % c_ = c_(1:n+1,1:n+1);
% % 
% % p_ = c_*(d);
% % 
% % p_'*toeplitz(cov)*p_

