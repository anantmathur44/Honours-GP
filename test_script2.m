
n=1000; a=0; b=1;
x=linspace(a,b,n+1);
mu = 0;
sigma_sq = 1;
%% Genereating covariance from kernal 
h = [0.0395 2.2635];
h_1 = [0.0395    2.2609];
cov = generate_cov_row(x,h);
circ = [cov; cov(end-1:-1:2)];


Y = circ_simulate(x,circ,mu,sigma_sq);
Y = Y(1:n+1)';

circ1 = circ;
plot(Y);


lambda = fft(circ);
Z = []

m = 5

sol = conjgrad_circ(lambda/(2*n),Y, n ,100);    
padded_x = [sol ; zeros(n-1,1)];
b = matrix_vector_multi(padded_x,lambda);
mode = b(n+2:end);

for i=[1:m]
    Z_sim = circ_simulate(x,circ,0,1);
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

Z;

[c_oo,c_uo,c_uu] = circ_partition(circ,length(x));
cond_var = (c_uu-c_uo*inv(c_oo)*c_uo');
cond_varf = zeros(2*n);
cond_varf(n+2:end,n+2:end) = cond_var;

h = h_1
cov = generate_cov_row(x,h);
circ = [cov; cov(end-1:-1:2)];
lambda = ifft(circ);
mode = [Y;mode];
quad = sqrt(lambda).^-1.*ifft(mode);
quad_term = (norm(quad)^2)
plot(mode)
prods = [];
est = zeros(2*n);


for i = Z
    comp_z = [Y;i];
    %comp_z'*inv(toeplitz(circ))*comp_z
    %est = est + (comp_z-mode)*(comp_z-mode)';
    quad = sqrt(lambda).^-1.*ifft(comp_z-mode);
    prods = [prods (norm(quad)^2)];
    %prods = [prods quadratic_product(lambda,comp_z,1)/n];
end



prods2 = [];
for i = Z
    comp_z = [Y;i];
%     comp_z'*inv(toeplitz(circ))*comp_z
    quad = sqrt(lambda).^-1.*ifft(comp_z);
    prods2 = [prods2 (norm(quad)^2)];
    %prods = [prods quadratic_product(lambda,comp_z,1)/n];
end



lambda = ifft(circ1);
prods_3 = [];


for i = Z
    comp_z = [Y;i];
    quad = sqrt(lambda).^-1.*ifft(comp_z-mode);
    prods_3 = [prods_3 (norm(quad)^2)];
    %prods = [prods quadratic_product(lambda,comp_z,1)/n];
end




A = inv(toeplitz(circ));
trace_var = 2*trace((A*cond_varf)^2)


trace_term = trace(A*cond_varf);


trace_term;
mean(prods-prods_3+(n-1));


var_2 = (prods-prods_3+(n-1));

% 
var_2;
prods;
prods2;
trace_term


var_new = mean(var_2)
var_o1d = mean(prods)
unobserved_dim_est = mean(prods_3)

var_red_2 = sum((var_2 - trace_term).^2)/m
var_red = sum((prods - trace_term).^2)/m
var_norm = sum((prods2 - ((trace_term)+ quad_term)).^2)/m

emp_diff_var = var_norm/var_red;
diff_var = (var_red+(4*mode'*A*cond_varf*A*mode))/var_red;

red = (sum(prods)/m);

trace_term;

normal = (sum(prods2)/m)

reduced = red + quad_term

reduced_2 = reduced + ((-mean(prods_3)+(n-1)))

zero_term = ((-mean(prods_3)+(n-1)))


red-trace_term


actual = trace_term+quad_term

figure(1);
plot(Y);
figure(2);
plot(mode);

