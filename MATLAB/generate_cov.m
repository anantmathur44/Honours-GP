function C = generate_cov(x,par)


k = @(x) exp(-sqrt(x'*x)/par);
n = length(x);
C = zeros(n,1);

j = 1;

for i = 1:n
    C(i,j) = k(x(i)-x(j));
end

C = toeplitz(C);

