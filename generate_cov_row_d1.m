function C = generate_cov_row_dl(x,par)


% k = @(x) exp(-((x/par(1))^par(2)));
%k = @(x) par(2)*exp(-((x/par(1))));

k = @(x) (par(2)*x*exp(-((x/par(1)))))/par(1)^2;
n = length(x);
C = zeros(n,1);

j = 1;

for i = 1:n
    C(i,j) = k(abs(x(i)-x(j)));
end


end