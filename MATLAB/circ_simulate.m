function Y = circ_simulate(x,c,mu,sigma_sq)

c = c.*sigma_sq;
n = length(x)-1;
lambda=fft(c);
if min(lambda) < 0
    fprintf('\n:NEGATIVE EIGENVALUE %d\n',min(lambda))
end
eta=sqrt(lambda./(2*n));
Z=randn(2*n,1)+sqrt(-1).*randn(2*n,1);
Zeta= Z.*eta;
X2n=fft(Zeta');
A=X2n;
Y=real(A);
Y  = Y + mu;

end
