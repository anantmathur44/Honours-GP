function [L,b]=LGL(gam)
% factorizes L*G*L'=diag(b), where G=toeplitz(gam) and 
% gam(1)*[1,rho(1:end)'], rho=gam(2:end)/gam(1);
% is the first row of the autocovariance matrix G,
% using Durbin's algorithm.



n=length(gam);
rho=gam(2:end)/gam(1);
L=diag(ones(n,1));
f=zeros(n,1);f(1)=rho(1);b(1)=1;
for k=1:n-1
    rk=rho(1:k);
    b(k+1)=1-rk'*f(1:k);
    fr=f(k:-1:1);
    L(k+1,1:k)=-fr;
    if k<n-1
        a=(rho(k+1)-rk'*fr)/b(k+1);
        f(1:k)=f(1:k)-a*fr;
        f(k+1)=a;
    end
end
b=b*gam(1);