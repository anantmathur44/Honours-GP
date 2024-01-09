
clear all
n=100; a=0; b=1;
a = 1/(n+1);
x=linspace(a,b,n+1);
theta = linspace(0.1-0.05,0.1,n+1)
par = [0.1 1];
hutch = []
bot = []
for i = theta
    cov = generate_cov_row(x,[i, 1]);
    col = generate_cov_row_d1(x,[i, 1]);
    a = toeplitz(col)*inv(toeplitz(cov));
    hutch = [hutch 2*trace(a^2)-2*sum(diag(a).^2)];
    bot = [bot trace(a^2) - trace(a)^2/(n+1)];
end

figure()
plot(theta,bot)
hold on;
plot(theta,hutch)