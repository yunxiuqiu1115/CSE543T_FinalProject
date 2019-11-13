% wiener process
n = 1000;
x = linspace(0,10,n);
muFunc = @(x) 0;
mu = arrayfun(muFunc, x);
kernel = @wienerKernel;

covariance = kernelMatrix(x, x, kernel);
figure(1);
myplot(x, mu, covariance, 3);
hold off;