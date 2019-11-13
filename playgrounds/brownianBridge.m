% brownian bridge
n = 1000;
x = linspace(0,1,n);
muFunc = @(x) 0;
mu = arrayfun(muFunc, x);
kernel = @brownianBridgeKernel;

covariance = kernelMatrix(x, x, kernel);
figure(1);
myplot(x, mu, covariance, 3);
hold off;