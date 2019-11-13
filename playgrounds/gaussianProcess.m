% gaussian process Figure 2.2
n = 1000;
x = linspace(-5,5,n);
muFunc = @(x) 0;
numData = 10;
mu = arrayfun(muFunc, x);
noise = 0.0;
kernel = @squareExponentialKernel;
data = [-4,-3,-1,0,2];
y = [-2,0,1,2,-1];
% data = -5 + rand(1,numData)*10;
% y = mvnrnd(linspace(0,0,numData), kernelMatrix(data, data, kernel), 1);

covariance = kernelMatrix(x, x, kernel);
figure(1);
myplot(x, mu, covariance, 3);
plot(data, y.', 'x', 'MarkerSize', 20);
hold off;

[mu, covariance] = bayesianRegression(data, y, kernel, noise, x);
figure(2);
myplot(x, mu, covariance, 3);
plot(data, y.', 'x', 'MarkerSize', 20);
hold off;
