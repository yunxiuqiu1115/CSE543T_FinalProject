% bayesian regression

function [ystar, covariancestar] = bayesianRegression(x, y, kernel, noise, xstar)
    n = length(x);
    Kstar = kernelMatrix(x, xstar, kernel);
    K = kernelMatrix(x, x, kernel);
    L = chol(K + noise*eye(n), 'lower');
    a = L.'\(L\y.');
    ystar = Kstar.'*a;
    v = L\Kstar;
    covariancestar = kernelMatrix(xstar, xstar, kernel) - v.'*v;
end
