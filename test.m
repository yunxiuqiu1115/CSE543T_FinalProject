rng('default');
x = (1:5)';
y = randn(size(x));
xs = (6:7)';
ys = randn(size(xs));
mean = {@meanConst};
inf = @infExact;
cov = {@covSEiso};
lik = @likGauss;

theta.mean = 1;
theta.cov  = [log(1); log(1)];
theta.lik  = log(0.01);
[~, nlZ, dnlZ] = gp_last(theta, inf, mean, cov, lik, x, y, xs, ys);
d = 1e-6;
new_theta = theta;
new_theta.mean = theta.mean + d;
[~, new_nlZ] = gp_last(new_theta, inf, mean, cov, lik, x, y, xs, ys);
fprintf('dnlZ.mean: %0.3f, ((new_nlZ - nlZ) / d): %0.3f\n', dnlZ.mean, (new_nlZ - nlZ) / d);

d = 1e-6;
new_theta = theta;
new_theta.cov(1) = theta.cov(1) + d;
[~, new_nlZ] = gp_last(new_theta, inf, mean, cov, lik, x, y, xs, ys);
fprintf('dnlZ.cov(1): %0.3f, ((new_nlZ - nlZ) / d): %0.3f\n', dnlZ.cov(1), (new_nlZ - nlZ) / d);

d = 1e-6;
new_theta = theta;
new_theta.cov(2) = theta.cov(2) + d;
[~, new_nlZ] = gp_last(new_theta, inf, mean, cov, lik, x, y, xs, ys);
fprintf('dnlZ.cov(2): %0.3f, ((new_nlZ - nlZ) / d): %0.3f\n', dnlZ.cov(2), (new_nlZ - nlZ) / d);

d = 1e-6;
new_theta = theta;
new_theta.lik = theta.lik + d;
[~, new_nlZ] = gp_last(new_theta, inf, mean, cov, lik, x, y, xs, ys);
fprintf('dnlZ.lik: %0.4f, ((new_nlZ - nlZ) / d): %0.4f\n', dnlZ.lik, (new_nlZ - nlZ) / d);
