% addpath("/Users/yahoo/Documents/WashU/CSE515T/Code/gpml-matlab-v4.2-2018-06-11");
addpath("/Users/yahoo/Documents/WashU/CSE515T/Code/gpml-matlab-v3.6-2015-07-07");
startup;
n = 1e2;
x = linspace(0,10,n).';

% define mean, kernel, likelihood and inference
mc = {@meanConst};
cse = {@covSEiso};
md = {@meanDiscrete, n};
ml = {@meanLinear};
msk = [true, false];
mask = [false, true];
meanfunc = {@meanMask, msk, mc};
cd = {@covDiag, md};
covfunc = cd;
likfunc = @likGauss;
inffunc = @infExact;

% define hyperparameters
offset       = 0;
% length_scale = 1;
% output_scale = 1;
noise = 0.1;
noise_vec = ones(n,1)*noise^2;
hyp.mean = [offset]; 
hyp.lik = log(noise);
hyp.cov = log(sqrt(noise_vec));
% hyp.cov = log([length_scale, output_scale]);

% generate fake data
mu = feval(meanfunc{:}, hyp.mean, x);
K_f = feval(covfunc{:}, hyp.cov, x);
K_y = K_f + diag(noise_vec);
y = mvnrnd(mu, K_y).';
xs = linspace(-1,11,n).';
% unknown xs does not have index
% so just set xs = x

[~, ~, fmu, fs2] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, x, y, xs);
plot_posterior(fmu, fs2, x, y, xs, 1);

% tuning hyperparameters
% hyp.mean = 1; 
% hyp.cov = log(1);
% hyp.lik = log(1);

pg = {@priorGauss, 0, 1}; pc = {@priorClamped};
prior.mean = {pg};
prior.cov  = {pg};
prior.lik  = {pg};
im = {@infPrior, inffunc, prior};
par = {im, meanfunc,covfunc,likfunc, x, y}; mfun = @minimize;
hyp = feval(mfun, hyp,  @gp, -100, par{:});

% perform regression
[~, ~, fmu, fs2] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, x, y, xs);
plot_posterior(fmu, fs2, x, y, xs, 2);

