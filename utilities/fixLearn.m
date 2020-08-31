function [hyperparameters, nlZs] = fixLearn(hyperparameters, ...
          inference_method, mean_function, covariance_function, ...
          likelihood, xs, ys, iter, parms)

% flag = True if this hyp is fixed      
  num_samples = size(xs,1);
  nlZs = [];
  
  infP = inference_method{1}; 
  inffunc = inference_method{2}; 
  prior = inference_method{3};
  p.method = 'LBFGS';
  p.mem = 100;
  mfun = @minimize;
  liksize = size(hyperparameters.lik, 1);
  covsize = size(hyperparameters.cov,1);
%   nfirm = size(hyperparameters.mean,1) - 2*num_samples;
  
  p.verbosity = 0;
  p.length = -iter;
  hyp.cov = hyperparameters.cov;
  hyp.lik = hyperparameters.lik;
  mask = false(size(unwrap(hyperparameters)));
  mask(1:4) = 1;
  mask(liksize+covsize) = 1;
%   mask(end-nfirm+1:end) = 1;
  i = parms.i;
  means = reshape(hyperparameters.mean, num_samples,2);
%   means = means(i,:);
  prior = rmfield(prior, 'mean');
  [hyp,nlZ] = feval(mfun, hyp, @gp_independent, -p.length, prior, means, inffunc, ...
          mean_function, covariance_function, likelihood, xs, ys);
%   [hyp,nlZ] = feval(mfun, hyp, @gp_independent_mask, p, hyperparameters,...
%       inference_method, mean_function, covariance_function,...
%       likelihood, xs, ys, mask, parms);
  nlZs = [nlZs,nlZ];
  u = unwrap(hyperparameters);
  u(mask) = unwrap(hyp);
  hyperparameters = rewrap(hyperparameters, u);
end