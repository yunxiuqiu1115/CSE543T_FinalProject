function [nlZ, dnlZ] = gp_independent_mask(hyp, hyperparameters, ...
          inference_method, mean_function, covariance_function, ...
          likelihood, xs, ys, mask)

  num_samples = numel(xs);
%   infP = inference_method{1}; 
%   inffunc = inference_method{2}; 
%   priors = inference_method{3};
  nfirm = size(hyperparameters.mean, 1) - 2*num_samples;
  
  u = unwrap(hyperparameters);
  u(mask) = unwrap(hyp);
  hyperparameters = rewrap(hyperparameters, u);
  liksize = size(hyperparameters.lik, 1);
  covsize = size(hyperparameters.cov,1);

  % initialize nlZ and, optionally, dnlZ struct
  nlZ = 0;
  gradient = true;
  if nargout == 1, gradient = false; end
  if (gradient)
    dnlZ = rewrap(hyp, 0 * unwrap(hyp));
  end
  
  dnlZs = cell(num_samples, 1);
  parfor i = 1:num_samples
    % im = {infP, inffunc, priors{i}};
    hyp_race = full2one(hyperparameters, i, num_samples, nfirm);
    unsharedflag = true(size(unwrap(hyp_race)));
    unsharedflag(liksize+covsize+1) = 0;
    unsharedflag(liksize+covsize+2) = 0;
    if (~gradient)
      this_nlZ = gp_mask(hyp, hyp_race, inference_method, mean_function, ...
             covariance_function, likelihood, xs{i}, ys{i}, unsharedflag);
    else
      [this_nlZ, this_dnlZ] = gp_mask(hyp, hyp_race, inference_method,...
          mean_function, covariance_function, likelihood, xs{i}, ys{i}, unsharedflag);
    end

    % accumulate likelihoods and derivatives
    nlZ = nlZ + this_nlZ;
    if (gradient)
      dnlZs{i} = this_dnlZ;
    end
  end
  
  if (gradient)
      tmp = unwrap(dnlZ);
      for i=1:num_samples
          tmp = tmp + unwrap(dnlZs{i})/num_samples;
      end
      dnlZ = rewrap(dnlZ, tmp);
  end
end