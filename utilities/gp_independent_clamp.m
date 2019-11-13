function [nlZ, dnlZ] = gp_independent_clamp(hyperparameters, ...
          inference_method, mean_function, covariance_function, ...
          likelihood, xs, ys)

  num_samples = numel(xs);
  
  infP = inference_method{1}; 
  inffunc = inference_method{2}; 
  priors = inference_method{3};
  nfirm = size(hyperparameters.mean, 1) - 2*num_samples;

  % initialize nlZ and, optionally, dnlZ struct
  nlZ = 0;
  if (nargout == 2)
    dnlZ = rewrap(hyperparameters, 0 * unwrap(hyperparameters));
  end

  for i = 1:num_samples
    id = xs{i}(1,6);
    im = {infP, inffunc, priors{id}};
    hyp = full2one(hyperparameters, id, num_samples, nfirm);
    unsharedflag = false(size(unwrap(hyp)));
    unsharedflag(1:2) = 1;
    if (nargout == 1)
      this_nlZ = ...
          gp_clamp(hyp, im, mean_function, ...
             covariance_function, likelihood, xs{i}, ys{i}, unsharedflag);
    else
      [this_nlZ, this_dnlZ] = ...
          gp_clamp(hyp, im, mean_function, ...
             covariance_function, likelihood, xs{i}, ys{i}, unsharedflag);
        this_dnlZ = dnlZ2full(this_dnlZ, id, num_samples, nfirm);
    end

    % accumulate likelihoods and derivatives
    nlZ = nlZ + this_nlZ;
    if (nargout == 2)
      dnlZ = rewrap(dnlZ, unwrap(dnlZ) + unwrap(this_dnlZ));
    end
  end

end