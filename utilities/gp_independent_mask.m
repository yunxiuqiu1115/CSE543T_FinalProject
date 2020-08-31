function [nlZ, dnlZ] = gp_independent_mask(hyp, hyperparameters, ...
          inference_method, mean_function, covariance_function, ...
          likelihood, xs, ys, mask, parms)

  num_samples = size(xs,1);
  nfirm = size(hyperparameters.mean, 1) - 2*num_samples;
  
  infP = inference_method{1}; 
  inffunc = inference_method{2}; 
  prior = inference_method{3};
  
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
  N = 0;
  for i = 1:num_samples
    if isempty(xs{i}), continue; end
    hyp_race = full2one(hyperparameters, i, num_samples, nfirm);
    unsharedflag = false(size(unwrap(hyp_race)));
    unsharedflag(covsize+liksize) = 1;
    unsharedflag(1:2) = 1;
    
    im = {infP, inffunc, prior};
    N = N + size(xs{i},1);
        
    if (~gradient)
        if i==1
            this_nlZ = gp_mask(hyp, hyp_race, im, mean_function, ...
             covariance_function, likelihood, xs{i}, ys{i}, unsharedflag, parms, i, parms.mode);
        else
            this_nlZ = gp_mask(hyp, hyp_race, inffunc, mean_function, ...
             covariance_function, likelihood, xs{i}, ys{i}, unsharedflag, parms, i, parms.mode);
        end
      
    else
        if i==1 
            [this_nlZ, this_dnlZ] = gp_mask(hyp, hyp_race, im,...
          mean_function, covariance_function, likelihood, xs{i}, ys{i}, unsharedflag, parms, i, parms.mode);
        else           
            [this_nlZ, this_dnlZ] = gp_mask(hyp, hyp_race, inffunc,...
          mean_function, covariance_function, likelihood, xs{i}, ys{i}, unsharedflag, parms, i, parms.mode);
        end
      tmp = unwrap(this_dnlZ);
%       /size(xs{i},1)
      dnlZs{i} = rewrap(this_dnlZ, tmp);
    end

    % accumulate likelihoods and derivatives
    nlZ = nlZ + this_nlZ;
%     /num_samples
  end
  nlZ = nlZ/N;
  
  if (gradient)
      tmp = unwrap(dnlZ);
      for i=1:num_samples
          if isempty(dnlZs{i}), continue; end
          tmp = tmp + unwrap(dnlZs{i})/N;
%           /num_samples
      end
      dnlZ = rewrap(dnlZ, tmp);
  end
end