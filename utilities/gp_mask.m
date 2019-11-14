function [nlZ, dnlZ] = gp_mask(hyp, hyperparameters,...
          inference_method, mean_function, covariance_function, ...
          likelihood, x, y, mask)
      
  u = unwrap(hyperparameters);
  u(mask) = unwrap(hyp);
  hyperparameters = rewrap(hyperparameters, u);
  if (nargout == 1)
    nlZ = gp(hyperparameters, inference_method, mean_function, ...
           covariance_function, likelihood, x, y);
  else
    [nlZ, dnlZ] = gp(hyperparameters, inference_method, mean_function, ...
           covariance_function, likelihood, x, y);
  end

  if (nargout == 2)
    u = unwrap(dnlZ);
    u = u(mask);
    dnlZ = rewrap(hyp, u);
  end

end