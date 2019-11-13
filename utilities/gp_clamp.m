function [nlZ, dnlZ] = gp_clamp(hyp, ...
          inference_method, mean_function, covariance_function, ...
          likelihood, x, y, clampflag)
  if (nargout == 1)
    nlZ = gp(hyp, inference_method, mean_function, ...
           covariance_function, likelihood, x, y);
  else
    [nlZ, dnlZ] = gp(hyp, inference_method, mean_function, ...
           covariance_function, likelihood, x, y);
  end

  if (nargout == 2)
    dnlZ = rewrap(dnlZ, ~clampflag.*unwrap(dnlZ));
  end

end