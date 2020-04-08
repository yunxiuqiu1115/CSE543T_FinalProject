function [nlZ, dnlZ] = gp_mask(hyp, hyperparameters,...
          inference_method, mean_function, covariance_function, ...
          likelihood, x, y, mask, parms, i)
  if numel(x)==0
     nlZ = 0;
     if (nargout == 2)
        dnlZ = rewrap(hyperparameters, unwrap(hyperparameters)*0);
     end
     return 
  end
  u = unwrap(hyperparameters);
  u(mask) = unwrap(hyp);
  hyperparameters = rewrap(hyperparameters, u);
  xstar = [0,0,1,parms.nfirm,x(1,end)];
  if (nargout == 1)
    if (parms.mode == "last")
        ystar = parms.vs(i);
         nlZ = gp_last(hyperparameters, inference_method, mean_function, ...
           covariance_function, likelihood, x, y, xstar, ystar);
    else
        nlZ = gp(hyperparameters, inference_method, mean_function, ...
           covariance_function, likelihood, x, y);
    end
  else
    if (parms.mode == "last")
        ystar = parms.vs(i);
        [nlZ, dnlZ] = gp_last(hyperparameters, inference_method, mean_function, ...
           covariance_function, likelihood, x, y, xstar, ystar);
    else
        [nlZ, dnlZ] = gp(hyperparameters, inference_method, mean_function, ...
           covariance_function, likelihood, x, y);
    end
  end

  if (nargout == 2)
    u = unwrap(dnlZ);
    u = u(mask);
    dnlZ = rewrap(hyp, u);
  end

end