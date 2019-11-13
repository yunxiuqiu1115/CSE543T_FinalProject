function hyperparameters = fixLearn(hyperparameters, ...
          inference_method, mean_function, covariance_function, ...
          likelihood, xs, ys, iter)

% flag = True if this hyp is fixed      
  num_samples = numel(xs);
  
  infP = inference_method{1}; 
  inffunc = inference_method{2}; 
  priors = inference_method{3};
  nfirm = size(hyperparameters.mean, 1) - 2*num_samples;

  mfun = @minimize_v2;
  p.method = 'LBFGS';
  p.length = -100;
  p.mem = 100;
  
  for it=1:iter % repeat
      % fix shared parameters
      % loop for each race
      % learn unshared parameters independently
      disp("fix shared parameters...")
      for i = 1:num_samples
        if mod(i, 10) == 0, disp("trainning sample: " + i); end
        p.verbosity = 0;
        id = xs{i}(1,6);
        im = {infP, inffunc, priors{id}};
        hyp = full2one(hyperparameters, id, num_samples, nfirm);
        sharedflag = true(size(unwrap(hyp)));
        sharedflag(1:2) = 0;
        hyp = feval(mfun, hyp, @gp_clamp, p, im,...
            mean_function, covariance_function, likelihood, xs{i}, ys{i}, sharedflag);
        hyperparameters.mean(id) = hyp.mean(1);
        hyperparameters.mean(id + num_samples) = hyp.mean(2);
      end
      
      % fix unshared parameters
      % learn shared parameters independently
      p.verbosity = 1;
      hyperparameters = feval(mfun, hyperparameters, @gp_independent_clamp, p,...
          inference_method, mean_function, covariance_function,...
          likelihood, xs, ys);
      
  end
end