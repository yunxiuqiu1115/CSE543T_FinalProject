function [hyperparameters, nlZ] = fixLearn(hyperparameters, ...
          inference_method, mean_function, covariance_function, ...
          likelihood, xs, ys, iter, parms)

% flag = True if this hyp is fixed      
  num_samples = size(xs,1);
  as = parms.a;
  
  infP = inference_method{1}; 
  inffunc = inference_method{2}; 
  prior = inference_method{3};
  p.method = 'LBFGS';
  p.mem = 100;
  mfun = @minimize_v2;
  liksize = size(hyperparameters.lik, 1);
  covsize = size(hyperparameters.cov,1);
  nfirm = size(hyperparameters.mean,1) - 2*num_samples;
  
  for it=1:iter % repeat
      % fix shared parameters
      % loop for each race
      % learn unshared parameters independently
      disp("fix shared parameters...")
      p.verbosity = 0;
      p.length = -100;
      trained_hyp = zeros(num_samples, 2);
      parfor i = 1:num_samples
        if mod(i, 50) == 0, disp("training iter" + it + " trainning sample: " + i); end
        hyp = struct;
        hyp_race = full2one(hyperparameters, i, num_samples, nfirm);
        hyp.mean(1) = hyperparameters.mean(i);
        hyp.mean(2) = hyperparameters.mean(i+num_samples);
        if isempty(xs{i}), trained_hyp(i,:) = hyp.mean; continue; end
        mask = false(size(unwrap(hyp_race)));
        mask(liksize+covsize+1) = 1;
        mask(liksize+covsize+2) = 1;
        im = {infP, inffunc, prior};
        im{3}.mean{2}{2} = as(i);
        [hyp] = feval(mfun, hyp, @gp_mask, p, hyp_race, im,...
            mean_function, covariance_function,...
            likelihood, xs{i}, ys{i}, mask, parms, i, "all");
        trained_hyp(i,:) = hyp.mean;
      end
      hyperparameters.mean(1:num_samples) = trained_hyp(:,1);
      hyperparameters.mean(num_samples+1:2*num_samples) = trained_hyp(:,2);
      
      % fix unshared parameters
      % learn shared parameters independently
      disp("learn shared parameters");
%       parms.mode = 'last';
      p.verbosity = 1;
      p.length = -20;
      hyp = hyperparameters;
      hyp.mean = hyperparameters.mean(2*num_samples+1:end);
      mask = false(size(unwrap(hyperparameters)));
      mask(1:liksize+covsize) = 1;
      mask(end-nfirm+1:end) = 1;
      [hyp,nlZ] = feval(mfun, hyp, @gp_independent_mask, p, hyperparameters,...
          inference_method, mean_function, covariance_function,...
          likelihood, xs, ys, mask, parms);
      u = unwrap(hyperparameters);
      u(mask) = unwrap(hyp);
      hyperparameters = rewrap(hyperparameters, u);
  end
end