function hyperparameters = fixLearn(hyperparameters, ...
          inference_method, mean_function, covariance_function, ...
          likelihood, xs, ys, iter)

% flag = True if this hyp is fixed      
  num_samples = numel(xs);
  
%   infP = inference_method{1}; 
%   inffunc = inference_method{2}; 
%   priors = inference_method{3};
  p.method = 'LBFGS';
  p.mem = 100;
  p.length = -100;
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
      trained_hyp = zeros(num_samples, 2);
      disp("training iter" + it);
      parfor i = 1:num_samples
        % if mod(i, 10) == 0, disp("trainning sample: " + i); end
        % disp("trainning sample: " + i);
        % im = {infP, inffunc, priors{i}};
        hyp = struct;
        hyp_race = full2one(hyperparameters, i, num_samples, nfirm);
        hyp.mean(1) = hyperparameters.mean(i);
        hyp.mean(2) = hyperparameters.mean(i+num_samples);
        mask = false(size(unwrap(hyp_race)));
        mask(liksize+covsize+1) = 1;
        mask(liksize+covsize+2) = 1;
        hyp = feval(mfun, hyp, @gp_mask, p, hyp_race, inference_method,...
            mean_function, covariance_function,...
            likelihood, xs{i}, ys{i}, mask);
        trained_hyp(i,:) = hyp.mean;
      end
      hyperparameters.mean(1:num_samples) = trained_hyp(:,1);
      hyperparameters.mean(num_samples+1:2*num_samples) = trained_hyp(:,2);
      
      % fix unshared parameters
      % learn shared parameters independently
      p.verbosity = 1;
%       p.length = -20;
      hyp = hyperparameters;
      hyp.mean = hyperparameters.mean(2*num_samples+1:end);
      mask = false(size(unwrap(hyperparameters)));
      mask(1:liksize+covsize) = 1;
      mask(end-nfirm+1:end) = 1;
      hyp = feval(mfun, hyp, @gp_independent_mask, p, hyperparameters,...
          inference_method, mean_function, covariance_function,...
          likelihood, xs, ys, mask);
      u = unwrap(hyperparameters);
      u(mask) = unwrap(hyp);
      hyperparameters = rewrap(hyperparameters, u);
  end
end