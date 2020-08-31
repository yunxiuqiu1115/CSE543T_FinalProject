function [nlZ, dnlZ] = hypposterior(hyp, hyperparameters, ...
          inference_method, mean_function, covariance_function, ...
          likelihood, xs, ys, parms)
      num_samples = size(xs,1); 
      tmp.cov = [hyp(1:2);log(500);log(0.1)];
      tmp.lik = hyp(3);
      disp(hyp);
      hyp = tmp;
      means = reshape(hyperparameters.mean, num_samples,2);
      prior = inference_method{3};
      prior = rmfield(prior, 'mean');
      [nlZ, dnlZ] = gp_independent(hyp, prior, means, inference_method, ...
          mean_function, covariance_function, likelihood, xs, ys);
      nlZ = -nlZ;
      dnlZ = -1*unwrap(dnlZ);
      tmp = [dnlZ(1:2);dnlZ(5)];
      dnlZ = tmp;
%       prior = inference_method{3};
%       fn = fieldnames(hyp);
%       for i=1:numel(fn)
%         nhyp = numel(hyp.(fn{i}));
%         for j=1:nhyp
%             dist = prior.(fn{i}){j};
%             x = hyp.(fn{i}){j};
%             [a, b] = feval(dist{:}, x);
%             nlZ = nlZ+a;
%             dnlZ.(fn{i}){j} = dnlZ.(fn{i}){j}+b;
%         end
%      end
end