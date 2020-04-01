function [post, nlZ, dnlZ] = gp_last(hyp, inf, mean, cov, lik, x, y, xs, ys)
    [m, v, ~, ~, lp] = gp(hyp, inf, mean, cov, lik, x, y, xs, ys);
    nlZ = -sum(lp);
    sn2 = exp(2*hyp.lik);
    ns = size(xs,1);
    % Vinv = diag(1./v);
    
    if sn2<1e-6                        % very tiny sn2 can lead to numerical trouble
      L = chol(v+sn2*eye(ns)); sl = 1;   % Cholesky factor of covariance with noise
      pL = -solve_chol(L,eye(ns));                            % L = -inv(K+inv(sW^2))
    else
      L = chol(v/sn2+eye(ns)); sl = sn2;                       % Cholesky factor of B
      pL = L;                                           % L = chol(eye(n)+sW*sW'.*K)
    end
    alpha = solve_chol(L,ys-m)/sl;
    % alpha = Vinv*(ys-m);
    post.alpha = alpha;
    post.sW = ones(n,1)/sqrt(v);
    post.L = pL;

    dnlZ = hyp;
    Q = solve_chol(L,eye(ns))/sl - alpha*alpha';
    
    for i = 1:numel(hyp.mean)
      dnlZ.mean(i) = -feval(mean{:}, hyp.mean, xs, i)'*alpha;
    end
    
    for i = 1:numel(hyp.cov)
      % dKi = feval(cov{:}, hyp.cov, xs, [], i);
      % dnlZ.cov(i) = alpha.'*dKi*alpha/2 + trace(Vinv*dKi)/2;
      dnlZ.cov(i) = sum(sum(Q.*feval(cov{:}, hyp.cov, xs, [], i)))/2;
    end
    
    % Q = trace(Vinv) - alpha.'*alpha;
    dnlZ.lik = sn2*trace(Q);
end

