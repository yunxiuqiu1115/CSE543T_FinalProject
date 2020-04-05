function [post, nlZ, dnlZ] = gp_last(hyp, inf, mean, cov, lik, x, y, xs, ys)
    % obtain ys_mu, ys_sn2, post
    [m, ~, ~, ~, ~, post] = gp(hyp, inf, mean, cov, lik, x, y, xs, ys);
    
    % define posterior derivative
    alpha = post.alpha;
    L = post.L;
    dnlZ = hyp;
    
    % log likelihood of predictive distribution
    kxs = feval(cov{:}, hyp.cov, xs);
    kxsx = feval(cov{:}, hyp.cov, xs, x);
    
    % get posterior variance
    if (is_chol(L))
      % posterior.L contains chol(sqrt(W) * K * sqrt(W) + I)
      tmp = bsxfun(@times, kxsx', post.sW)';
      b = solve_chol(L, tmp');
      pv = kxs - tmp * b;
    else
      % posterior.L contains -inv(K + inv(W))
      pv = kxs + kxsx * L * kxsx;
    end
    
    Ls = chol(pv);
    ns = size(xs,1);
%     m = feval(mean{:}, hyp.mean, xs) + feval(cov{:}, hyp.cov, xs, x)*alpha;
%     pv = kxs - kxsx*V_inv*kxsx';
    
    bs = solve_chol(Ls, (ys-m));
    nlZ = (ys-m).'*bs/2 + log(det(pv)) + ns*log(2*pi)/2;
%    nlZ = -sum(lp);
    
    
    % V = feval(cov{:}, hyp.cov, x)+exp(2*hyp.lik);
    % V_inv = inv(V);
%     bt = kxsx*V_inv;
%     bt = solve_chol(L, kxsx')';
    bt = bsxfun(@times, b, post.sW)';

    for i = 1:numel(hyp.mean)
      dmxi = feval(mean{:}, hyp.mean, x, i)';
      dmxsi = feval(mean{:}, hyp.mean, xs, i)';
      dnlZ.mean(i) = -bs'*(dmxsi'-bt*dmxi');
    end
    
    sn2 = exp(2*hyp.lik);
    dsms = -bt*alpha*2*sn2;
    dsvs = bt*bt'*2*sn2;
    Q = solve_chol(Ls, dsvs);
    dnlZ.lik = -bs'*dsms - bs'*dsvs*bs/2 + trace(Q);
    
    for i = 1:numel(hyp.cov)
      dKxi = feval(cov{:}, hyp.cov, x, x, i);
      dKxsxi = feval(cov{:}, hyp.cov, xs, x, i);
      dKxsi = feval(cov{:}, hyp.cov, xs, xs, i);
      dKms = (dKxsxi-bt*dKxi)*alpha;
      dvxs = dKxsi-bt*(2*dKxsxi'-dKxi*bt');
      Q = solve_chol(Ls, dvxs);
      dnlZ.cov(i) = -bs'*dKms - bs'*dvxs*bs/2 + trace(Q);
    end
    
end
