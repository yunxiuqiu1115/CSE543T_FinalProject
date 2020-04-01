function [post, nlZ, dnlZ] = infLast(hyp, mean, cov, lik, x, y)
    xs = x(end,:);
    ys = y(end,:);
    [post, nlZ, dnlZ] = gp_last(hyp, @infExact, mean, cov, lik, x, y, xs, ys);
    if xs(1)<=-2*7
       nlZ = 0;
       dnlZ = rewrap(hyp, unwrap(dnlZ)*0);
    end
end

