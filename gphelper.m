function [ymuvec ys2vec] = gphelper(hyp, inf, mean, cov, lik, xarray, yarray, xpredict)
    ymuvec = [];
    ys2vec = [];
    for i=1:numel(xarray)
        if numel(xarray{i}(:,1)) == 0
            continue;
        end
%         [nlZ, dnlZ] = gp(hyp, @infExact, mean, cov, lik, xarray{i}(:,1), yarray{i});
        [ymu ys2 fmu fs2] = gp(hyp, @infGaussLik, mean, cov, lik, xarray{i}(:,1), yarray{i}, xpredict{i}(:,1));
%         count = count + 1
    end
%     nlZsum = nlZsum / count;
%     dnlZsum = struct('mean', [], 'cov', [dnlZsumvec(1) / count;dnlZsumvec(2) / count], 'lik', dnlZsumvec(3));
end
