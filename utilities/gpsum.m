function [nlZsum, dnlZsum] = gpsum(hyp, inf, mean, cov, lik, xarray, yarray)
    nlZsum = 0.0;
    dnlZsumvec = [0;0;0];
%     count = 0;
    for i=1:numel(xarray)
        if numel(xarray{i}(:,1)) == 0
            continue;
        end
        [nlZ, dnlZ] = gp(hyp, @infExact, mean, cov, lik, xarray{i}(:,1), yarray{i});
%         count = count + 1;
        nlZsum = nlZsum + nlZ;
        dnlZsumvec = dnlZsumvec + unwrap(dnlZ);
    end
%     nlZsum = nlZsum / count;
    dnlZsum = struct('mean', [], 'cov', [dnlZsumvec(1);dnlZsumvec(2)], 'lik', dnlZsumvec(3));
%     dnlZsum = struct('mean', [], 'cov', [dnlZsumvec(1) / count;dnlZsumvec(2) / count], 'lik', dnlZsumvec(3));
end
