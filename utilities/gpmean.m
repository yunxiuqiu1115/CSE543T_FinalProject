function [nlZsum, dnlZsum] = gpmean(hyp, inf, mean, cov, lik, xarray, yarray)
    nlZsum = 0.0;
    nlZvec = [];
    dnlZsumvec = [0;0;0;0];     
    counter = 1;
    for i=1:numel(xarray)
        if numel(xarray{i}(:,1)) == 0
            continue;
        end
        [nlZ, dnlZ] = gp(hyp, @infExact, mean, cov, lik, xarray{i}(:,1), yarray{i});
        nlZvec(counter) = nlZ;
        counter = counter + 1;
        nlZsum = nlZsum + nlZ;
        dnlZsumvec = dnlZsumvec + unwrap(dnlZ);
    end
    disp("nlZsum");disp(nlZsum);
    dnlZsumvec = dnlZsumvec / counter;
%      dnlZsum = struct('mean', [], 'cov', [dnlZsumvec(1);dnlZsumvec(2)], 'lik', dnlZsumvec(3));
%       dnlZsum = struct('mean', [], 'cov', [dnlZsumvec(1)], 'lik', dnlZsumvec(2));
      dnlZsum = struct('mean', [], 'cov', [dnlZsumvec(1); dnlZsumvec(2); dnlZsumvec(3)], 'lik', dnlZsumvec(4));
end