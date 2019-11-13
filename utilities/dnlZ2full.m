function dnlZ_full = dnlZ2full(dnlZ, id, counter, nfirm)
    dnlZ_full = dnlZ;
    dnlZ_full.mean = zeros(2*counter + nfirm, 1);
%     dnlZ_full.cov = zeros(2*counter + nfirm, 1);
%     dnlZ_full.lik = zeros(counter, 1);

    for j=1:nfirm
        dnlZ_full.mean(2*counter+j) = dnlZ.mean(2+j);
%         dnlZ_full.cov(2*counter+j) = dnlZ.cov(2+j);
    end
    dnlZ_full.mean(id) = dnlZ.mean(1);
    dnlZ_full.mean(id + counter) = dnlZ.mean(2);
%     dnlZ_full.cov(id) = dnlZ.cov(1);
%     dnlZ_full.cov(id + counter) = dnlZ.cov(2);
%     dnlZ_full.lik(id) = dnlZ.lik;
end