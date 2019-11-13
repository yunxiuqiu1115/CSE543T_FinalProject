function hyp = full2one(hyperparameters, id, counter, nfirm)
    hyp = hyperparameters;
    hyp.mean = zeros(2+nfirm,1);
%     hyp.cov = zeros(2+nfirm,1);
    hyp.mean(1) = hyperparameters.mean(id);
    hyp.mean(2) = hyperparameters.mean(id+counter);
%     hyp.cov(1) = hyperparameters.cov(id);
%     hyp.cov(2) = hyperparameters.cov(id+counter);
    for j=1:nfirm
        hyp.mean(2+j) = hyperparameters.mean(2*counter+j);
%         hyp.cov(2+j) = hyperparameters.cov(2*counter+j);
    end
%     hyp.lik = hyperparameters.lik(id);
end