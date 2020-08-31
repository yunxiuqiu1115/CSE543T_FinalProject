function hyp = sample_separate_prior(prior, parms, n, seed)
    rng("default");
    rng(seed);
%     hyps = cell(n, 1);
    
%     for i=1:n
%         prior.mean{2}{2} = parms.a(i);
%         hyps{i} = sample_prior(prior);
%     end
    
%     hyp = hyps{1};
    hyp.mean = zeros(2*n+parms.nfirm,1);
%     hyp.cov = zeros(2*n+parms.nfirm,1);
%     hyp.lik = zeros(n,1);
    for i=1:n
        hyp.mean(i) = 0;
        hyp.mean(i+n) = parms.a(i);
%         hyp.cov(i) = hyps{i}.cov(1);
%         hyp.cov(i+n) = hyps{i}.cov(2);
%         hyp.lik(i) = hyps{i}.lik;
    end
    
    sigma_ml = 0.002;
    sigma_mc = 0.1;
    
    hyp.cov = [feval(prior.cov{1}{:});feval(prior.cov{2}{:});-log(sigma_ml);log(sigma_mc)];
    hyp.lik = [feval(prior.lik{1}{:})];
    
    for i=1:parms.nfirm
       hyp.mean(i+2*n) = prior.mean{2+i}{2};
    end
end