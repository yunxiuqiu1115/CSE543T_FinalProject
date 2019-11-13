function hyp = sample_prior(prior)
    fn = fieldnames(prior);
    for i=1:numel(fn)
        nhyp = numel(prior.(fn{i}));
        sampled_hyp = zeros(nhyp,1);
        for j=1:nhyp
            dist = prior.(fn{i}){j};
            sampled_hyp(j) = feval(dist{:});
        end
        hyp.(fn{i}) = sampled_hyp;
    end
end