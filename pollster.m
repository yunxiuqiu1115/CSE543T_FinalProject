% pollster biases
% pbs = hyp.mean(760*2+1:end);
% pstds = exp(hyp.cov(3:end));
function pbs = pollster(besthyp, xs, ys, raceinfos, parms)
    [meanfunc, covfunc, likfunc, inffunc, prior] = model(parms);
    im = {@infPrior, inffunc, prior};
    n = numel(xs);
    p.method = 'LBFGS';
    p.mem = 100;
    p.verbosity = 0;
    p.length = -100;
    mfun = @minimize_v2;
    pbs = cell(41,1);
    for i=1:41
       pbs{i} = []; 
    end

    for i=1:n
        pids = xs{i}(:,4);
        pvi =raceinfos{i}{5};
        experienced =raceinfos{i}{6};
        republican = raceinfos{i}{7};
        priorb = prior.mean{2};

        xstar = xs{i};
        if i<=760
            % training forecasting
            % im{3}.mean{2}{2} = parms.a(i);
            hyp = full2one(besthyp, i, parms.ncandidates, parms.nfirm);
        else
            % testing forecasting
            % im{3}.mean{2}{2} = computePrior(pvi, experienced, republican);
            hyp = full2one(besthyp, 1, parms.ncandidates, parms.nfirm);
            hyp.mean(1) = prior.mean{1}{2};
            priorb{2} = computePrior(pvi, experienced, republican);
            hyp.mean(2) = feval(priorb{:});
            hyp.mean(2) = priorb{2};
            im{3}.mean{2}{2} = priorb{2};
            parms.mode = "all";
            mask = false(size(unwrap(hyp)));
            liksize = size(hyp.lik, 1);
            covsize = size(hyp.cov,1);
            mask(liksize+covsize+1) = 1;
            mask(liksize+covsize+2) = 1;
            hypab.mean = hyp.mean(1:2);
            hypab = feval(mfun, hypab, @gp_mask, p, hyp, im,...
                    meanfunc, covfunc,...
                    likfunc, xs{i}, ys{i}, mask, parms, i, "all");
            hyp.mean(1:2) = hypab.mean;
        end

        xstar(:,5)=0;
        [~, ~, fmu, ~] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, xs{i}, ys{i}, xstar);

        for j=1:length(pids)
            pid = pids(j);
            b = ys{i}(j)-fmu(j);
            pbs{pid} = [pbs{pid},b];
        end
    end

    for i=1:41
       pbs{i} = sum(pbs{i})/length(pbs{i}); 
    end
    
    scatter(besthyp.mean(760*2+1:end),cell2mat(pbs)); xlabel('learnt pollster biases'); ylabel("calculated pollster biases");
end