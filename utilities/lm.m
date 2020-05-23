function [allRaces,fts,s2s] = lm(xs, raceinfos, parms)
    % iterate cycle/state race
    [~,~,~,~, prior] = model(parms);
    allRaces = struct;
    n = numel(xs);
    fts = zeros(n,1);
    s2s = zeros(n,1);
    for i = 1:n
        year = raceinfos{i}{1};
        state = raceinfos{i}{2}{1};
        candidateName = raceinfos{i}{3};
        trueVote = raceinfos{i}{4};
        pvi =raceinfos{i}{5};
        experienced =raceinfos{i}{6};
        republican = raceinfos{i}{7};
        fn = char(state+""+year);
        fn = fn(~isspace(fn));
        if numel(xs{i})==0
            fts(i) = computePrior(pvi, experienced, republican);
            s2s(i) = prior.mean{2}{3};
        else 
            SIGMA_inv = [1/prior.mean{1}{3}, 0; 0, 1/prior.mean{2}{3}];
            mu_b = [prior.mean{1}{2};computePrior(pvi, experienced, republican)];
            x = xs{i};
            ps = x(:,2);
            idx = find(ps==0);
            if numel(idx)~=0
                x(idx, :) = [];
            end
            ps = x(:,2);
            ns = x(:,3);
            ts = x(:,1);
            LAMBDA_inv = diag(ns./ps./(1-ps));
            x = [ts, ones(size(ts,1),1)];
            M = x'*LAMBDA_inv*x+SIGMA_inv;
            b = x'*LAMBDA_inv*ps + SIGMA_inv*mu_b;
            xstar = [0, 1];
            fts(i) = xstar*(M\b);
            s2s(i) = xstar*(M\xstar');
        end
      
        if ~isfield(allRaces, fn)
            allRaces.(fn) = [fts(i), trueVote];
        else
            allRaces.(fn) = [allRaces.(fn), fts(i), trueVote];
        end
    end
end