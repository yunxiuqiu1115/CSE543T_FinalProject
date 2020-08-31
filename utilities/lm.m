function [allRaces,fts,s2s] = lm(hyp, xs, ys, raceinfos, parms)
    % iterate cycle/state race
    allRaces = struct;
    n = numel(xs);
    fts = zeros(n,1);
    s2s = zeros(n,1);
    for i = 767:767
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
            s2s(i) = 0.1;
        else 
            SIGMA_inv = [1/0.002, 0; 0, 1/0.1];
            SIGMA = [0.002, 0; 0, 0.1];
            noise = exp(2*hyp.lik);
            mu_b = [0;computePrior(pvi, experienced, republican)];
            x = xs{i};
            y = ys{i};
            ps = x(:,2);
            idx = find(ps==0);
            if numel(idx)~=0
                x(idx, :) = [];
                y(idx, :) = [];
            end
            ps = x(:,2);
            ns = x(:,3);
            ts = x(:,1);
            LAMBDA_inv = diag(ns./ps./(1-ps));
            LAMBDA = diag(ps.*(1-ps)./ns);
            x = [ts, ones(size(ts,1),1)];
            
            nz = 200;
            xstar = [linspace(x(1,1),0,nz).', ones(nz,1)];
            V = x*SIGMA*x'+LAMBDA+noise;
            mu_p = mu_b + SIGMA*x'*(V\(y-x*mu_b));
            s2_p = SIGMA - SIGMA*x'*(V\(x*SIGMA));
            
            fmu =xstar*mu_p;
            fs2=xstar*s2_p*xstar';
            fts(i) = fmu(end);
            s2s(i) = fs2(end);
            
            plot_posterior(fmu, fs2, x(:,1), y, linspace(x(1,1),0,nz).', 0, i);
            
%             M = x'*LAMBDA_inv*x+SIGMA_inv + noise;
%             b = x'*LAMBDA_inv*ps + SIGMA_inv*mu_b;
            
%             fts(i) = xstar*(M\b);
%             s2s(i) = xstar*(M\xstar');
        end
      
        if ~isfield(allRaces, fn)
            allRaces.(fn) = [fts(i), trueVote];
        else
            allRaces.(fn) = [allRaces.(fn), fts(i), trueVote];
        end
    end
end