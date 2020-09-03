function [allRaces,fts,s2s] = lm(besthyp, xs, ys, raceinfos, plot_path, parms)
    
    % iterate cycle/state race
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
            fts(i) = computePrior(pvi, experienced, republican, parms);
            s2s(i) = 0.1;
        else 
            SIGMA_inv = [1/0.002^2, 0; 0, 1/0.1^2];
            SIGMA = [0.002^2, 0; 0, 0.1^2];
            noise = exp(2*besthyp.lik);
            mu_b = [0;computePrior(pvi, experienced, republican, parms)];
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
            V = x*SIGMA*x'+LAMBDA+diag(noise*ones(size(ts,1),1));
            mu_p = mu_b + SIGMA*x'*(V\(y-x*mu_b));
            s2_p = SIGMA - SIGMA*x'*(V\(x*SIGMA));
            
            fmu = xstar*mu_p;
            fs2 = diag(xstar*s2_p*xstar');
            fts(i) = fmu(end);
            s2s(i) = fs2(end);
            
            
%             plot_posterior(fmu, fs2, x(:,1), y, linspace(x(1,1),0,nz).', 0, i);

            x = xs{i};
            y = ys{i};
            xstar = [linspace(x(1,1),0,nz).',zeros(1,nz)',ones(1,nz)',...
                            parms.nfirm*ones(1,nz)',republican*zeros(1,nz)'];
                        
            [meanfunc, covfunc, likfunc, inffunc, prior] = model(parms);
            covfunc{2} = covfunc{2}(2:end);
            hyp.mean(1) = 0;
            hyp.mean(2) = computePrior(pvi, experienced, republican, parms);
            hyp.cov = besthyp.cov(3:4);
            hyp.lik = besthyp.lik;
            
            [~,~, fmu, fs2] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, x, y, xstar);
            
            predPoll = fmu(end);
            fts(i) = predPoll;
            s2s(i) = fs2(end);
%             parms.prior = [computePrior(pvi, experienced, republican, parms), 0.1];
%             fig = plot_posterior(fmu, fs2, xs{i}(:,1), ys{i}, xstar(:,1), trueVote/100, parms);
% 
%             plot_title = year + " " + state + " " + candidateName;
%             title(plot_title);
%             yearFolder = fullfile(plot_path, num2str(year));
%             stateFolder = fullfile(yearFolder, state);
%             if ~exist(plot_path, 'dir')
%                 mkdir(plot_path);
%             end
%             if ~exist(yearFolder, 'dir')
%                 mkdir(yearFolder);
%             end
%             if ~exist(stateFolder, 'dir')
%                 mkdir(stateFolder);
%             end
%             filename = fullfile(stateFolder, plot_title + num2str(parms.j) + ".jpg");
%             saveas(fig, filename);
%             close;
            
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