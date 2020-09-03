function [allRaces,fts,s2s] = forcastAllRaces(besthyp, xs, ys, raceinfos, plot_path, parms)
    % iterate cycle/state race
    [meanfunc, covfunc, likfunc, inffunc, prior] = model(parms);
    im = {@infPrior, inffunc, prior};
    allRaces = struct;
    n = numel(xs);
    fts = zeros(n,1);
    s2s = zeros(n,1);
    nz = 200;
    p.method = 'LBFGS';
    p.mem = 100;
    p.verbosity = 0;
    p.length = -100;
%     a_mu = sum(besthyp.mean(1:760))/760;
%     a_std = std(besthyp.mean(1:760));
%     b_std = std(besthyp.mean(761:760*2));
    mfun = @minimize_v2;
    for i = 1:n
        year = raceinfos{i}{1};
        state = raceinfos{i}{2}{1};
        candidateName = raceinfos{i}{3};
        trueVote = raceinfos{i}{4};
        pvi = raceinfos{i}{5};
        experienced = raceinfos{i}{6};
        republican = raceinfos{i}{7};
        fn = char(state+""+year);
        fn = fn(~isspace(fn));
        priorb = prior.mean{2};
        if numel(xs{i})==0
            % if there is no data avaiable
            % use prior
            if i<=760
                % training forecasting
                priorb{2} = besthyp.mean(parms.ncandidates+i);
                predPoll = besthyp.mean(parms.ncandidates+i);
            else 
                % testing forecasting
                priorb{2} = computePrior(pvi, experienced, republican, parms);
                % + besthyp.mean(end)*republican;
                predPoll = priorb{2};
            end
            if ~isfield(allRaces, fn)
                allRaces.(fn) = [predPoll, trueVote];
            else
                allRaces.(fn) = [allRaces.(fn), predPoll, trueVote];
            end
            fts(i) = predPoll;
            s2s(i) = exp(2*besthyp.lik);
        else
            % if there is new data avaiable
            % use MAP
            % republican = xs{i}(1,5);
            xstar = [linspace(xs{i}(1,1),0,nz).',zeros(1,nz)',ones(1,nz)',...
                parms.nfirm*ones(1,nz)',republican*ones(1,nz)'];
            if i<=760
                % training forecasting
                % im{3}.mean{2}{2} = parms.a(i);
%                 hyp = full2one(besthyp, i, parms.ncandidates, parms.nfirm);
                hyp.mean(1) = 0;
                hyp.mean(2) = computePrior(pvi, experienced, republican, parms);
                hyp.cov = besthyp.cov;
                hyp.lik = besthyp.lik;
%                 hyp.cov(1)=cov1(i);
%                 hyp.cov(2)=cov2(i);
%                 hyp.lik=lik1(i);
            else
                % testing forecasting
                % im{3}.mean{2}{2} = computePrior(pvi, experienced, republican);
%                 hyp = full2one(besthyp, i, parms.ncandidates, parms.nfirm);
                hyp.mean(1) = 0;
               
                hyp.mean(2) = computePrior(pvi, experienced, republican, parms);
                if strcmp(state,'Georgia Special')
                    hyp.mean(2) = 0.25;
                end
                hyp.cov = besthyp.cov;
                hyp.lik = besthyp.lik;
%                 hyp.cov(1)=cov1(i);
%                 hyp.cov(2)=cov2(i);
%                 hyp.lik=lik1(i);
                im{3}.mean{2}{2} = computePrior(pvi, experienced, republican, parms);
%                 im{3}.mean{2}{3} = b_std;
%                 im{3}.mean{1}{2} = a_mu;
%                 im{3}.mean{1}{3} = a_std;
%                 parms.mode = "all";
%                 mask = false(size(unwrap(hyp)));
%                 liksize = size(hyp.lik, 1);
%                 covsize = size(hyp.cov,1);
%                 mask(liksize+covsize+1) = 1;
%                 mask(liksize+covsize+2) = 1;
%                 hypab.mean = hyp.mean(1:2);
%                 hypab = feval(mfun, hypab, @gp_mask, p, hyp, im,...
%                         meanfunc, covfunc,...
%                         likfunc, xs{i}, ys{i}, mask, parms, i, "all");
%                 hyp.mean(1:2) = hypab.mean;
            end
            
%             im{3}.mean{2}{2} = hyp.mean(2);
%             im{3}.mean{1}{2} = hyp.mean(1);
            
            xstar = [linspace(xs{i}(1,1),0,nz).',zeros(1,nz)',ones(1,nz)',...
                            parms.nfirm*ones(1,nz)',republican*zeros(1,nz)'];
            [~,~, fmu, fs2] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, xs{i}, ys{i}, xstar);
            
            predPoll = fmu(end);
            fts(i) = predPoll;
            s2s(i) = fs2(end);
%             parms.prior = [computePrior(pvi, experienced, republican, parms), 0.1];
%             fig = plot_posterior(fmu, fs2, xs{i}(:,1), ys{i}, xstar(:,1), trueVote/100, parms);
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
%             filename = fullfile(stateFolder, plot_title + num2str(parms.j)+ ".jpg");
%             saveas(fig, filename);
%             close;
%             disp(plot_title + " predicted winning rate: " + predPoll);
%             disp(plot_title + " actual votes won: " + trueVote + newline);
            if ~isfield(allRaces, fn)
                allRaces.(fn) = [predPoll, trueVote];
            else
                allRaces.(fn) = [allRaces.(fn), predPoll, trueVote];
            end
        end
    end
end