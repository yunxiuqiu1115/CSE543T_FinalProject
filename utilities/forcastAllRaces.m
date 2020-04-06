function [allRaces,fts,s2s] = forcastAllRaces(besthyp, xs, ys, raceinfos, plot_path, parms)
    % iterate cycle/state race
    [meanfunc, covfunc, likfunc, inffunc, prior] = model(parms);
    im = {@infPrior, inffunc, prior};
    allRaces = struct;
    n = numel(xs);
    fts = zeros(n,1);
    s2s = zeros(n,1);
    nz = 200;
    for i = 1:n
        year = raceinfos{i}{1};
        state = raceinfos{i}{2}{1};
        candidateName = raceinfos{i}{3};
        trueVote = raceinfos{i}{4};
        pvi =raceinfos{i}{5};
        experienced =raceinfos{i}{6};
        republican = (raceinfos{i}{7}+1)/2;
        fn = char(state+""+year);
        fn = fn(~isspace(fn));
        if numel(xs{i})==0
            if i<=parms.ncandidates
                predPoll = besthyp.mean(param.ncandidates+i) + besthyp.mean(end)*republican;
            else 
                predPoll = sample_prior(prior).mean(1) + besthyp.mean(end)*republican;
            end
            if ~isfield(allRaces, fn)
                allRaces.(fn) = [predPoll, trueVote];
            else
                allRaces.(fn) = [allRaces.(fn), predPoll, trueVote];
            end
            fts(i) = predPoll;
            s2s(i) = (exp(hyp.lik) + exp(besthyp.cov(end)))^2;
        else
            % republican = xs{i}(1,5);
            xstar = [linspace(xs{i}(1,1)-10,0,nz).',zeros(1,nz)',ones(1,nz)',...
                parms.nfirm*ones(1,nz)',republican*ones(1,nz)'];
            if i<=parms.ncandidates
                im{3}.mean{2}{2} = parms.a(i);
                hyp = full2one(besthyp, i, parms.ncandidates, parms.nfirm);
            else
                im{3}.mean{2}{2} = computePrior(pvi, experienced, republican);
                hyp = full2one(besthyp, 1, parms.ncandidates, parms.nfirm);
                hyp.mean(1:2) = sample_prior(prior).mean(1:2);
            end
            
            [~, ~, fmu, fs2] = gp(hyp, im, meanfunc, covfunc, likfunc, xs{i}, ys{i}, xstar);

            fig = plot_posterior(fmu, fs2, xs{i}(:,1), ys{i}, xstar(:,1), i);
            predPoll = fmu(end);
            fts(i) = predPoll;
            s2s(i) = fs2(end);
            plot_title = year + " " + state + " " + candidateName;
            title(plot_title);
            yearFolder = fullfile(plot_path, num2str(year));
            stateFolder = fullfile(yearFolder, state);
            if ~exist(plot_path, 'dir')
                mkdir(plot_path);
            end
            if ~exist(yearFolder, 'dir')
                mkdir(yearFolder);
            end
            if ~exist(stateFolder, 'dir')
                mkdir(stateFolder);
            end
            filename = fullfile(stateFolder, plot_title + ".jpg");
            saveas(fig, filename);
            close;
            disp(plot_title + " predicted winning rate: " + predPoll);
            disp(plot_title + " actual votes won: " + trueVote + newline);
            if ~isfield(allRaces, fn)
                allRaces.(fn) = [predPoll, trueVote];
            else
                allRaces.(fn) = [allRaces.(fn), predPoll, trueVote];
            end
        end
    end
end