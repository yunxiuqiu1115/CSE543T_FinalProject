function [allRaces,fts,s2s] = forcastAllRaces(besthyp, xs, ys, raceinfos, plot_path, parms)
    % iterate cycle/state race
    [meanfunc, covfunc, likfunc, inffunc, prior] = model(parms);
    mu_ml = prior.slope(1);
    sigma_mc = prior.intercept(2);
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
        pvi = raceinfos{i}{5};
        experienced = raceinfos{i}{6};
        republican = raceinfos{i}{7};
        fn = char(state+""+year);
        fn = fn(~isspace(fn));
        mu_b = computePrior(pvi, experienced, republican, parms);
        if numel(xs{i})==0
            % if there is no data avaiable
            % use prior
            if ~isfield(allRaces, fn)
                allRaces.(fn) = [mu_b, trueVote];
            else
                allRaces.(fn) = [allRaces.(fn), mu_b , trueVote];
            end
            fts(i) = mu_b;
            s2s(i) = sigma_mc^2 + exp(besthyp.cov(2));
        else
            % if there is data avaiable     
            hyp.mean(1) = mu_ml;
            hyp.mean(2) = mu_b;
            if strcmp(state,'Georgia Special')
                hyp.mean(2) = 0.25;
            end
            hyp.cov = besthyp.cov;
            hyp.lik = besthyp.lik;
            xstar = [0,0,1];
            [~, ~, fmu, fs2] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, xs{i}, ys{i}, xstar);
            predPoll = fmu(end);
            fts(i) = predPoll;
            s2s(i) = fs2(end);
            
            if parms.plot==1
                xstar = [linspace(xs{i}(1,1),0,nz).',zeros(1,nz)',ones(1,nz)'];
                [~, ~, fmu, fs2] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, xs{i}, ys{i}, xstar);
                parms.prior = [mu_b, sigma_mc];
                fig = plot_posterior(fmu, fs2, xs{i}(:,1), ys{i}, xstar(:,1), trueVote/100, parms);
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
                filename = fullfile(stateFolder, plot_title + num2str(parms.j)+ ".jpg");
                saveas(fig, filename);
                close;
            end   
             
            if ~isfield(allRaces, fn)
                allRaces.(fn) = [predPoll, trueVote];
            else
                allRaces.(fn) = [allRaces.(fn), predPoll, trueVote];
            end
        end
    end
end