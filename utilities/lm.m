function [allRaces,fts,s2s] = lm(besthyp, xs, ys, raceinfos, plot_path, parms)

%  obtain posterior belief of voter preferences on day 0 using lm model
%  input:
%    - besthyp: hyperparameters
%    - x: [tau, polling proportion, sample size]
%    - y: [polling proportion]
%    - raceinfos: [year, state, candidatename, actual vote share, pvi, experienced, isRep]
%    - plot_path: path of plots if parms.plot==1
%    - parms: parameters related to functionaity of this function.
%
% WARNING: This function is not for the purpose of being directly called. 
% Function calling this function is responsible for specifying parms for
% required functionalities.

    % define model
    [meanfunc, covfunc, likfunc, inffunc, prior] = model(parms);
    mu_ml = prior.slope(1);
    sigma_mc = prior.intercept(2);
    allRaces = struct;
    n = numel(xs);
    fts = zeros(n,1);
    s2s = zeros(n,1);
    
    % iterate every race
    for i = 1:n
        % obtain meta data
        year = raceinfos{i}{1};
        state = raceinfos{i}{2}{1};
        candidateName = raceinfos{i}{3};
        trueVote = raceinfos{i}{4};
        pvi =raceinfos{i}{5};
        experienced =raceinfos{i}{6};
        republican = raceinfos{i}{7};
        fn = char(state+""+year);
        fn = fn(~isspace(fn));
        mu_b = computePrior(pvi, experienced, republican, parms);
        if numel(xs{i})==0
            % if no polls are avaiable
            % use prior
            fts(i) = mu_b;
            s2s(i) = sigma_mc^2;
        else 
            x = xs{i};
            y = ys{i};
            % obtain lm posterior
            if parms.plot==0
                xstar = [0,0,1];
                covfunc{2} = covfunc{2}(2:end);
                hyp.mean(1) = mu_ml;
                hyp.mean(2) = mu_b;
                hyp.cov = besthyp.cov(3:4);
                hyp.lik = besthyp.lik;

                [~, ~, fmu, fs2] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, x, y, xstar);

                predPoll = fmu(end);
                fts(i) = predPoll;
                s2s(i) = fs2(end);
            % plot upon requests
            else
                xstar = [linspace(xs{i}(1,1),0,nz).',zeros(1,nz)',ones(1,nz)'];
                [~, ~, fmu, fs2] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, xs{i}, ys{i}, xstar);
                parms.prior = [mu_b, simga_mc];
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
        end
      
        if ~isfield(allRaces, fn)
            allRaces.(fn) = [fts(i), trueVote];
        else
            allRaces.(fn) = [allRaces.(fn), fts(i), trueVote];
        end
    end
end