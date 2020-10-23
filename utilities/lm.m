function [allRaces,fts,s2s] = lm(hyperparameter, xs, ys, raceinfos, plot_path, parms)

%  obtain posterior belief of voter preferences on election day using lm model
%  input:
%    - hyperparameter: hyperparameter of mean/cov/lik
%    - xs: cell array of polling data in form of [tau, polling proportion, sample size]
%    - ys: cell array of [polling proportion]
%    - raceinfos: election metadata [year, state, candidatename, actual vote share, pvi, experienced, isRep]
%    - plot_path: path of plots if parms.plot==1
%    - parms: function keyword arguments
%       - tau: forecasting horizon
%       - j: current index of searching sequence
%       - type: prior model, 'GP' or 'LM'
%       - plot:
%           - 1: generate plots of underlying voter preference for each election race
%           - 0: just obtain posterior belief of voter preference on election day
%       - test_year: validation year in LOYO process or test year in testing process
%       - coefs: the coefs of prior linear model for the linear trend intercept
%
% Caller of this function should specify the keyword arguments parms.
%

    % define model
    [meanfunc, covfunc, likfunc, inffunc, prior] = model();
    mu_ml = prior.slope(1);
    sigma_mc = prior.intercept(2);
    
    % AllRaces is a struct with key year+state.
    % Value of allRaces struct is an array of [model posterior mean, actual
    % vote, ...] of all candidates for the year/state race
    allRaces = struct;
    n = numel(xs);
    
    % fts/s2s are the array of posterior mean/var on election day of all races
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
        party = raceinfos{i}{7};
        
        % key for allRaces struct
        fn = char(state+""+year);
        fn = fn(~isspace(fn));
        
        % precompute prior mean of linear trend intercept
        mu_b = computePrior(pvi, experienced, party, parms);
        
        % obtain lm posterior
        if numel(xs{i})==0
            % if no polls are avaiable
            % use prior
            fts(i) = mu_b;
            s2s(i) = sigma_mc^2;
        else 
            % if there is data avaiable   
            x = xs{i};
            y = ys{i};
            if parms.plot==0
                % obtain gp posterior on election day only
                % test position is [day 0, 0 polling porportion, 1 samplesize]
                xstar = [0,0,1];
                
                % redefine covfunc since it will be overwrited
                % lm model does not have Matérn covariance.
                [~, covfunc, ~, ~, ~] = model();
                covfunc{2} = covfunc{2}(2:end);
                
                % define hyp for one specific race
                % hyp is exactly the same of hyperparameter except potentially
                % different prior mean of linear trend intercept and an
                % absense of Matérn covariance
                hyp.mean = [mu_b];
                hyp.cov = hyperparameter.cov(3:4);
                hyp.lik = hyperparameter.lik;

                [~, ~, fts(i), s2s(i)] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, x, y, xstar);
            else
                % plot upon required
                % define 200 test location from the earliest poll day to
                % election day
                nz = 200;
                xstar = [linspace(xs{i}(1,1),0,nz).',zeros(1,nz)',ones(1,nz)'];
                [~, ~, fmu, fs2] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, xs{i}, ys{i}, xstar);
                fts(i) = fmu(end);
                s2s(i) = s2s(end);
                parms.prior = [mu_b, simga_mc];
                fig = plot_posterior(fmu, fs2, xs{i}(:,1), ys{i}, xstar(:,1), parms);
                plot_title = year + " " + state + " " + candidateName;
%                 title(plot_title);
                
                % save plot to files
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
      
        % update allRaces struct
        if ~isfield(allRaces, fn)
            allRaces.(fn) = [fts(i), trueVote];
        else
            allRaces.(fn) = [allRaces.(fn), fts(i), trueVote];
        end
    end
end