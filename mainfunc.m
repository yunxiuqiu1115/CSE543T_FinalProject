function [varout] = mainfunc(TYPE, mode, tau, plot)

%  main function of obtaining gp posteriors for each election race
%  input:
%    - TYPE: 'GP' or 'LM', specifying the prior model
%    - mode: 
%      -: 1 for doing LOYO on year 1992 to 2016
%      -: 2 for testing 2018 races
%      -: 3 for forecasting 2020 races
%    - tau: forecasting horizon
%    - plot: whether to plot gp
%

    % define output
    varout = "";
    
    % load gmpl, utilities and data
    addpath("gpml-matlab-v3.6-2015-07-07");
    addpath("utilities");
    addpath("data");
    startup;

    % define horizons
    % 8/6/4/3/2/1/0 weeks
    taus = [0,7, 14,21, 28,42, 56];

    % define search space
    search_size = 100;
    if strcmp(TYPE, "GP")==1
        p = sobolset(3, 'Skip', 1);
    else 
        search_size = 20;
        p = linspace(0,0.1,search_size);
    end
    
    if mode>=2
        % no search needed
        % forecasting year 2018 or 2020
            for j=1:search_size
                if strcmp(TYPE, "GP")==1
                    ls = p(j,1)*(56-7)+7; % 7-56
                    os = p(j,2)/20; % 0%-5%
                    lik = p(j,3)/20; % 0%-5%
                    myrun(tau,TYPE, ls, os, lik, j, mode, plot);
                else 
                    % linear model does not have ls/os
                    ls = 0; os = 0;
                    lik = p(j); % 0%-10%
                    myrun(tau,TYPE, ls, os, lik, j, mode, plot);
                end
            end
%         end
        else
        % performing a sobel sequence search
%         for i=1:numel(taus)
            for j=1:search_size
                if strcmp(TYPE, "GP")==1
                    ls = p(j,1)*(56-7)+7; % 7-56
                    os = p(j,2)/20; % 0%-5%
                    lik = p(j,3)/20; % 0%-5%
                    myrun(tau,TYPE, ls, os, lik, j, mode, plot);
                else 
                    % linear model does not have ls/os
                    ls = 0; os = 0;
                    lik = p(j); % 0%-10%
                    myrun(tau,TYPE, ls, os, lik, j, mode, plot);
                end
            end
%         end
    end
end

function myrun(tau,type, ls, os, lik, j, mode, plot)
%  load 1992-2016 data
%  feature of data includes: 
%       - cycle: election year
%       - state: state that holds election in the election year
%       - Candidateidentier: include election year, state and name of candidate
%       - Percentage_of_Vote_won_x: actual vote share, 0 for future elections
%       - samplesize: number of respondences in the poll
%       - daysLeft: days prior to election when the poll was conducted
%       - numberSupport: number of respondences that supports the candidate
%       - pollster: name of the polling company/institution
%       - Republican: 1 if the candidate is republican, 0 otherwise
%       - Democrat: 1 if the candidate is democratic, 0 otherwise
%       - pvi: Cook partisan voting index of the election state
%       - experienced: 1 if the candidate has served in any political office, 0 otherwise

    CNNdata = readData("data/CNNdata1992-2016.csv");
    CNNdata2018 = readData("data/CNNData2018.csv");
    CNNdata2018(:, ["candidate_name"]) = [];
    
    if mode==2
        % test on 2018 data
        % load 2018 data
        CNNdata2018 = readData("data/CNNData2018.csv");
        CNNdata2018(:, ["candidate_name"]) = [];
        CNNdata = vertcat(CNNdata, CNNdata2018);
    elseif mode==3
        % forecast 2020 races
        CNNdata2020 = readData("data/CNNData2020.csv");
        CNNdata2020(:, ["candidate_name"]) = [];
        CNNdata2020.Percentage_of_Vote_won_x = zeros(size(CNNdata2020,1),1);
        CNNdata = vertcat(CNNdata, CNNdata2020);
    end

    % split data into cell arrays
    % obtain features (days before election) and regressor (polling proportions)
    years = unique(CNNdata.cycle);
    states = unique(CNNdata.state);
    [xs, ys, raceinfos] = buildTrainCellArrays(CNNdata, years, states);
    years2 = unique(CNNdata2018.cycle);
    states2 = unique(CNNdata2018.state);
    [xp, yp, raceinfos2] = buildTrainCellArrays(CNNdata2018, years2, states2);

    % shrink data according to the forecasting horizons
    counter = size(xs,1);
    for i=1:counter
        idx = xs{i}(:,1) <= -tau;
        xs{i} = xs{i}(idx,:);
        ys{i} = ys{i}(idx);
    end
    
    counter = size(xp,1);
    for i=1:counter
        idx = xp{i}(:,1) <= -tau;
        xp{i} = xp{i}(idx,:);
        yp{i} = yp{i}(idx);
    end

    % define hyp (hyperparameter) struct
    %   - hyp.mean:
    %       - prior mean of the linear trend slope
    %       - prior mean of the linear trend intercept
    %   - hyp.cov:
    %       - log of length scale in Matérn kernel
    %       - log of outpur scale in Matérn kernel
    %       - 1 / log of prior std of the linear trend slope
    %       - log of prior std of the linear trend intercept
    %   - hyp.lik: log of observation noise std
    if mode == 1
        meanfunc = [];
        covfunc = @covSEiso;
        likfunc = @likGauss;
        hyp = struct('mean', [], 'cov', [0 0], 'lik', -1);
        xp = [0];
        hyp = minimize(hyp, @gpsum, -100, @infExact, meanfunc, covfunc, likfunc, xs, ys);
        [mu s2] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, xs{1}(:,1), ys{1}, xp);
%         [mu s2] = gphelper(hyp, @infGaussLik, meanfunc, covfunc, likfunc, xs, ys, xp);
        f = [mu+2*sqrt(s2); flipdim(mu-2*sqrt(s2),1)];
        fill([xp; flipdim(xp,1)], f, [7 7 7]/8)
        hold on; plot(xp, mu);
    else
        hyp.cov(1) = log(ls);
        hyp.cov(2) = log(os);

        % model() defines prior distribution of linear trends
        [~,~,~,~, prior] = model();
        sigma_ml = prior.slope(2);
        sigma_mc = prior.intercept(2);
        hyp.cov(3) = log(1/sigma_ml);
        hyp.cov(4) = log(sigma_mc);
        hyp.lik = log(lik);
    end

    % parms struct specifies keyword arguments
    %   - tau: forecasting horizon
    %   - j: current index of searching sequence
    %   - type: prior model, 'GP' or 'LM'
    %   - plot:
    %       - 1: generate plots of underlying voter preference for each election race
    %       - 0: just obtain posterior belief of voter preference on election day
    %   - test_year: validation year in LOYO process or test year in testing process
    %   - coefs: the coefs of prior linear model for the linear trend intercept
    parms.tau = tau;
    parms.j = j;
    parms.type = type;
    parms.plot = plot;
    
    % plot days bin
    parms.BIN = 30;
    parms.test_year = 2018;
    parms.coefs = priorModel(CNNdata2018, parms.test_year);
    plot_path = "plots/" + type + "MargLinTre" + num2str(parms.test_year)+"_"+num2str(tau);
%     [allRaces, fts, s2s] = gpm(hyp, xs, ys, raceinfos, plot_path, parms);
    % Do the forecasting
end
