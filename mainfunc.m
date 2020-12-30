function [varout] = mainfunc(tau)
    % import packages
    addpath("gpml-matlab-v3.6-2015-07-07");
    addpath("utilities");
    addpath("data");
    startup;
    
    % import data
    CNNdata = readData("data/CNNdata1992-2016.csv");
    CNNdata2018 = readData("data/CNNData2018.csv");
    CNNdata2018(:, ["candidate_name"]) = [];
    CNNdata = vertcat(CNNdata, CNNdata2018);
    CNNdata2020 = readData("data/CNNData2020.csv");
    CNNdata2020(:, ["candidate_name"]) = [];
    CNNdata2020.Percentage_of_Vote_won_x = zeros(size(CNNdata2020,1),1);
    CNNdata = vertcat(CNNdata, CNNdata2020);
    
    % preprocess the data
    years = unique(CNNdata.cycle);
    states = unique(CNNdata.state);
    [xs, ys, raceinfos] = buildTrainCellArrays(CNNdata, years, states);
    
    % shrink data according to the forecasting horizons
    counter = size(xs,1);
    for i=1:counter
        idx = xs{i}(:,1) <= -tau;
        xs{i} = xs{i}(idx,:);
        ys{i} = ys{i}(idx);
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
    meanfunc = [];
    meanmask = [true, false, false];
    % mask of polling porportion and sample size
    covmask = [false, true, true];
    cm = {@covMaterniso, 3};
    mb = {@logsqrtbinom};
    cdb = {@covDiag, mb};
    cmd = {@covMask, {covmask, cdb}};
    cmm = {@covMask, {meanmask, cm}};
    cs = {@covLINiso};
    ci = {@covConst};
    cms = {@covMask, {meanmask, cs}};     
    covfunc = {@covSum, {cmm, cmd, cms, ci}};
    sf = 1; ell = 0.7;
    likfunc = @likGauss;
    sn = 0.2;
    hyp = struct('mean', [], 'cov', [3.44998754583159 -3.68887945411394 6.21460809842219 -2.30258509299405], 'lik', -3.68887945411394);
    hyp = minimize(hyp, @gpsum, -100, @infExact, meanfunc, covfunc, likfunc, xs(1:873,:), ys(1:873,:));
    disp(hyp);
    % plotting parameters
    parms.tau = tau;
    parms.j = 1;
    parms.type = "GP";
    parms.plot = 1;
    
    % plot days bin
    parms.BIN = 30;
    % forecasting 2018 races
    parms.test_year = 2020;
    % precompute coefs of prior linear model of the linear trend intercept
    parms.coefs = priorModel(CNNdata, parms.test_year);
    plot_path = "plots/" + parms.type + "MargLinTre"+num2str(parms.test_year)+"_"+num2str(tau);
    [allRaces, fts, s2s] = gpm(hyp, xs, ys, raceinfos, plot_path, parms);
end