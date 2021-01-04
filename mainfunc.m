function [varout] = mainfunc(tau)
    global xList nlZsumList iter;
    xList=[];
    nlZsumList=[];
    iter=0;
    warning off;
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
    CNNdata2020 = readData("data/CNNdata2020.csv");
    CNNdata2020(:, ["candidate_name"]) = [];
    % TX, Mich, NC
    %CNNdata2020.Percentage_of_Vote_won_x = zeros(size(CNNdata2020,1),1);
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
    
%     meanfunc = [];
%     covfunc = {@covSEiso};
     mu = 1.0; s2 = 0.01^2;
     prior.mean = {{@priorGauss, mu, s2}; {'priorLaplace', mu, s2}};
     prior.lik = {{@priorDelta}};
    inffunc = {@infPrior, @infGaussLik, prior};
%     hyp = struct('mean', [], 'cov', [0 0], 'lik', -1);
    meanfunc = [];
    meanmask = [true, false, false];
    % mask of polling porportion and sample size
    
    
    covmask = [false, true, true];
    cm = {@covMaterniso, 1};
    mb = {@logsqrtbinom};
    cdb = {@covDiag, mb};
    cmd = {@covMask, {covmask, cdb}};
    cmm = {@covMask, {meanmask, cm}};
    cs = {@covLINiso};
    ci = {@covConst};
    cms = {@covMask, {meanmask, cs}};     
    covfunc = {@covSum, {cmm, cms}};
    sf = 1; ell = 0.7;
    likfunc = @likGauss;
    sn = 0.2;
    hyp = struct('mean', [], 'cov', [0 0 0], 'lik', -1);
    hyp = minimize(hyp, @gpmean, -100, inffunc, meanfunc, covfunc, likfunc, xs(1:200,:), ys(1:200,:));
    disp(hyp);
    figure
    plot(xList,nlZsumList);
    xlabel('iteration');
    ylabel('sum of negative log marginal likelihood');
    saveas(gcf,'nlZSum.png');
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
    posttrain(raceinfos,fts,s2s,allRaces,hyp, tau, parms)
end