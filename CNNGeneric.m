% CNN Forecasting generic model
function CNNGeneric(pollthres,iter,seed)
    if nargin < 1, pollthres = 50; iter=1;seed = 1; end
    if nargin < 2, iter=5;seed = 1; end
    if nargin < 3, seed = 1; end
    addpath("gpml-matlab-v3.6-2015-07-07");
    addpath("utilities");
    addpath("data");
    startup;

    % read race data
    CNNdata = readData("data/CNNData.csv");
    [CNNdata,pollster2idx] = indexPollster(CNNdata, pollthres);
    LAST_TIME = 0; % positive 20160Thres50Seed1
    jobname = "AllMargLinTre2016" + num2str(LAST_TIME) + "Thres" + pollthres + "Iter" + iter +  "Seed" + seed;
    disp(jobname);

    plot_path = "plots/" + jobname;
    parms.mode = "all";
    
    method = parms.mode;
    % pollsters having less than threshold of polls will be indexed by nfirm
    parms.nfirm = max(CNNdata.pollsteridx);
    parms.days = min(CNNdata.daysLeft);
    years = unique(CNNdata.cycle);
    states = unique(CNNdata.state);
    
    % try no firm
    parms.nfirm = 0;

    % build training cell arrays
    [xs, ys, raceinfos] = buildTrainCellArrays(CNNdata, years(1:end-1), states);
    counter = size(xs,1);
    parms.ncandidates = counter;
    parms.a = ones(counter,1)*0.5;
    
    vs = zeros(counter,1);
    for i=1:counter
        pvi =raceinfos{i}{5};
        experienced =raceinfos{i}{6};
        republican = xs{i}(1,5);
        parms.a(i) = computePrior(pvi, experienced, republican);
        vs(i) = raceinfos{i}{4}/100;
    end
    
    parms.vs = vs;

    % define model
    [meanfunc, covfunc, likfunc, inffunc, prior] = model(parms);
    im = {@infPrior, inffunc, prior};
    par = {meanfunc,covfunc,likfunc, xs, ys};
    

    % training
    disp("start training...");
    cov1 = zeros(counter,1);
    cov2 = zeros(counter,1);
    lik1 = zeros(counter,1);
    for i=1:counter
       hyp = sample_separate_prior(prior, parms, counter, 1);
       parms.i = i;
       [hyp, nlZ] = fixLearn(hyp, im, par{:}, iter, parms);
       cov1(i) = hyp.cov(1);
       cov2(i) = hyp.cov(2);
       lik1(i) = hyp.lik;
    end
    [f,xi] = ksdensity(cov1); 
    figure;
    plot(xi,f);
    legend('log(ls)');
    hold on;
    [f,xi] = ksdensity(cov2); 
    figure;
    plot(xi,f);
    legend('log(os)');
    hold on;
    [f,xi] = ksdensity(lik1); 
    figure;
    plot(xi,f);
    legend('log(lik)');
    
    
%     best_nlZ = 0;
%     best_hyp = 0;
%     for i=1:1
%         hyp = sample_separate_prior(prior, parms, counter, i);
%         [hyp, nlZs] = fixLearn(hyp, im, par{:}, iter, parms);
%         nlZ = nlZs(end);
%         disp("hyp");
%         disp(hyp);
%         disp("nlzs");
%         disp(nlZs);
%         if nlZ<best_nlZ
%            best_nlZ = nlZ;
%            best_hyp = hyp;
%         end
%     end
    
%     hyp = best_hyp;

    hyperparameters = sample_separate_prior(prior, parms, counter, 1);
    num_samples = size(xs,1); 
    hyp.mean = hyp.mean(2*num_samples+1:end);
    hyp.cov = hyperparameters.cov(1:2);
    hyp.lik = hyperparameters.lik;
    means = reshape(hyperparameters.mean, num_samples,2);
    
    logpdf = @(hyp)hypposterior(hyp, hyperparameters, im, par{:}, parms);
    smp = hmcSampler(logpdf,unwrap(hyp));
    [chain, endpoint, accratio] = drawSamples(smp, 'NumSamples', 100);
    
    hyperparameters.cov(1:2) = hyp.cov(1:2);
    hyperparameters.lik = hyp.lik;
    
    [xs, ys, raceinfos] = buildTrainCellArrays(CNNdata, years, states);
    counter = size(xs,1);

    for i=1:counter
        if raceinfos{i}{1}>=2016
            idx = xs{i}(:,1) <= -LAST_TIME;
            xs{i} = xs{i}(idx,:);
            ys{i} = ys{i}(idx);
        end
    end
    
    [allRaces, fts, s2s] = forcastAllRaces(hyp, xs, ys, raceinfos, plot_path, parms);
     posttrain(raceinfos,fts,s2s,allRaces,hyp, LAST_TIME, method);
     
     save(jobname + ".mat");

end