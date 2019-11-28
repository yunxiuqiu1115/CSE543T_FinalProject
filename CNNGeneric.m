% CNN Forecasting generic model
function CNNGeneric(pollthres,iter,seed)
    if nargin < 1, pollthres = 50; iter=5;seed = 1; end
    if nargin < 2, iter=5;seed = 1; end
    if nargin < 3, seed = 1; end
    addpath("gpml-matlab-v3.6-2015-07-07");
    addpath("utilities");
    addpath("data");
    startup;

    % read race data
    CNNdata = readData("data/CNNData.csv"); 
    CNNdata = indexPollster(CNNdata, pollthres, "data/CNNDataidx.csv");
    jobname = "CNNValid2016Thres" + pollthres + "Iter" + iter +  "Seed" + seed;
    plot_path = "plots/" + jobname;

    parms.mode = true;
    % pollsters having less than threshold of polls will be indexed by nfirm
    parms.nfirm = max(CNNdata.pollsteridx);
    parms.days = min(CNNdata.daysLeft);
    years = unique(CNNdata.cycle);
    states = unique(CNNdata.state);

    % build training cell arrays
    [xs, ys, raceinfos] = buildTrainCellArrays(CNNdata, years, states);
    counter = size(xs,1);
    parms.ncandidates = counter;
    parms.a = ones(counter,1)*0.5;
    
    for i=1:counter
        pvi =raceinfos{i}{5};
        experienced =raceinfos{i}{6};
        republican = xs{i}(1,5);
        parms.a(i) = computePrior(pvi, experienced, republican);
        if raceinfos{i}{1}>=2016
            idx = xs{i}(:,1) <= -90;
            xs{i} = xs{i}(idx,:);
            ys{i} = ys{i}(idx);
        end
    end

    % define model
    [meanfunc, covfunc, likfunc, inffunc, prior] = model(parms);
    im = {@infPrior, inffunc, prior};
    par = {meanfunc,covfunc,likfunc, xs, ys};

    % training
    disp("start training...");
    hyp = sample_separate_prior(prior, parms, counter, seed);
    hyp = fixLearn(hyp, im, par{:}, iter, parms);
    allRaces = forcastAllRaces(hyp, xs, ys, raceinfos, plot_path, parms);
    
%     data2018 = readData("data/CNNData2018.csv");
%     data2018 = indexPollster(data2018, pollthres, "data/CNNData2018idx.csv");
    [validxs, validys, validraceinfos] = buildTrainCellArrays(CNNdata, (2016), states);
    parms.valididx = size(xs,1) - size(validxs,1);
    [fts, s2s] = performForcasting(hyp, validxs, validys, validraceinfos, plot_path, parms);
    
    posttrain(CNNdata, allRaces, hyp);
    save(jobname + ".mat");
end