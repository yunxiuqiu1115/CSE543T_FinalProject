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
    CNNdata = readData("data/CNNdata.csv");
    CNNdata = indexPollster(CNNdata, pollthres, "data/CNNdataidx.csv");
    jobname = "CNN2018Thres" + pollthres + "Iter" + iter +  "Seed" + seed;
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

    % define model
    [meanfunc, covfunc, likfunc, inffunc, prior] = model(parms);
    im = {@infPrior, inffunc, prior};
    par = {meanfunc,covfunc,likfunc, xs, ys};

    % training
    disp("start training...");
    hyp = sample_separate_prior(prior, parms, counter, seed);
    besthyp = fixLearn(hyp, im, par{:}, iter);

    % iterate cycle/state race
    allRaces = struct;
    nz = 200;
    for i = 1:numel(xs)
        republican = xs{i}(1,5);
        xstar = [linspace(xs{i}(1,1)-10,0,nz).',zeros(1,nz)',ones(1,nz)',...
            parms.nfirm*ones(1,nz)',republican*ones(1,nz)', i*republican*ones(1,nz)'];

        hyp = full2one(besthyp, i, counter, parms.nfirm);
        [~, ~, fmu, fs2] = gp(hyp, im, meanfunc, covfunc, likfunc, xs{i}, ys{i}, xstar);
        
%         ndays = size(xs{1},1);
%         ft = zeros(1,ndays);
%         s2t = zeros(1,ndays);
%         nlZt = zeros(1,ndays);
%         for j = 1:ndays
%             xt = x(1:j,:);
%             yt = y(1:j);
%             par = {meanfunc,covfunc,likfunc, xt, yt};
%             hyp = feval(mfun, hyp,  @gp, -100, im, par{:}); 
%             [~, ~, mu, s2] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, xt, yt, [0,0,1]);
%             [nlZ, ~] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, x, y);
%             ft(i) = mu; s2t(i) = s2; nlZt(i) = nlZ;
%         end
        
        fig = plot_posterior(fmu, fs2, xs{i}(:,1), ys{i}, xstar(:,1), i);
        predPoll = fmu(end);
        year = raceinfos{i}{1};
        state = raceinfos{i}{2}{1};
        candidateName = raceinfos{i}{3};
        trueVote = raceinfos{i}{4};
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

        fn = char(state+""+year);
        fn = fn(~isspace(fn));
        if ~isfield(allRaces, fn)
            allRaces.(fn) = [predPoll, trueVote];
        else
            allRaces.(fn) = [allRaces.(fn), predPoll, trueVote];
        end
    end
    
    posttrain(CNNdata, allRaces, besthyp);
    save(jobname + ".mat");
end