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
    jobname = "CNNCutoff2018Thres" + pollthres + "Iter" + iter +  "Seed" + seed;
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

%     % iterate cycle/state race
%     allRaces = struct;
%     nz = 200;
%     for i = 1:numel(xs)
%         republican = xs{i}(1,5);
%         xstar = [linspace(xs{i}(1,1)-10,0,nz).',zeros(1,nz)',ones(1,nz)',...
%             parms.nfirm*ones(1,nz)',republican*ones(1,nz)', i*republican*ones(1,nz)'];
% 
%         hyp = full2one(besthyp, i, counter, parms.nfirm);
%         [~, ~, fmu, fs2] = gp(hyp, im, meanfunc, covfunc, likfunc, xs{i}, ys{i}, xstar);
%         
%         fig = plot_posterior(fmu, fs2, xs{i}(:,1), ys{i}, xstar(:,1), i);
%         predPoll = fmu(end);
%         year = raceinfos{i}{1};
%         state = raceinfos{i}{2}{1};
%         candidateName = raceinfos{i}{3};
%         trueVote = raceinfos{i}{4};
%         plot_title = year + " " + state + " " + candidateName;
%         title(plot_title);
%         yearFolder = fullfile(plot_path, num2str(year));
%         stateFolder = fullfile(yearFolder, state);
%         if ~exist(plot_path, 'dir')
%             mkdir(plot_path);
%         end
%         if ~exist(yearFolder, 'dir')
%             mkdir(yearFolder);
%         end
%         if ~exist(stateFolder, 'dir')
%             mkdir(stateFolder);
%         end
%         filename = fullfile(stateFolder, plot_title + ".jpg");
%         saveas(fig, filename);
%         close;
%         disp(plot_title + " predicted winning rate: " + predPoll);
%         disp(plot_title + " actual votes won: " + trueVote + newline);
% 
%         fn = char(state+""+year);
%         fn = fn(~isspace(fn));
%         if ~isfield(allRaces, fn)
%             allRaces.(fn) = [predPoll, trueVote];
%         else
%             allRaces.(fn) = [allRaces.(fn), predPoll, trueVote];
%         end
%     end
    
    data2018 = readData("data/CNNData2018.csv");
    data2018 = indexPollster(data2018, pollthres, "data/CNNData2018idx.csv");
    [xs2018, ys2018, raceinfos2018] = buildTrainCellArrays(data2018, [2018], states);
    p.method = 'LBFGS';
    p.mem = 100;
    p.verbosity = 0;
    p.length = -100;
    liksize = size(besthyp.lik, 1);
    covsize = size(besthyp.cov,1);
    
    for i=1:numel(xs2018)
        xs2018{i}(:,4) = parms.nfirm;
        republican = xs2018{i}(1,5);
        ndays = size(xs2018{i},1);
        ft = zeros(1,ndays);
        s2t = zeros(1,ndays);
        for j = 1:ndays
            xt = xs2018{i}(1:j,:);
            yt = ys2018{i}(1:j);
            par = {meanfunc,covfunc,likfunc, xt, yt};
            hyp = struct; hyp.mean(1:2) = [0; 0.5];
            hyp_race = full2one(besthyp, 1, counter, parms.nfirm); hyp_race.mean(1:2) = [0; 0.5];
            mask = false(size(unwrap(hyp_race)));
            mask(liksize+covsize+1) = 1;
            mask(liksize+covsize+2) = 1;
            hyp = minimize_v2(hyp, @gp_mask, p, hyp_race, inffunc, par{:}, mask);
            hyp_race.mean(1:2) = hyp.mean;
            [~, ~, mu, s2] = gp(hyp_race, inffunc, par{:}, [0,0,1,parms.nfirm,republican]);
            ft(j) = mu; s2t(j) = s2;
        end
        
        ts = xs2018{i}(:,1);
        fig = figure('visible', 'off');
        f = [ft+2*sqrt(s2t); flip(ft-2*sqrt(s2t),1)];
        fill([ts.'; flip(ts.',1)], f, [7 7 7]/8);
        hold on; plot(ts, ft);
        
        year = raceinfos2018{i}{1};
        state = raceinfos2018{i}{2}{1};
        candidateName = raceinfos2018{i}{3};
        trueVote = raceinfos2018{i}{4}/100;
        plot([0], [trueVote], 'rx');
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
    end
    
    posttrain(CNNdata, allRaces, besthyp);
    save(jobname + ".mat");
end