taus = [0,14,28,42,90,120];

for i=1:numel(taus)
    % read race data
    CNNdata = readData("data/CNNData.csv");
%     CNNdata = readData("data/CNNData1992to2018.csv");
    pollthres = 50;
    seed = 1;
    [CNNdata,pollster2idx] = indexPollster(CNNdata, pollthres);
    LAST_TIME = taus(i); % positive
    jobname = "AllMargLinTre2016" + num2str(LAST_TIME) + "Thres" + pollthres +  "Seed" + seed;
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
        idx = xs{i}(:,1) <= -LAST_TIME;
        xs{i} = xs{i}(idx,:);
        ys{i} = ys{i}(idx);

        vs(i) = raceinfos{i}{4}/100;
    end

    parms.vs = vs;

    % define model

    meanmask = [true, false, false, false, false];
    covmask = [false, true, true, false, false];
    mc = {@meanConst};
    ml = {@meanLinear};
    mml = {@meanMask, meanmask, ml};
    mmc = {@meanMask, meanmask, mc};

    % try no firm
    meanfunc = {@meanSum, {mml, mmc}};
    cm = {@covMaterniso, 3};
    mb = {@logsqrtbinom};
    cdb = {@covDiag, mb};
    cmd = {@covMask, {covmask, cdb}};
    cmm = {@covMask, {meanmask, cm}};

    cs = {@covPoly, 1};
    ci = {@covConst};
    cms = {@covMask, {meanmask, cs}};

    % try no firm
    covfunc = {@covSum, {cmm, cmd, cms, ci}};
    likfunc = @likGauss;
    inffunc = @infExact;

    mu_ml = 0.0; sigma_ml = 0.002;
    mu_mc = 0.5; sigma_mc = 0.1;
    pc_ml = {@priorClamped};
    pc_mc = {@priorClamped};
    prior.mean = {pc_ml, pc_mc};


    mu_ls = log(30); sigma_ls = log(90/30)/2;
    mu_os = log(1/40); sigma_os = log(2);
    pg_ls = {@priorGauss, mu_ls, sigma_ls^2};
    pg_os = {@priorGauss, mu_os, sigma_os^2};

    pg_ml = {@priorClampedMulti};
    pg_mc = {@priorClamped};

    prior.cov = {pg_ls, pg_os, pg_ml, pg_mc};

    mu_lik = log(0.03); sigma_lik = log(2);
    pg_lik = {@priorGauss, mu_lik, sigma_lik^2};
    prior.lik = {pg_lik};


    [xs, ys, raceinfos] = buildTrainCellArrays(CNNdata, years, states);
    counter = size(xs,1);

    allRaces = struct;
    n = numel(xs);
    fts = zeros(n,1);
    s2s = zeros(n,1);
    for i=1:counter
%     for i=803
        if raceinfos{i}{1}>=0
            idx = xs{i}(:,1) <= -LAST_TIME;
            xs{i} = xs{i}(idx,:);
            ys{i} = ys{i}(idx);
        
            year = raceinfos{i}{1};
            state = raceinfos{i}{2}{1};
            candidateName = raceinfos{i}{3};
            trueVote = raceinfos{i}{4};
            pvi =raceinfos{i}{5};
            experienced =raceinfos{i}{6};
            republican = raceinfos{i}{7};
            fn = char(state+""+year);
            fn = fn(~isspace(fn));
            b = computePrior(pvi, experienced, republican);

            im = {@infPrior, inffunc, prior};

            hyp.mean = [0;b];
            hyp.cov = [mu_ls;mu_os;log(0);log(sigma_ml);log(sigma_mc)];
            hyp.lik = [mu_lik];
            hyp.cov = besthyp.cov;
            hyp.lik = besthyp.lik;
            nz = 200;
            if numel(xs{i})==0
                predPoll = b;
                fts(i) = predPoll;
                s2s(i) = exp(2*mu_lik)+sigma_mc^2;
            else        
            xstar = [linspace(xs{i}(1,1),0,nz).',zeros(1,nz)',ones(1,nz)',...
                                parms.nfirm*ones(1,nz)',republican*zeros(1,nz)'];
            [~, ~, fmu, fs2] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, xs{i}, ys{i}, xstar);
            predPoll = fmu(end);
            fts(i) = predPoll;
            s2s(i) = fs2(end);
%             fig = plot_posterior(fmu, fs2, xs{i}(:,1), ys{i}, xstar(:,1), trueVote/100, i,hyp.mean(1),hyp.mean(2));
%             plot_title = year + " " + state + " " + candidateName + " ";
%             title(plot_title);
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
            end
            
            if ~isfield(allRaces, fn)
                allRaces.(fn) = [predPoll, trueVote];
            else
                allRaces.(fn) = [allRaces.(fn), predPoll, trueVote];
            end
        end
    end

    posttrain(raceinfos,fts,s2s,allRaces,hyp, LAST_TIME, method);
end