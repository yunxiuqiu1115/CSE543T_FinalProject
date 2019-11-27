function [fts,s2s] = performForcasting(besthyp, xs, ys, raceinfos, plot_path, parms)
    [meanfunc, covfunc, likfunc, inffunc, prior] = model(parms);
    im = {@infPrior, inffunc, prior};
    DAYS = [-90,-60,-30,-14,-7,-1];
    
    fts = cell(numel(xs),1);
    s2s = cell(numel(xs),1);
    for i=1:numel(xs)
        xs{i}(:,4) = parms.nfirm;
        republican = xs{i}(1,5);
        fts{i} = zeros(1,numel(DAYS));
        s2s{i} = zeros(1,numel(DAYS));
        for j = 1:numel(DAYS)
            idx = xs{i}(:,1) < DAYS(j);
            xt = xs{i}(idx,:);
            yt = ys{i}(idx,:);
            
            im{3}.mean{2}{2} = parms.a(i + parms.valididx);
            par = {meanfunc,covfunc,likfunc, xt, yt};
            hyp_race = full2one(besthyp, i + parms.valididx, parms.ncandidates, parms.nfirm);
            nz = max([-xs{i}(1,1),-DAYS(j)]);
            xstar = [(-nz+1:0)',zeros(nz,1),ones(nz,1),parms.nfirm*ones(nz,1),republican*ones(nz,1)];
            
            if isempty(xt)
                fmu = ones(nz,1)*im{3}.mean{2}{2};
                fs2 = ones(nz,1)*im{3}.mean{2}{3}^2;
            else
               [~, ~, fmu, fs2] = gp(hyp_race, im, par{:}, xstar); 
            end
            fts{i}(j) = fmu(end); s2s{i}(j) = fs2(end);
            
            fig = plot_posterior(fmu, fs2, xt(:,1), yt, xstar(:,1), i);
            year = raceinfos{i}{1};
            state = raceinfos{i}{2}{1};
            candidateName = raceinfos{i}{3};
            trueVote = raceinfos{i}{4}/100;
            plot([0], [trueVote], 'rx');
            plot_title = year + " " + state + " " + candidateName + " before days " + -DAYS(j);
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
    end

end