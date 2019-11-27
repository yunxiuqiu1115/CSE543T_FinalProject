function performForcasting(besthyp, xs, ys, raceinfos, plot_path, parms)
    [meanfunc, covfunc, likfunc, inffunc, prior] = model(parms);
    im = {@infPrior, inffunc, prior};
    p.method = 'LBFGS';
    p.mem = 100;
    p.verbosity = 0;
    p.length = -100;
%     liksize = size(besthyp.lik, 1);
%     covsize = size(besthyp.cov,1);
    
    DAYS = [-90,-60,-30,-14,-7,-1];
    
    for i=1:numel(xs)
        xs{i}(:,4) = parms.nfirm;
        republican = xs{i}(1,5);
%         ndays = size(xs{i},1);
%         ft = zeros(1,ndays);
%         s2t = zeros(1,ndays);
        for j = 1:numel(DAYS)
            idx = xs{i}(:,1) < DAYS(j);
            xt = xs{i}(idx,:);
            yt = ys{i}(idx);
            
            im{3}.mean{2}{2} = parms.a(i);
            par = {meanfunc,covfunc,likfunc, xt, yt};
            hyp_race = full2one(besthyp, 1, parms.ncandidates, parms.nfirm);
%             hyp = struct; hyp.mean(1:2) = [0; parms.a(i+parms.valididx)];
%             hyp_race.mean(1:2) = [0; parms.a(i+parms.valididx)];
%             mask = false(size(unwrap(hyp_race)));
%             mask(liksize+covsize+1) = 1;
%             mask(liksize+covsize+2) = 1;
%             hyp = minimize_v2(hyp, @gp, p, hyp_race, im, par{:}, mask);
%             hyp_race.mean(1:2) = hyp.mean;
            
            xstar = [(DAYS(j)+1:0)',zeros(1,-DAYS(j))',ones(1,-DAYS(j))',...
                parms.nfirm*ones(1,-DAYS(j))',republican*ones(1,-DAYS(j))'];
            
            if isempty(xt)
                [~, ~, fmu, fs2] = gp(hyp_race, inffunc, par{:}, xstar);
            else
               [~, ~, fmu, fs2] = gp(hyp_race, im, par{:}, xstar); 
            end
            % ft(j) = mu; s2t(j) = s2;
            
            fig = plot_posterior(fmu, fs2, xs{i}(:,1), ys{i}, xstar(:,1), i);
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
%         ts = xs{i}(:,1);
%         fig = figure('visible', 'off');
%         f = [ft+2*sqrt(s2t); flip(ft-2*sqrt(s2t),1)];
%         fill([ts.'; flip(ts.',1)], f, [7 7 7]/8);
%         hold on; plot(ts, ft);
    end

end