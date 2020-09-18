function [varout]=test(TYPE, CV)

%  Perform LOYO or testing for GP/LM model
%  input:
%    - TYPE: 'GP' or 'LM'
%    - CV: 1 if doing LOYO, 0 if test 2018
% 
    varout = "";
    addpath("gpml-matlab-v3.6-2015-07-07");
    addpath("utilities");
    addpath("data");
    startup;

    % define horizons
    taus = [0,7, 14,28,42,90,120];

    % define search space
    search_size = 100;
    if strcmp(TYPE, "GP")==1
        p = sobolset(3);
    else 
        p = sobolset(1);
    end
    
    % best cv index
    if CV==0
        ts = [32,32,94,46,46,36,36];
%         ts = [7, 79, 34,21, 95, 41, 93];

        for i=1:numel(taus)
            j = ts(i);
            if strcmp(TYPE, "GP")==1
                ls = p(j,1)*max(taus(i),30)+3; % 3-tau;
                os = p(j,2)/10; % 0%-10%
                lik = p(j,3)/10; % 0%-10%
                myrun(taus(i),TYPE, ls, os, lik, j, CV);
            else 
                % linear model does not have ls/os
                ls = 0;
                os = 0;
                lik = p(1,3)/10; % 0%-10%
                myrun(taus(i),TYPE, ls, os, lik, j, CV);
            end
        end
    else

        for i=1:numel(taus)
            for j=1:search_size
                if strcmp(TYPE, "GP")==1
                    ls = p(j,1)*max(taus(i),30)+3; % 3-tau;
                    os = p(j,2)/10; % 0%-10%
                    lik = p(j,3)/10; % 0%-10%
                    myrun(taus(i),TYPE, ls, os, lik, j, CV);
                else 
                    % linear model does not have ls/os
                    ls = 0;
                    os = 0;
                    lik = p(1,3)/10; % 0%-10%
                    myrun(taus(i),TYPE, ls, os, lik, j, CV);
                end
            end
        end
    end
end

function myrun(tau,type, ls, os, lik, j, CV)
    CNNdata = readData("data/CNNData.csv");
    
    if CV==0
        CNNdata2018 = readData("data/CNNData2018.csv");
        CNNdata2018(:, ["candidate_name"]) = [];
        CNNdata = vertcat(CNNdata, CNNdata2018);
        parms.test_year = 2018;
        parms.coefs = priorModel(CNNdata, parms.test_year);
    end

%     CNNdata2020 = readData("data/CNNData2020.csv");
%     CNNdata2020(:, ["candidate_name"]) = [];
%     CNNdata2020.Percentage_of_Vote_won_x = zeros(size(CNNdata2020,1),1);
%     CNNdata = vertcat(CNNdata, CNNdata2020);

    
    parms.type = type;
    parms.days = min(CNNdata.daysLeft);
    years = unique(CNNdata.cycle);
    states = unique(CNNdata.state);
    [xs, ys, raceinfos] = buildTrainCellArrays(CNNdata, years, states);

    counter = size(xs,1);
    for i=1:counter
        idx = xs{i}(:,1) <= -tau;
        xs{i} = xs{i}(idx,:);
        ys{i} = ys{i}(idx);
    end

    hyp.cov(1) = log(ls);
    hyp.cov(2) = log(os);
    [~,~,~,~, prior] = model(parms);
    sigma_ml = prior.slope(2);
    sigma_mc = prior.intercept(2);
    hyp.cov(3) = log(1/sigma_ml);
    hyp.cov(4) = log(sigma_mc);
    hyp.lik = log(lik);

    disp(type);
    disp("tau: "+tau);
    parms.tau = tau;
    parms.j = j;
    parms.plot = 0;
    
    if CV==0
        plot_path = "plots/" + type + "MargLinTre"+num2str(parms.test_year)+"_"+num2str(tau);
        if strcmp(type, "GP")==1
            [allRaces, fts, s2s] = forcastAllRaces(hyp, xs, ys, raceinfos, plot_path, parms);
        else
            [allRaces, fts, s2s] = lm(hyp, xs, ys, raceinfos, plot_path, parms);
        end
        posttrain(raceinfos,fts,s2s,allRaces,hyp, tau, parms);   
    else
        for y=1:numel(years)
            parms.test_year = years(y);
            parms.coefs = priorModel(CNNdata, parms.test_year);
            parms.type = type;
            plot_path = "plots/" + type + "MargLinTre"+num2str(parms.test_year)+"_"+num2str(tau);
            if strcmp(type, "GP")==1
                [allRaces, fts, s2s] = forcastAllRaces(hyp, xs, ys, raceinfos, plot_path, parms);
            else
                [allRaces,fts,s2s] = lm(hyp, xs, ys, raceinfos, plot_path, parms);
            end        
            posttrain(raceinfos,fts,s2s,allRaces,hyp, tau, parms);
        end
    end
end
