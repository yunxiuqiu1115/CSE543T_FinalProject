% rng('default');
% x = (1:5)';
% y = randn(size(x));
% xs = (6)';
% ys = randn(size(xs));
% mean = {@meanConst};
% inf = @infExact;
% cov = {@covSEiso};
% lik = @likGauss;
% 
% theta.mean = 1;
% theta.cov  = [log(1); log(1)];
% theta.lik  = log(0.01);
% [nlZ, dnlZ] = gp_last(theta, inf, mean, cov, lik, x, y, xs, ys);
% d = 1e-6;
% new_theta = theta;
% new_theta.mean = theta.mean + d;
% [new_nlZ] = gp_last(new_theta, inf, mean, cov, lik, x, y, xs, ys);
% fprintf('dnlZ.mean: %0.6f, ((new_nlZ - nlZ) / d): %0.6f\n', dnlZ.mean, (new_nlZ - nlZ) / d);
% 
% d = 1e-6;
% new_theta = theta;
% new_theta.cov(1) = theta.cov(1) + d;
% [new_nlZ] = gp_last(new_theta, inf, mean, cov, lik, x, y, xs, ys);
% fprintf('dnlZ.cov(1): %0.6f, ((new_nlZ - nlZ) / d): %0.6f\n', dnlZ.cov(1), (new_nlZ - nlZ) / d);
% 
% d = 1e-6;
% new_theta = theta;
% new_theta.cov(2) = theta.cov(2) + d;
% [new_nlZ] = gp_last(new_theta, inf, mean, cov, lik, x, y, xs, ys);
% fprintf('dnlZ.cov(2): %0.6f, ((new_nlZ - nlZ) / d): %0.6f\n', dnlZ.cov(2), (new_nlZ - nlZ) / d);
% 
% d = 1e-6;
% new_theta = theta;
% new_theta.lik = theta.lik + d;
% [new_nlZ] = gp_last(new_theta, inf, mean, cov, lik, x, y, xs, ys);
% fprintf('dnlZ.lik: %0.6f, ((new_nlZ - nlZ) / d): %0.6f\n', dnlZ.lik, (new_nlZ - nlZ) / d);

% diary('Log.txt');

taus = [0,14,28,42,90,120];

for i=1:numel(taus)
    myrun(taus(i),"model");
   myrun(taus(i),"last");
end

% diary('off');

function myrun(tau,type)
%     taus = [42,28,14,0];
%     length_scales = [4.5005,4.4987,4.0310,3.5473];
    if strcmp(type, "model")==1
%         load("models/all" + tau + ".mat");
        load("models/All2016InitEstiThres50Iter20Seed1.mat");
        method = 'all';
%         load('models/All2016InitEstiThres50Iter20Seed1.mat');
%         method = 'initest';
    else
        load(MPLV(tau));
        method = 'last';
    end
%     idx = taus==tau;
%     hyp.cov(1) = length_scales(idx);
    CNNdata = readData("data/CNNData.csv");
    [~,pollster2idx] = indexPollster(CNNdata, pollthres);
    CNNdata = readData("data/CNNData1992to2018.csv");
    CNNdata = indexPollster(CNNdata, pollster2idx);
    years = unique(CNNdata.cycle);
    [xs, ys, raceinfos] = buildTrainCellArrays(CNNdata, years, states);
    counter = size(xs,1);
    for i=1:counter
        if raceinfos{i}{1}>=2016
            idx = xs{i}(:,1) <= -tau;
            xs{i} = xs{i}(idx,:);
            ys{i} = ys{i}(idx);
        end
    end

    disp(type);
%     plot_path = strcat("plots/All2016InitEst",num2str(tau));
    
    disp("tau: "+tau);
    parms.days = min(CNNdata.daysLeft);
    [allRaces, fts, s2s] = forcastAllRaces(hyp, xs, ys, raceinfos, plot_path, parms);
    
    posttrain(raceinfos,fts,s2s,allRaces,hyp, tau, method);
end

function f=MPLV(t)
    f = "models/last" + t + "-1.mat";
end

