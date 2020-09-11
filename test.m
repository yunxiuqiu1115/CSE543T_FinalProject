% rng('default');
% x = (1:5)';
% [meanfunc, covfunc, likfunc, inffunc, prior] = model(parms);
% prior = rmfield(prior, 'mean');
% im = {@infPrior, inffunc, prior};
% par = {meanfunc,covfunc,likfunc, xs, ys};
% hyperparameters = sample_separate_prior(prior, parms, counter, 1);
% num_samples = size(xs,1); 
% hyp.mean = hyp.mean(2*num_samples+1:end);
% hyp.cov = hyperparameters.cov;
% hyp.lik = hyperparameters.lik;
% means = reshape(hyperparameters.mean, num_samples,2);
% 
% d = 1e-6;
% [nlZ, dnlZ] = gp_independent(hyp, prior, means, im, par{:});
% 
% new_hyp = hyp;
% new_hyp.lik = new_hyp.lik + d;
% [new_nlZ, new_dnlZ] = gp_independent(new_hyp, prior, means, im, par{:});
% 
% fprintf('dnlZ.lik: %0.6f, ((new_nlZ - nlZ) / d): %0.6f\n', dnlZ.lik, (new_nlZ - nlZ) / d);

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

% means = [1;-3];
% standevs = [1;2];
% logpdf = @(theta)normalDistGrad(theta,means,standevs);
% 
% startpoint = randn(2,1);
% 
% smp = hmcSampler(logpdf,startpoint);

addpath("gpml-matlab-v3.6-2015-07-07");
addpath("utilities");
addpath("data");
startup;

taus = [0, 7, 14,28,42,90,120];

p = sobolset(3,'Skip',1e3);
search_size = 100;

% ts = [65,21,63,86,36,36];
ts = [89,89,66,12,12,36,36];

taus = [42];
ts = [12];


for i=1:numel(taus)
    j = ts(i);
    ls = p(j,1)*max(taus(i),30)+3; % 3-tau;
    os = p(j,2)/10; % 0%-10%
    lik = p(j,3)/10; % 0%-10%
    myrun(taus(i),"GP", ls, os, lik, j);
end


% for i=1:numel(taus)
%     for j=1:search_size
%         ls = p(j,1)*max(taus(i),30)+3; % 3-tau;
%         os = p(j,2)/10; % 0%-10%
%         lik = p(j,3)/10; % 0%-10%
%         myrun(taus(i),"GP", ls, os, lik, j);
% %         myrun(taus(i),"LM", ls, os, lik, j);
%     end
% end


% for i=1:numel(taus)
%     LL = LLs(i,:);
%     lss = [];
%     oss = [];
%     liks = [];
%     for j=1:search_size
%         ls = p(j,1)*max(taus(i),30)+3; % 3-tau;
%         os = p(j,2)/10; % 0%-10%
%         lik = p(j,3)/10; % 0%-10%
%         lss= [lss, ls];
%         oss = [oss,os];
%         liks = [liks,lik];
%     end
%     fig = figure(i*100+j);
%     scatter(lss, oss, 30, LL,'filled');
%     title(strcat('os vs ls at horizon',num2str(taus(i))));
%     colorbar;
%     saveas(fig, "plots/CV/ls-"+num2str(taus(i))+".jpg");
%     close;
%     
%     fig = figure();
%     scatter(oss, liks, 30,  LL,'filled');
%     title(strcat('lik vs os at horizon',num2str(taus(i))));
%     colorbar;
%     saveas(fig, "plots/CV/os-"+num2str(taus(i))+".jpg");
%     close;
%     
%     fig = figure();
%     scatter(liks, lss, 30, LL,'filled');
%     colorbar;
%     title(strcat('ls vs lik at horizon',num2str(taus(i))));
%     saveas(fig, "plots/CV/lik-"+num2str(taus(i))+".jpg");
%     close;
% end

% diary('off');

function myrun(tau,type, ls, os, lik, j)
    load("models/model.mat");
%     CNNdata = readData("data/CNNData1992to2018.csv");
%     [CNNdata, pollster2idx] = indexPollster(CNNdata, pollthres);
    CNNdata = readData("data/CNNData.csv");
    [CNNdata,pollster2idx] = indexPollster(CNNdata, 50);
    CNNdata2018 = readData("data/CNNData2018.csv");
    
    CNNdata2018 = indexPollster(CNNdata2018, pollster2idx);
    CNNdata2018(:, ["candidate_name"]) = [];
    CNNdata = vertcat(CNNdata, CNNdata2018);

    CNNdata2020 = readData("data/CNNData2020.csv");
    
    CNNdata2020 = indexPollster(CNNdata2020, pollster2idx);
    CNNdata2020(:, ["candidate_name"]) = [];
    CNNdata2020.Percentage_of_Vote_won_x = zeros(size(CNNdata2020,1),1);
    CNNdata = vertcat(CNNdata, CNNdata2020);
    
    parms.test_year = 2020;
    parms.coefs = priorModel(CNNdata, parms.test_year);
    parms.type = type;

    parms.nfirm = max(CNNdata.pollsteridx);
    parms.days = min(CNNdata.daysLeft);
    parms.nfirm = 0;
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
    hyp.lik = log(lik);

    disp(type);

    disp("tau: "+tau);
    parms.days = min(CNNdata.daysLeft);
    parms.tau = tau;
    parms.j = j;
    
%     NOPOLLSTATES = ['Arkansas'; 'Delaware'; 'Idaho'; 'Louisiana' 'Massachusetts';...
%         'Minnesota'; 'Nebraska'; 'Oregon'; 'Rhode Island';...
%         'South Dakota'; 'Virginia'; 'West Virginia'; 'Wyoming'];
    
    
    
    plot_path = "plots/GPMargLinTre"+num2str(parms.test_year)+"_"+num2str(tau);
    [allRaces, fts, s2s] = forcastAllRaces(hyp, xs, ys, raceinfos, plot_path, parms);
    posttrain(raceinfos,fts,s2s,allRaces,hyp, tau, parms);
    
%      fprintf('ls: %0.4f, os: %0.4f, lik: %0.4f\n',ls, os, lik);
   
%     for y=1:numel(years)
%         parms.test_year = years(y);
%         parms.coefs = priorModel(CNNdata, parms.test_year);
%         parms.type = type;
%         
%         if strcmp(type, "GP")==1
%             plot_path = "plots/GPMargLinTre"+num2str(years(y))+"_"+num2str(tau);
%             [allRaces, fts, s2s] = forcastAllRaces(hyp, xs, ys, raceinfos, plot_path, parms);
%         else
%             plot_path = "plots/LMMargLinTre"+num2str(years(y))+"_"+num2str(tau);
%             [allRaces,fts,s2s] = lm(hyp, xs, ys, raceinfos, plot_path, parms);
%         end        
%         
% %         'Arkansas'; 'Delaware'; 'Idaho'; 'Louisiana'; 'Minnesota'; 
% %         'Nebraska'; 'Oregon'; 'Rhode Island'; 'South Dakota'; 'Virginia'; 'West Virginia'; 'Wyoming';
%         posttrain(raceinfos,fts,s2s,allRaces,hyp, tau, parms);
%     end
end


function [lpdf,glpdf] = normalDistGrad(X,Mu,Sigma)
Z = (X - Mu)./Sigma;
lpdf = sum(-log(Sigma) - .5*log(2*pi) - .5*(Z.^2));
glpdf = -Z./Sigma;
end


% plot(1:numel(nlZs),nlZs,'k','LineWidth',1);
% hold on;
% xlabel("number of iterations");
% ylabel("Averaged nlz");
% g_x=[1:1:numel(nlZs)]; % user defined grid Y [start:spaces:end]
% g_y=[nlZs(end)-0.01:0.005:nlZs(1)+0.01]; % user defined grid X [start:spaces:end]
% for i=1:length(g_x)
%    plot([g_x(i) g_x(i)],[g_y(1) g_y(end)],'k:'); %y grid lines
%    hold on;  
% end
% for i=1:length(g_y)
%    plot([g_x(1) g_x(end)],[g_y(i) g_y(i)],'k:'); %x grid lines
%    hold on;  
% end
