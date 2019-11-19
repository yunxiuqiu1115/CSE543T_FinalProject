% CNN Forecasting Toy Example
% addpath("/Users/yahoo/Documents/WashU/CSE515T/Code/gpml-matlab-v4.2-2018-06-11");
addpath("gpml-matlab-v3.6-2015-07-07");
addpath("utilities");
startup;

% reading race data from all years and states
% CNNdata = readData("CNNData.csv");
% pollthres = 40;
% CNNdata = indexPollster(CNNdata, pollthres, "Gaussian Process/CNNDataidx.csv");

CNNdata = readData("CNNDataidx.csv");

parms.mode = true;
parms.nfirm = max(CNNdata.pollsteridx);
years = unique(CNNdata.cycle);
states = unique(CNNdata.state);

% iterate cycle/state race
N = 0; nsuc = 0;
allPollVote = cell(numel(years),numel(states));
for i = 1:numel(years)
   for j = 1:numel(states)
      % success = -1 if there is race data for year/state
      % success = 0 if failling to predict, otherwise 1
      [success, predPolls, trueVotes] = buildOneRaceModel(CNNdata, years(i), states(j), parms);
      N = N + (success>=0); nsuc = nsuc + max([0, success]);
      if success>=0, allPollVote{i, j} = {predPolls, trueVotes, success}; end
   end
end

disp(N + " races run.");
disp(nsuc + " successful predictions.");
disp("Prediction rate: " + nsuc/N);

% post training
a = [];
b = [];
errors = [];
for i = 1:numel(years)
   for j = 1:numel(states)
      if isempty(allPollVote{i,j}), continue; end
      predPolls = allPollVote{i, j}{1};
      trueVotes = allPollVote{i, j}{2};
      a = [a, predPolls'];
      b = [b, trueVotes'/100];
      errors = [errors, (predPolls - trueVotes/100).'];
   end
end

histogram(abs(errors));
xlabel("absolute error");
ylabel("freq");

wp = [3,1,2,2,0,3,3,0,1,2,1,2,3];
cp = [34,33,33,34,33,33,34,33,33,34,33,33,34] - wp;
h = bar(years, [wp;cp]);
l{1} = 'wrong';
l{2} = 'correct';
legend(h,l);
xlabel("year");
ylabel("#predictions");
title("#predictions over year");
