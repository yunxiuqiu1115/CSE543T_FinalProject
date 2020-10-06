# setwd('/Users/yahoo/Documents/WashU/CSE515T/Code/Gaussian Process/')

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is argument: if not, use default 2018
if (length(args)==0) {
  input_str = '0'
  test_year = 2018
}
if (length(args)==1){
  # define the test year
  input_str = args[1]
  test_year = 2018
}
if (length(args)==2){
  # define the test year
  input_str = args[1]
  test_year = as.double(args[2])
}

library(rstan)

output_path = 'results/stan_dynamic18_'
output_file <- paste(output_path, input_str,'.csv',sep='')

computeh = function(df, test_year){
  if(test_year==2018){
    data = df[df$cycle!=test_year, ]
  }
  else{
    data = df[df$cycle!=test_year & df$cycle!=2018, ]
  }
  model <- lm(vote ~ pvi + experienced + Democrat + Republican + pvi*Republican, data) 
  return(model)
}


# read data
df <- read.csv("data/CNNdata1992-2018.csv")
colnames(df)[colnames(df) == "Percentage_of_Vote_won_x"] = 'vote'
df$vote = df$vote/100
# Candidateidentifier = c()
# Vote = c()
# for (c in unique(df$Candidateidentifier)){
#   Candidateidentifier = c(Candidateidentifier, c)
#   Vote = c(Vote, df[df$Candidateidentifier==c, c('vote')][1])
# }
# c2v = data.frame(Candidateidentifier, Vote)
# c2v = df %>% group_by(Candidateidentifier) %>% summarise_at(vars(vote), list(mean = mean))
df$polls <- df$numberSupport/df$samplesize

priorModel = computeh(df, test_year)

df <- df[df$cycle==test_year, ]
df <- df[df$daysLeft<=as.numeric(input_str),]

CYCLE <- c()
STATE <- c()
CANDIDATE <- c()
PMEAN <- c()
VOTE <- c()
RMSE <- c()
LOWER95 <- c()
PRIOR <- c()
UPPER95 <- c()
WIN <- c()
WINNERS <- c()
MEDIAN <- c()
NLZ <- c()

correct_predictions <- 0
Nout_test <- 0
Nout <- 0

# slicing parameter
W = 4

n_candidate = length(unique(df$Candidateidentifier))
N_polls = c()
for(state in unique(df$state)){
  cs = unique(df[df$state==state,'Candidateidentifier'])
  for (c in cs){
    data = df[as.character(df$Candidateidentifier)==c,]
    N_polls = c(N_polls, nrow(data))
  }
}
N_max = max(N_polls)
stan_days = matrix(0,n_candidate,N_max)
stan_ns = matrix(0,n_candidate,N_max)
stan_ys = matrix(0,n_candidate,N_max)
stan_h = rep(0, n_candidate)
stan_J = rep(0, n_candidate)

counter = 1
# iterate over candidate
for(state in unique(df$state)){
  cs = unique(df[df$state==state,'Candidateidentifier'])
  for (c in cs) {
    data = df[as.character(df$Candidateidentifier)==c,]
    cycle = data$cycle[1]
    state = as.character(data$state[1])
    vote = data$vote[1]
    republican = data$Republican[1]
    democratic = data$Democrat[1]
    pvi = data$pvi[1]
    experienced = data$experienced[1]
    n_poll = nrow(data)
    days = as.array(ceiling(-data$daysLeft/W))
    ns = as.array(data$samplesize)
    ys = as.array(data$polls)
    h = predict(priorModel, data[1,])[[1]]
    J = max(days)
    stan_days[counter,1:length(days)] = days
    stan_ns[counter,1:length(ns)] = ns
    stan_ys[counter,1:length(ys)] = ys
    stan_h[counter] = h
    stan_J[counter] = J
    counter = counter + 1
  }
}

J_max = max(stan_J)

stan_data <- list(n_candidate=n_candidate,
                  N_max=N_max,
                  N=N_polls,
                  days=stan_days,
                  ns=stan_ns,
                  ys=stan_ys,
                  h=stan_h,
                  J=stan_J,
                  J_max=J_max)

# define stan model
model <- stan_model("dynamicBayesian.stan")

# train stan model
fit <- stan(file = "dynamicBayesian.stan",
            data = stan_data, 
            warmup = 2000, 
            iter = 10000, 
            chains = 1, 
            cores = 1, 
            thin = 4,
            control=list(adapt_delta=.99, max_treedepth = 20),
            seed = 1,
            refresh=0
)

fit_params <- as.data.frame(fit)

exit;

pred = fit_params$y
preds = c(preds, pred)
u = quantile(pred,probs=c(0.975),names = FALSE)
l = quantile(pred,probs=c(0.025),names = FALSE)
m = mean(pred)
s = sd(pred)
if (vote>=u | vote<=l){
  Nout_test = Nout_test + 1
}
CYCLE <- c(CYCLE, cycle)
STATE <- c(STATE, state)
PRIOR <- c(PRIOR, h)
CANDIDATE <- c(CANDIDATE,substring(c, 7))
PMEAN <- c(PMEAN, m)
VOTE <- c(VOTE, vote)
MEDIAN <- c(MEDIAN, median(pred))
LOWER95 <- c(LOWER95, l)
UPPER95 <- c(UPPER95, u)
NLZ <- c(NLZ, (vote/100-m)^2/2/s^2 + log(s) + log(2*pi)/2)


correct_predictions = 0
n_race = 0
for (cycle in unique(CYCLE)) {
  for (state in unique(STATE)) {
    idx = intersect(which(STATE %in% c(state)),  which(CYCLE %in% c(cycle)))
    n_race = n_race + 1
    if (which.max(PMEAN[idx])==which.max(VOTE[idx])){
      correct_predictions = correct_predictions + 1
    }
    else{
      print("Wrong prediction:")
      print(cycle)
      print(state)
    }
  }
}

print(paste("Correct predictions: ",correct_predictions,correct_predictions/n_race))

print(paste("Correlation: ",cor(PMEAN, VOTE)))

print(paste("RSME: ",sqrt(mean((PMEAN- VOTE/100)^2))))

print(paste("Ratio in 95% : ",1-Nout_test/73))

print(paste("Predictive averaged nlZ: ",mean(NLZ)))

print(paste("Mean of predictive std: ",mean(PSTD)))

print(paste("Median of predictive std: ",median(PSTD)))

print(paste("Std of predictive std: ",sd(PSTD)))

# write results to csv
result <- data.frame(CYCLE,
                     STATE,
                     CANDIDATE,
                     PRIOR,
                     PMEAN,
                     PSTD,
                     VOTE,
                     LOWER95,
                     UPPER95,
                     MEDIAN)

names(result) <- tolower(names(result))

# write.csv(result,output_file)