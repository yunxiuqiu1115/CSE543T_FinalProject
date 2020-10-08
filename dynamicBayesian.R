# setwd('/Users/yahoo/Documents/WashU/CSE515T/Code/Gaussian Process/')

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is argument: if not, use default 2018
if (length(args)==0) {
  input_str = '21'
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

output_path = 'results/brw_'
output_file <- paste(output_path, test_year, "day_", input_str,'.csv',sep='')

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

CYCLE <- c()
STATE <- c()
CANDIDATE <- c()
PMEAN <- c()
PRIOR <- c()
VOTE <- c()
RMSE <- c()
LOWER95 <- c()
UPPER95 <- c()
WIN <- c()
WINNERS <- c()
MEDIAN <- c()
NLZ <- c()

# slicing parameter
W = 4

# iterate over candidate
for(state in unique(df$state)){
  cs = unique(df[df$state==state,'Candidateidentifier'])
  votes = c()
  preds = c()
  for (c in cs) {
    data = df[as.character(df$Candidateidentifier)==c,]
    cycle = data$cycle[1]
    state = as.character(data$state[1])
    vote = data$vote[1]
    votes = c(votes, vote)
    republican = data$Republican[1]
    democratic = data$Democrat[1]
    pvi = data$pvi[1]
    experienced = data$experienced[1]
    h = predict(priorModel, data[1,])[[1]]
    PRIOR = c(PRIOR, h)
    
    data = data[data$daysLeft<=-as.numeric(input_str),]
    n_poll = nrow(data)
    if(n_poll>0 & max(data$daysLeft)<0){
      ns = as.array(data$samplesize)
      ys = as.array(data$polls)
      days = as.array(ceiling(-data$daysLeft/W))
      J = max(days)
      stan_data <- list(N=n_poll,
                        days=days,
                        ns=ns,
                        ys=ys,
                        h=h,
                        v=vote,
                        J=J)
      
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
      pred = fit_params$y
    }
    else{
      pred = rnorm(10000, mean = h, sd = 0.1)
    }
    preds = c(preds, pred)
    u = quantile(pred,probs=c(0.975),names = FALSE)
    l = quantile(pred,probs=c(0.025),names = FALSE)
    m = mean(pred)
    s = sd(pred)
    CYCLE <- c(CYCLE, cycle)
    STATE <- c(STATE, state)
    CANDIDATE <- c(CANDIDATE,substring(c, 7))
    PMEAN <- c(PMEAN, m)
    RMSE = c(RMSE, sqrt(mean((pred - rep(vote,length(pred)))^2)))
    VOTE <- c(VOTE, vote)
    MEDIAN <- c(MEDIAN, median(pred))
    LOWER95 <- c(LOWER95, l)
    UPPER95 <- c(UPPER95, u)
    if(n_poll){
      NLZ <- c(NLZ,  -log(mean(exp(fit_params$ll))))
    }
    else{
      NLZ <- c(NLZ,  -dnorm(vote,h,0.1, log=TRUE))
    }
    
  }
  
  preds <- matrix(preds, nrow = length(cs), byrow = TRUE)
  win_rates = rep(0, length(cs))
  for(k in 1:ncol(preds)){
    idx = which.max(preds[,k])
    win_rates[idx] = win_rates[idx] + 1
  }
  win_rates = win_rates / sum(win_rates)
  
  WIN <- c(WIN, win_rates)
  winners = rep(0, length(cs))
  winners[which.max(votes)] = 1
  WINNERS  <- c(WINNERS, winners)
}

# write results to csv
result <- data.frame(CYCLE,
                     STATE,
                     CANDIDATE,
                     PRIOR,
                     PMEAN,
                     VOTE,
                     LOWER95,
                     UPPER95,
                     MEDIAN,
                     WIN,
                     WINNERS,
                     RMSE,
                     NLZ)

names(result) <- tolower(names(result))

write.csv(result,output_file)