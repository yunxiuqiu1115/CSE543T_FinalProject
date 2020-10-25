# setwd('/Users/yahoo/Documents/WashU/CSE515T/Code/Gaussian Process/')

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is argument: if not, use default 2018
if (length(args)==0) {
  input_str = '0'
  test_year = 2016
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
df <- read.csv("data/CNNdata1992-2016.csv")
colnames(df)[colnames(df) == "Percentage_of_Vote_won_x"] = 'vote'
df$vote = df$vote/100

polls = unique(df[,c("cycle", "state", "samplesize", "daysLeft","pollster")])
for (i in 1:nrow(polls)) {
  poll = polls[i,]
  cycle = poll$cycle[1]
  state = poll$state[1]
  samplesize = poll$samplesize[1]
  daysLeft = poll$daysLeft[1]
  pollster = poll$pollster[1]
  mask = df$cycle==cycle & df$state==state & df$samplesize==samplesize & df$daysLeft==daysLeft & df$pollster==pollster
  if(daysLeft==1){
    mask = df$cycle==cycle & df$state==state & df$samplesize==samplesize & df$daysLeft==daysLeft
  }
  if(length(unique(df[mask,c("Candidateidentifier")]))!=length(df[mask,c("samplesize")])){
    print(df[mask,])
  }
  df[mask,c("samplesize")] = sum(df[mask,c("numberSupport")])
  df[mask, c("vote")] = df[mask, c("vote")]/sum(df[mask, c("vote")])
}

df2018 <- read.csv("data/CNNdata2018.csv")
colnames(df2018)[colnames(df2018) == "Percentage.of.Vote.won.x"] = 'vote'
df2018$vote = df2018$vote/100
df2018 = subset(df2018, select=-c(candidate_name))

polls = unique(df2018[,c("cycle", "state", "samplesize", "daysLeft","pollster")])
for (i in 1:nrow(polls)) {
  poll = polls[i,]
  cycle = poll$cycle[1]
  state = poll$state[1]
  samplesize = poll$samplesize[1]
  daysLeft = poll$daysLeft[1]
  pollster = poll$pollster[1]
  mask = df2018$cycle==cycle & df2018$state==state & df2018$samplesize==samplesize & df2018$daysLeft==daysLeft & df2018$pollster==pollster
  if(daysLeft==1){
    mask = df2018$cycle==cycle & df2018$state==state & df2018$samplesize==samplesize & df2018$daysLeft==daysLeft
  }
  df2018[mask, c("vote")] = df2018[mask, c("vote")]/sum(df2018[mask, c("vote")])
}

df = rbind(df, df2018)

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
PRECISION <- c()

# slicing parameter
W = 4
precision = 1/0.1^2

# iterate over candidate
# unique(df$state)
states = unique(df$state)
I = 0
NS <- c()
YS <- c()
STATES <- c()
DAYS <- c()
K = 0
for(state in states){
  cs = unique(df[df$state==state,'Candidateidentifier'])
  I = I + 1
  for (c in cs) {
    data = df[as.character(df$Candidateidentifier)==c,]
    cycle = data$cycle[1]
    state = as.character(data$state[1])
    vote = data$vote[1]
    republican = data$Republican[1]
    democratic = data$Democrat[1]
    if(democratic==0){
      next
    }
    if(c=="2016CASanchez"){
      next
    }
    
    pvi = data$pvi[1]
    experienced = data$experienced[1]
    h = predict(priorModel, data[1,])[[1]]

    data = data[data$daysLeft<=-as.numeric(input_str),]
    n_poll = nrow(data)
    
    if(n_poll>0 & max(data$daysLeft)<0){
      ns = as.array(data$samplesize)
      ys = as.array(data$numberSupport)
      days = as.array(ceiling(-data$daysLeft/W))
      k = length(ns)
      NS <- c(NS, ns)
      YS <- c(YS, ys)
      DAYS <- c(DAYS, days)
      K = K + k
      STATES <- c(STATES, rep(I, k))
      J = max(days)
      PRIOR = c(PRIOR, h)
      PRECISION = c(PRECISION, precision)
      VOTE <- c(VOTE, vote)
    }
   
    # else{
    #   pred = rnorm(10000, mean = h, sd = 0.1)
    # }
    # preds = c(preds, pred)
    # u = quantile(pred,probs=c(0.975),names = FALSE)
    # l = quantile(pred,probs=c(0.025),names = FALSE)
    # m = mean(pred)
    # s = sd(pred)
    # CYCLE <- c(CYCLE, cycle)
    # STATE <- c(STATE, state)
    # CANDIDATE <- c(CANDIDATE,substring(c, 7))
    # PMEAN <- c(PMEAN, m)
    # RMSE = c(RMSE, sqrt(mean((pred - rep(vote,length(pred)))^2)))
    # MEDIAN <- c(MEDIAN, median(pred))
    # LOWER95 <- c(LOWER95, l)
    # UPPER95 <- c(UPPER95, u)
    # if(n_poll){
    #   NLZ <- c(NLZ,  -log(mean(exp(fit_params$ll))))
    # }
    # else{
    #   NLZ <- c(NLZ,  -dnorm(vote,h,0.1, log=TRUE))
    # }
    # 
  }
}

J = max(DAYS)
DAYS = J + 1 - DAYS

stan_data <- list(I=I,
                  J=J,
                  K=K,
                  n=NS,
                  y=YS,
                  state=STATES,
                  day=DAYS,
                  h=PRIOR,
                  tau=PRECISION,
                  v=VOTE)

# define stan model
# model <- stan_model("dynamicBayesian.stan")

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

# write.csv(result,output_file)


# tmp = c()
# upper = c()
# lower = c()
# for (i in n_poll:1){
#   pred = fit_params[[paste("ps[", as.character(i) , "]", sep="")]]
#   tmp = c(tmp, quantile(pred, 0.5, na.rm = TRUE))
#   upper = c(upper, quantile(pred, 0.975, na.rm = TRUE))
#   lower = c(lower, quantile(pred, 0.025, na.rm = TRUE))
# }
# 
# plot(sort(-W*days), ys, pch=19)
# lines(sort(-W*days), tmp, type="l", ylim=c(0.2,0.8))
# lines(sort(-W*days), upper,type="l")
# lines(sort(-W*days), lower,type="l")

for (i in 1:I){
  tmp = c()
  upper = c()
  lower = c()
  for (j in 1:J){
    pred = fit_params[[paste("beta[",as.character(i) ,",", as.character(j) , "]", sep="")]]
    tmp = c(tmp, quantile(pred, 0.5, na.rm = TRUE))
    upper = c(upper, quantile(pred, 0.975, na.rm = TRUE))
    lower = c(lower, quantile(pred, 0.025, na.rm = TRUE))
  }
  jpeg(paste("plots/brw/",as.character(i), ".jpg", sep=""), width = 1000, height = 1000)
  plot(-(J:1),1/(1+exp(-tmp)), type="l",ylim=c(0,1))
  lines(-(J:1),1/(1+exp(-upper)),type="l")
  lines(-(J:1),1/(1+exp(-lower)),type="l")
  points(-(J+1-DAYS[STATES==i]), rev(YS[STATES==i]/NS[STATES==i]))
  points(0, VOTE[i], col="blue",pch=20, cex=2)
  dev.off()
}


