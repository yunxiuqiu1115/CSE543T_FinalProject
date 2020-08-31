setwd('/Users/yahoo/Documents/WashU/CSE515T/Code/Gaussian Process/')
library(rstan)
library(dplyr)

output_path = 'results/stan_dynamic18_'
input_strs = c('0',
               '14',
               '28',
               '42',
               '90',
               '120')

input_strs = c('7')

for (i in 1:length(input_strs)) {
  input_str = input_strs[i]
  print(input_str)
  output_file <- paste(output_path, input_str,'.csv',sep='')
  
  # read data
  df <- read.csv("data/CNNData1992to2018.csv")
  colnames(df)[colnames(df) == "Percentage.of.Vote.won.x"] = 'vote'
  c2v = df %>% group_by(Candidateidentifier) %>% summarise_at(vars(vote), list(mean = mean))
  df$polls <- df$numberSupport/df$samplesize
  
  df <- df[df$daysLeft<=as.numeric(input_str),]
  
  CYCLE <- c()
  STATE <- c()
  CANDIDATE <- c()
  PMEAN <- c()
  PSTD <- c()
  VOTE <- c()
  LOWER95 <- c()
  UPPER95 <- c()
  MEDIAN <- c()
  NLZ <- c()
  PRIOR <- c()
  Nout_test = 0
  
  # slicing parameter
  W = 4
  
  # iterate over candidate
  for (i in 761:761) {
    c = as.character(c2v$Candidateidentifier[i])
    data = df[as.character(df$Candidateidentifier)==c,]
    cycle = data$cycle[1]
    state = as.character(data$state[1])
    vote = c2v$mean[i]
    republican = data$Republican[1]
    democratic = data$Democrat[1]
    pvi = data$pvi[1]
    experienced = data$experienced[1]
    n_poll = nrow(data)
    days = as.array(ceiling(-data$daysLeft/W))
    ns = as.array(data$samplesize)
    ys = as.array(data$polls)
    h = (18.98808 + pvi*0.66579 + experienced*5.63548
         + republican*23.98806 + democratic*25.74237
         - pvi*republican*1.2485)/100;
    sigma_J = 0.1
    J = max(days)
    
    stan_data <- list(N=n_poll,
                      days=days,
                      ns=ns,
                      ys=ys,
                      h=h,
                      sigma_J=sigma_J,
                      J=J)
    
    # define stan model
    model <- stan_model("dynamicBayesian.stan")
    
    # train stan model
    fit <- stan(file = "dynamicBayesian.stan",
                data = stan_data, 
                warmup = 1000, 
                iter = 5000, 
                chains = 1, 
                cores = 1, 
                thin = 4,
                control=list(adapt_delta=.99, max_treedepth = 20),
                refresh=0
    )
    
    fit_params <- as.data.frame(fit)
    
    pred = fit_params$y
    u = quantile(pred,probs=c(0.975),names = FALSE)
    l = quantile(pred,probs=c(0.025),names = FALSE)
    m = mean(pred)
    s = sd(pred)
    if (vote/100>=u | vote/100<=l){
      Nout_test = Nout_test + 1
    }
    CYCLE <- c(CYCLE, cycle)
    STATE <- c(STATE, state)
    PRIOR <- c(PRIOR, h)
    CANDIDATE <- c(CANDIDATE,substring(c, 7))
    PMEAN <- c(PMEAN, m)
    PSTD <- c(PSTD, s)
    VOTE <- c(VOTE, vote)
    MEDIAN <- c(MEDIAN, median(pred))
    LOWER95 <- c(LOWER95, l)
    UPPER95 <- c(UPPER95, u)
    NLZ <- c(NLZ, (vote/100-m)^2/2/s^2 + log(s) + log(2*pi)/2)
  }
  
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
}