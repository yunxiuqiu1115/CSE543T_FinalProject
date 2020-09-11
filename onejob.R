#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input_str, cv_year, type).n", call.=FALSE)
}

input_str = args[1]
cv_year = args[2]
TYPE = args[3]

search_size = 10

# library(rstan)

averaged_nlZs = c() 

for (b in 1:search_size){
  averaged_nlZs = c(averaged_nlZs, b)
  
  input_file = paste('results/LOO', TYPE, '_' , cv_year, 'day', input_str, '_', b ,'.csv',sep='')
  output_file = paste('nlZs/', TYPE, '_' , cv_year, 'day', input_str, '_', b,'.csv',sep='')
  next;
  # loading data
  data <- read.csv(input_file)
  
  data <- data[data$cycle!=2016 | data$state!='Louisiana' | data$candidate!='Flemsing',]
  data <- data[data$cycle!=2020 | data$state!='Georgia' | data$candidate!='Loeffler',]
  data <- data[data$cycle!=2020 | data$state!='Georgia' | data$candidate!='Tarver',]
  
  data_test <- data[(data$cycle==cv_year),]
  data <- data[(data$cycle!=cv_year & data$cycle!=2018 & data$cycle!=2020),]
  
  # data <- data[(data$cycle!=2016),]
  
  cycles <- unique(data$cycle)
  states <- union(unique(data$state), unique(data_test$state))
  
  C <- 4
  
  # define variables
  metadata <- list()
  year_idx <- c()
  mu <- list()
  sigma <- list()
  y <- list()
  nc <- c()
  pvi <- list()
  party <- list()
  experienced <- list()
  counter <- 0
  
  # iterate over races
  for (cycle in cycles) {
    for (state in states) {
      pmu = data[data$cycle==cycle & data$state==state,c("posteriormean")]
      pstd = data[data$cycle==cycle & data$state==state,c("posteriorstd")]
      vote = data[data$cycle==cycle & data$state==state,c("vote")]
      pvi_ = data[data$cycle==cycle & data$state==state,c("pvi")]
      party_ = data[data$cycle==cycle & data$state==state,c("party")]
      experienced_ = data[data$cycle==cycle & data$state==state,c("experienced")]
      
      if(length(pmu)){
        counter <- counter + 1
        metadata[[counter]] = c(cycle, state)
        mu[[counter]] = pmu
        sigma[[counter]] = pstd
        y[[counter]] = vote / 100
        pvi[[counter]] = pvi_
        party[[counter]] = party_
        experienced[[counter]] = experienced_
        nc = c(nc, length(vote))
        year_idx = c(year_idx, (cycle-1990)/2)
      }
    }
  }
  
  # build stan data
  stan_mu <- matrix(0,counter,C)
  stan_sigma <- matrix(0.0001,counter,C)
  stan_y <- matrix(0,counter,C)
  stan_pvi <- matrix(0,counter,C)
  stan_party <- matrix(0,counter,C)
  stan_experienced <- matrix(0,counter,C)
  
  for (i in 1:counter) {
    stan_mu[i,1:nc[i]] <- mu[[i]]
    stan_sigma[i,1:nc[i]] = sigma[[i]]
    stan_y[i,1:nc[i]] = y[[i]]
    stan_y[i,] = stan_y[i,]/sum(stan_y[i,])
    stan_pvi[i,1:nc[i]] = pvi[[i]]
    stan_party[i,1:nc[i]] = party[[i]]
    stan_experienced[i, 1:nc[i]] = experienced[[i]]
  }
  
  idx2 <- c()
  idx3 <- c()
  idx4 <- c()
  
  for (i in 1:counter) {
    tmp = sum(stan_mu[i,]!=0)
    if(tmp==2) idx2 = c(idx2, i)
    if(tmp==3) idx3 = c(idx3, i)
    if(tmp==4) idx4 = c(idx4, i)
  }
  
  # test data
  test_metadata <- list()
  test_year_idx <- c()
  test_mu <- list()
  test_sigma <- list()
  test_y <- list()
  test_nc <- c()
  test_pvi <- list()
  test_party <- list()
  test_experienced <- list()
  test_counter <- 0
  
  # iterate over races
  for (cycle in unique(data_test$cycle)){
    for (state in states) {
      pmu = data_test[(data_test$state==state & data_test$cycle==cycle),c("posteriormean")]
      pstd = data_test[data_test$state==state & data_test$cycle==cycle,c("posteriorstd")]
      vote = data_test[data_test$state==state & data_test$cycle==cycle,c("vote")]
      pvi_ = data_test[data_test$state==state & data_test$cycle==cycle,c("pvi")]
      party_ = data_test[data_test$state==state & data_test$cycle==cycle,c("party")]
      experienced_ = data_test[data_test$state==state & data_test$cycle==cycle,c("experienced")]
      if(length(pmu)){
        test_counter = test_counter + 1
        test_metadata[[test_counter]] = c(cycle,state)
        test_mu[[test_counter]] = pmu
        test_sigma[[test_counter]] = pstd
        test_y[[test_counter]] = vote / 100
        test_pvi[[test_counter]] = pvi_
        test_party[[test_counter]] = party_
        test_experienced[[test_counter]] = experienced_
        test_nc = c(test_nc, length(vote))
        test_year_idx = c(test_year_idx, (cycle-1990)/2)
      }
    }
  }
  
  # build stan data
  test_stan_mu <- matrix(0,test_counter,C)
  test_stan_sigma <- matrix(0.0001,test_counter,C)
  test_stan_y <- matrix(0.0001,test_counter,C)
  test_stan_pvi <- matrix(0,test_counter,C)
  test_stan_party <- matrix(0,test_counter,C)
  test_stan_experienced <- matrix(0,test_counter,C)
  
  for (i in 1:test_counter) {
    test_stan_mu[i,1:test_nc[i]] = test_mu[[i]]
    test_stan_sigma[i,1:test_nc[i]] = test_sigma[[i]]
    test_stan_y[i,1:test_nc[i]] = test_y[[i]]
    test_stan_y[i,] = test_stan_y[i,]/sum(test_stan_y[i,])
    test_stan_pvi[i,1:test_nc[i]] = test_pvi[[i]]
    test_stan_party[i,1:test_nc[i]] = test_party[[i]]
    test_stan_experienced[i, 1:test_nc[i]] = test_experienced[[i]]
  }
  
  
  test_idx2 <- c()
  test_idx3 <- c()
  test_idx4 <- c()
  
  for (i in 1:test_counter) {
    tmp = sum(test_stan_mu[i,]!=0)
    if(tmp==2) test_idx2 = c(test_idx2, i)
    if(tmp==3) test_idx3 = c(test_idx3, i)
    if(tmp==4) test_idx4 = c(test_idx4, i)
  }
  
  # define stan data structure
  stan_data <- list(N2 = length(idx2), 
                    mu2 = stan_mu[idx2,1:2], 
                    sigma2 = stan_sigma[idx2,1:2],
                    y2 = stan_y[idx2,1:2],
                    pvi2 = stan_pvi[idx2,1:2],
                    party2 = stan_party[idx2,1:2],
                    experienced2 = stan_experienced[idx2,1:2],
                    year_idx2 = year_idx[idx2],
                    N3 = length(idx3), 
                    mu3 = stan_mu[idx3,1:3], 
                    sigma3 = stan_sigma[idx3,1:3],
                    y3 = stan_y[idx3,1:3],
                    pvi3 = stan_pvi[idx3,1:3],
                    party3 = stan_party[idx3,1:3],
                    experienced3 = stan_experienced[idx3,1:3],
                    year_idx3 = year_idx[idx3],
                    N4 = length(idx4), 
                    mu4 = stan_mu[idx4,1:4], 
                    sigma4 = stan_sigma[idx4,1:4],
                    y4 = stan_y[idx4,1:4],
                    pvi4 = stan_pvi[idx4,1:4],
                    party4 = stan_party[idx4,1:4],
                    experienced4 = stan_experienced[idx4,1:4],
                    year_idx4 = year_idx[idx4],
                    test_N2 = length(test_idx2), 
                    test_mu2 = test_stan_mu[test_idx2,1:2], 
                    test_sigma2 = test_stan_sigma[test_idx2,1:2],
                    test_pvi2 = test_stan_pvi[test_idx2,1:2],
                    test_party2 = test_stan_party[test_idx2,1:2],
                    test_experienced2 = test_stan_experienced[test_idx2,1:2],
                    test_year_idx2 = test_year_idx[test_idx2],
                    test_N3 = length(test_idx3), 
                    test_mu3 = matrix(test_stan_mu[test_idx3,1:3],ncol=3,byrow = FALSE), 
                    test_sigma3 = matrix(test_stan_sigma[test_idx3,1:3],ncol=3,byrow = FALSE),
                    test_pvi3 = matrix(test_stan_pvi[test_idx3,1:3],ncol=3,byrow = FALSE),
                    test_party3 = matrix(test_stan_party[test_idx3,1:3],ncol=3,byrow = FALSE),
                    test_experienced3 = matrix(test_stan_experienced[test_idx3,1:3],ncol=3,byrow = FALSE),
                    test_year_idx3 = array(test_year_idx[test_idx3]),
                    test_N4 = length(test_idx4), 
                    test_mu4 = matrix(test_stan_mu[test_idx4,1:4],ncol=4,byrow = FALSE), 
                    test_sigma4 = matrix(test_stan_sigma[test_idx4,1:4],ncol=4,byrow = FALSE),
                    test_pvi4 = matrix(test_stan_pvi[test_idx4,1:4],ncol=4,byrow = FALSE),
                    test_party4 = matrix(test_stan_party[test_idx4,1:4],ncol=4,byrow = FALSE),
                    test_experienced4 = matrix(test_stan_experienced[test_idx4,1:4],ncol=4,byrow = FALSE),
                    test_year_idx4 = array(test_year_idx[test_idx4]),
                    max_year_idx = max(c(year_idx, test_year_idx)))
  
  # define stan model
  model <- stan_model("model.stan")
  
  # train stan model
  fit <- stan(file = "model.stan",
              data = stan_data, 
              warmup = 500, 
              iter = 3000, 
              chains = 3, 
              cores = 3, 
              thin = 4,
              control=list(adapt_delta=.98, max_treedepth = 15),
              refresh=0
  )
  
  fit_params <- as.data.frame(fit)
  
  # within 95% CI
  flags <- matrix(0, test_counter, C)
  
  CYCLE <- c()
  STATE <- c()
  CANDIDATE <- c()
  POSTERIORMEAN <- c()
  POSTERIORSTD <- c()
  PMEAN <- c()
  PSTD <- c()
  VOTE <- c()
  NORM_VOTE <- c()
  LOWER95 <- c()
  UPPER95 <- c()
  WIN <- c()
  MEDIAN <- c()
  NLZ <- c()
  
  correct_predictions <- 0
  Nout_test <- 0
  Nout <- 0
  
  
  for(i in 1:length(test_idx2)) {
    cycle = test_metadata[[test_idx2[i]]][1]
    state = test_metadata[[test_idx2[i]]][2]
    pmu = data_test[data_test$state==state & data_test$cycle==cycle,c("posteriormean")]
    pstd = data_test[data_test$state==state & data_test$cycle==cycle,c("posteriorstd")]
    vote = data_test[data_test$state==state & data_test$cycle==cycle,c("vote")]
    vote = test_y[[test_idx2[i]]]*100/(sum(test_y[[test_idx2[i]]]))
    party = data_test[data_test$state==state & data_test$cycle==cycle,c("party")]
    candidates = data_test[data_test$state==state & data_test$cycle==cycle,c("candidate")]
    # preds <- sample_posterior(stan_mu[i,], stan_sigma[i,], nc[i], gs=10, ds=1000, fit_params=fit_params)
    preds= c()
    for(j in 1:2){
      tmp = paste('test_y2[',i,',',j,']',sep='')
      pred = fit_params[[tmp]]
      preds = c(preds, pred)
      u=quantile(pred,probs=c(0.975),names = FALSE)
      l=quantile(pred,probs=c(0.025),names = FALSE)
      m = mean(pred)
      s = sd(pred)
      
      CYCLE <- c(CYCLE, cycle)
      STATE <- c(STATE,state)
      CANDIDATE <- c(CANDIDATE,as.character(candidates[j]))
      POSTERIORMEAN <- c(POSTERIORMEAN,pmu[j])
      POSTERIORSTD <- c(POSTERIORSTD,pstd[j])
      PMEAN <- c(PMEAN, m)
      PSTD <- c(PSTD, s)
      VOTE <- c(VOTE, vote[j])
      MEDIAN <- c(MEDIAN, median(pred))
      LOWER95 <- c(LOWER95, l)
      UPPER95 <- c(UPPER95, u)
      NLZ <- c(NLZ, (vote[j]/100-m)^2/2/s^2 + log(s) + log(2*pi)/2)
    }
    preds <- matrix(preds, nrow = 2, byrow = TRUE)
    win_rates = rep(0, 2)
    for(k in 1:ncol(preds)){
      idx = which.max(preds[,k])
      win_rates[idx] = win_rates[idx] + 1
    }
    win_rates = win_rates / sum(win_rates)
    
    WIN <- c(WIN, win_rates)
  }
  
  if(length(test_idx3)){
    for(i in 1:length(test_idx3)) {
      cycle = test_metadata[[test_idx3[i]]][1]
      state = test_metadata[[test_idx3[i]]][2]
      pmu = data_test[data_test$state==state & data_test$cycle==cycle,c("posteriormean")]
      pstd = data_test[data_test$state==state & data_test$cycle==cycle,c("posteriorstd")]
      vote = data_test[data_test$state==state & data_test$cycle==cycle,c("vote")]
      vote= test_y[[test_idx3[i]]]*100/(sum(test_y[[test_idx3[i]]]))
      party = data_test[data_test$state==state & data_test$cycle==cycle,c("party")]
      candidates = data_test[data_test$state==state & data_test$cycle==cycle,c("candidate")]
      # preds <- sample_posterior(stan_mu[i,], stan_sigma[i,], nc[i], gs=10, ds=1000, fit_params=fit_params)
      preds= c()
      for(j in 1:3){
        tmp = paste('test_y3[',i,',',j,']',sep='')
        pred = fit_params[[tmp]]
        preds = c(preds, pred)
        u=quantile(pred,probs=c(0.975),names = FALSE)
        l=quantile(pred,probs=c(0.025),names = FALSE)
        m = mean(pred)
        s = sd(pred)
        CYCLE <- c(CYCLE, cycle)
        STATE <- c(STATE,state)
        CANDIDATE <- c(CANDIDATE,as.character(candidates[j]))
        POSTERIORMEAN <- c(POSTERIORMEAN,pmu[j])
        POSTERIORSTD <- c(POSTERIORSTD,pstd[j])
        PMEAN <- c(PMEAN, m)
        PSTD <- c(PSTD, s)
        VOTE <- c(VOTE, vote[j])
        MEDIAN <- c(MEDIAN, median(pred))
        LOWER95 <- c(LOWER95, l)
        UPPER95 <- c(UPPER95, u)
        NLZ <- c(NLZ, (vote[j]/100-m)^2/2/s^2 + log(s) + log(2*pi)/2)
      }
      preds <- matrix(preds, nrow = 3, byrow = TRUE)
      win_rates = rep(0, 3)
      for(k in 1:ncol(preds)){
        idx = which.max(preds[,k])
        win_rates[idx] = win_rates[idx] + 1
      }
      win_rates = win_rates / sum(win_rates)
      WIN <- c(WIN, win_rates)
    }
  }
  
  if (length(test_idx4)){
    for(i in 1:length(test_idx4)) {
      cycle = test_metadata[[test_idx4[i]]][1]
      state = test_metadata[[test_idx4[i]]][2]
      pmu = data_test[data_test$state==state & data_test$cycle==cycle ,c("posteriormean")]
      pstd = data_test[data_test$state==state & data_test$cycle==cycle ,c("posteriorstd")]
      vote = data_test[data_test$state==state & data_test$cycle==cycle,c("vote")]
      vote= test_y[[test_idx4[i]]]*100/(sum(test_y[[test_idx4[i]]]))
      party = data_test[data_test$state==state & data_test$cycle==cycle,c("party")]
      candidates = data_test[data_test$state==state & data_test$cycle==cycle,c("candidate")]
      # preds <- sample_posterior(stan_mu[i,], stan_sigma[i,], nc[i], gs=10, ds=1000, fit_params=fit_params)
      preds= c()
      for(j in 1:4){
        tmp = paste('test_y4[',i,',',j,']',sep='')
        pred = fit_params[[tmp]]
        preds = c(preds, pred)
        u=quantile(pred,probs=c(0.975),names = FALSE)
        l=quantile(pred,probs=c(0.025),names = FALSE)
        m = mean(pred)
        s = sd(pred)
        CYCLE <- c(CYCLE, cycle)
        STATE <- c(STATE, state)
        CANDIDATE <- c(CANDIDATE,as.character(candidates[j]))
        POSTERIORMEAN <- c(POSTERIORMEAN,pmu[j])
        POSTERIORSTD <- c(POSTERIORSTD,pstd[j])
        PMEAN <- c(PMEAN, m)
        PSTD <- c(PSTD, s)
        VOTE <- c(VOTE, vote[j])
        MEDIAN <- c(MEDIAN, median(pred))
        LOWER95 <- c(LOWER95, l)
        UPPER95 <- c(UPPER95, u)
        NLZ <- c(NLZ, (vote[j]/100-m)^2/2/s^2 + log(s) + log(2*pi)/2)
      }
      preds <- matrix(preds, nrow = 4, byrow = TRUE)
      win_rates = rep(0, 4)
      for(k in 1:ncol(preds)){
        idx = which.max(preds[,k])
        win_rates[idx] = win_rates[idx] + 1
      }
      win_rates = win_rates / sum(win_rates)
      WIN <- c(WIN, win_rates)
    }
  }
  
  # write results to csv
  result <- data.frame(CYCLE,
                       STATE,
                       CANDIDATE,
                       POSTERIORMEAN,
                       POSTERIORSTD,
                       LOWER95,
                       UPPER95,
                       MEDIAN,
                       WIN)
  
  names(result) <- tolower(names(result))

  averaged_nlZs = c(averaged_nlZs, mean(NLZ))
}

write.csv(averaged_nlZs,output_file)

cat(averaged_nlZs)
