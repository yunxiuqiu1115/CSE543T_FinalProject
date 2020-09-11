# setwd('/Users/yahoo/Documents/WashU/CSE515T/Code/Gaussian Process/')

library(rstan)
# library(MCMCpack)
# library(bayesplot)
# library(dplyr)
# library(ggridges)
# library(ggplot2)
# library(grid)
# library(pBrackets)


input_strs = c('0',
               '7',
               '14',
               '28',
               '42',
               '90',
               '120')

search_size = 100

cv_years = c(1992,1994,1996,1998,2000,2002,2004,2006,2008,2010,2012,2014,2016)

TYPE = 'GP'
# TYPE = 'LM'

# f
best_cv_idx = c(89,89, 66,12,12,36,36)

# fit_objs = c()

LLs = matrix(data=0, nrow=length(input_strs), ncol=search_size)

search_size = 1
input_strs = c('42')
best_cv_idx = c(12)
cv_years = c(2020)

for (a in 1:length(input_strs)) {
  for (b in 1:search_size){
    cv_LL = c()
    
    for (cv_year in cv_years) {
      input_file = paste('results/LOO', TYPE, '_' , cv_year, 'day', input_strs[a], '_', best_cv_idx[a] ,'.csv',sep='')
      output_file = paste('results/stan_LOO', TYPE, '_' , cv_year, 'day', input_strs[a], '_', best_cv_idx[a] ,'.csv',sep='')
      
      # loading data
      data <- read.csv(input_file)
      print(input_file)
      data <- data[data$cycle!=2016 | data$state!='Louisiana' | data$candidate!='Flemsing',]
      data <- data[data$cycle!=2020 | data$state!='Georgia' | data$candidate!='Loeffler',]
      data <- data[data$cycle!=2020 | data$state!='Georgia' | data$candidate!='Tarver',]
      
      # data %>%
      #   group_by(cycle, state) %>%
      #   summarise(count=n()) %>%
      #   filter(count >=4)
      
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
    # for (cycle in c(2016, 2018)){
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
    
    # print(summary(fit,c('alpha','beta','ppb','eb','year_sig'))$summary)
    fit_params <- as.data.frame(fit)
    
    # fit_objs = c(fit_objs, fit)
  
    
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
    
    # for(i in 1:length(idx2)) {
    #   state = metadata[[idx2[i]]]
    #   preds = list()
    #   for(j in 1:2){
    #     tmp = paste('p2[',i,',',j,']',sep='')
    #     pred = fit_params[[tmp]]
    #     preds[[j]] = pred
    #   }
    #   for(j in 1:2){
    #     pred = preds[[j]] / (preds[[1]]+preds[[2]])
    #     u=quantile(pred,probs=c(0.975),names = FALSE)
    #     l=quantile(pred,probs=c(0.025),names = FALSE)
    #     m = mean(pred)
    #     s = sd(pred)
    #     if (stan_y[idx2[i],j]>u & stan_y[idx2[i],j]<l){
    #       Nout = Nout + 1
    #     }
    #   }
    # }
    # 
    # for(i in 1:length(idx3)) {
    #   state = metadata[[idx3[i]]]
    #   preds = list()
    #   for(j in 1:3){
    #     tmp = paste('p3[',i,',',j,']',sep='')
    #     pred = fit_params[[tmp]]
    #     preds[[j]] = pred
    #   }
    #   for(j in 1:3){
    #     pred = preds[[j]] / (preds[[1]]+preds[[2]]+preds[[3]])
    #     u=quantile(pred,probs=c(0.975),names = FALSE)
    #     l=quantile(pred,probs=c(0.025),names = FALSE)
    #     m = mean(pred)
    #     s = sd(pred)
    #     if (stan_y[idx3[i],j]>u & stan_y[idx3[i],j]<l){
    #       Nout = Nout + 1
    #     }
    #   }
    # }
    # 
    # for(i in 1:length(idx4)) {
    #   state = metadata[[idx4[i]]]
    #   preds = list()
    #   for(j in 1:4){
    #     tmp = paste('p4[',i,',',j,']',sep='')
    #     pred = fit_params[[tmp]]
    #     preds[[j]] = pred
    #   }
    #   for(j in 1:4){
    #     pred = preds[[j]] / (preds[[1]]+preds[[2]]+preds[[3]]+preds[[4]])
    #     u=quantile(pred,probs=c(0.975),names = FALSE)
    #     l=quantile(pred,probs=c(0.025),names = FALSE)
    #     m = mean(pred)
    #     s = sd(pred)
    #     if (stan_y[idx4[i],j]>u & stan_y[idx4[i],j]<l){
    #       Nout = Nout + 1
    #     }
    #   }
    # }

    # print(paste("In-sample ratio in 95% :",1-Nout/760))
    
    
    COLNAMES = c('Posterior_Vote','Party','State','Type')
    posteriors = data.frame(matrix(ncol = length(COLNAMES), nrow = 0))
    colnames(posteriors) = COLNAMES

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
        
        one_posterior = data.frame(matrix(ncol = length(COLNAMES), nrow = length(pred)))
        colnames(one_posterior) = COLNAMES
        
        one_posterior$Posterior_Vote = pred*100
        one_posterior$State = state.abb[match(state,state.name)]
        if(party[j]==1){
          one_posterior$Party = 'REP'
          RepWin = mean(one_posterior$Posterior_Vote>50)
        }
        else{
          one_posterior$Party = 'DEM'
          RepWin = mean(one_posterior$Posterior_Vote<50)
        }
        if(RepWin>0.95){
          one_posterior$Type = "Safe R"
        }
        else if (RepWin>0.7) {
          one_posterior$Type = "Likely R"
        }
        else if (RepWin>0.55) {
          one_posterior$Type = "Lean R"
        }
        else if (RepWin>0.45) {
          one_posterior$Type = "Toss-up"
        }
        else if (RepWin>0.3) {
          one_posterior$Type = "Lean D"
        }
        else if (RepWin>0.05) {
          one_posterior$Type = "Likely D"
        }
        else {
          one_posterior$Type = "Safe D"
        }
        posteriors = rbind(posteriors, one_posterior)
        
        # if (test_stan_y[test_idx2[i],j]<=u & test_stan_y[test_idx2[i],j]>=l){
        #   flags[test_idx2[i],j] = 1
        # }
        # else{
        #   Nout_test = Nout_test + 1
        # }
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
      # if (which.max(win_rates)==which.max(vote)){
      #   correct_predictions = correct_predictions + 1
      # }
      # else{
      #   print("Wrong prediction:")
      #   print(test_metadata[[test_idx2[i]]])
      # }
    }
    
    
    # posteriors %>%
    #   ggplot(aes(y = State)) +
    #   geom_density_ridges(
    #     aes(x = Posterior_Vote, fill = paste(State, Party)),
    #     alpha = .8, color = 'grey', from = 0, to = 100
    #   ) +
    #   labs(
    #     x = "Posterior Vote (%)",
    #     y = "Election State",
    #     title = "Forecasting REP vs DEM Candidate Vote in 2020 US Senate elections"
    #   ) +
    #   scale_y_discrete(expand = c(0, 0)) +
    #   scale_x_continuous(expand = c(0, 0)) +
    #   scale_fill_cyclical(
    #     breaks = c("REP", "DEM"),
    #     labels = c(`REP` = "REP", `DEM` = "DEM"),
    #     values = c("#0000ff", "#ff0000","#ff8080", "#8080ff"),
    #     name = "Party", guide = "legend"
    #   ) +
    #   coord_cartesian(clip = "off") +
    #   theme_ridges(grid = FALSE)
    
    LEVELS = posteriors[posteriors$Party=='REP',] %>% 
      group_by(State) %>% 
      mutate(tmp = mean(Posterior_Vote)) %>% 
      select(State, tmp) %>%
      distinct(State, tmp) %>%
      arrange(desc(tmp)) %>%
      select(State)
    
    LEVELS = LEVELS$State
    
    # WINS = posteriors[posteriors$Party=='REP',] %>% 
    #   group_by(State) %>% 
    #   mutate(win = mean(Posterior_Vote>=50)) %>% 
    #   select(State, win) %>%
    #   distinct(State, win) %>%
    #   arrange(desc(win))
    
    posteriors$State <- factor(posteriors$State, levels = LEVELS)
    
    LIKENAMES = c('Safe R', 'Likely R','Lean R', 'Toss-up', 'Lean D','Likely D','Safe D')
    posteriors$Type <- factor(posteriors$Type, levels = LIKENAMES)

    # ggplot(posteriors, aes(x = Posterior_Vote, y = reorder(State, desc(State)), color = Party, fill = Party)) +
    #   geom_density_ridges(alpha=0.8) +
    #   scale_y_discrete(expand = c(0, 0), name = "") +
    #   scale_x_continuous(expand = c(0, 0), breaks = c(0,10,20,30,40,45,50,55,60,70,80,90,100),
    #                      labels = c('0','Safe D','20', 'Likely D','40', 'Lean D', '50', 'Lean R', '60', 'Likely R', '80', 'Safe R', '100'),
    #                      name = "Posterior Vote (%)") +
    #   theme(axis.text.x = element_text(angle=60),
    #         axis.ticks.x = element_line(size = c(.5,0,.5,0,.5,0,.5,0,.5,0,.5,0,.5),
    #                                     color = c("black", NA,"black", NA, "black", NA,"black", NA, "black", NA,"black", NA, "black")),
    #         panel.grid.minor = element_blank(),
    #         panel.grid.major.x = element_line(color = c("gray",
    #                                                     NA,"gray", NA, "gray", NA,"gray", NA, "gray", NA,"gray", NA, "gray"))) +
    #   scale_fill_manual(values = c("#34AAE0", "#E9141D"), labels = c("DEM", "REP")) +
    #   scale_color_manual(values = c("#D3D3D3","#D3D3D3"), guide = "none") +
    #   coord_cartesian(xlim = c(0, 100), clip='on') +
    #   guides(fill = guide_legend(
    #     override.aes = list(
    #       fill = c("#34AAE0", "#E9141D"),
    #       color = NA, point_color = NA)
    #   )
    #   ) +
    #   ggtitle("Posterior densities for vote shares for major-party candidates") +
    #   theme(plot.title = element_text(hjust=0.5),
    #         panel.background = element_rect(fill = 'white', colour = 'white'))
    
    # ggplot(posteriors, aes(x = Posterior_Vote, y = reorder(State, desc(State)), color = Party, fill = Party)) +
    #   geom_density_ridges(alpha=0.6) +
    #   scale_y_discrete(expand = c(0, 0), name = "") +
    #   # facet_wrap(Type ~ ., scale ="free") +
    #   scale_x_continuous(expand = c(0, 0), breaks = c(0,20,40,60,80,100),
    #                      name = "Posterior Vote (%)") +
    #   theme(panel.grid.minor = element_blank(),
    #        panel.grid.major.x = element_line(color = "gray")) +
    #   scale_fill_manual(values = c("blue", "red"), labels = c("DEM", "REP")) +
    #   scale_color_manual(values = c(NA,NA), guide = "none") +
    #   coord_cartesian(xlim = c(0, 100), clip='on') +
    #   guides(fill = guide_legend(
    #     override.aes = list(
    #       fill = c("blue", "red"),
    #       color = NA, point_color = NA)
    #   )
    #   ) +
    #   ggtitle("Posterior predictive density of vote share for major party candidates") +
    #   theme(plot.title = element_text(hjust=0.5),
    #         panel.background = element_rect(fill = 'white', colour = 'white'))


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
          # if (test_stan_y[test_idx3[i],j]<=u & test_stan_y[test_idx3[i],j]>=l){
          #   flags[test_idx3[i],j] = 1
          # }
          # else{
          #   Nout_test = Nout_test + 1
          # }
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
        # if (which.max(win_rates)==which.max(vote)){
        #   correct_predictions = correct_predictions + 1
        # }
        # else{
        #   print("Wrong prediction:")
        #   print(test_metadata[[test_idx3[i]]])
        # }
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
          # if (test_stan_y[test_idx4[i],j]<=u & test_stan_y[test_idx4[i],j]>=l){
          #   flags[test_idx4[i],j] = 1
          # }
          # else{
          #   Nout_test = Nout_test + 1
          # }
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
        # if (which.max(win_rates)==which.max(vote)){
        #   correct_predictions = correct_predictions + 1
        # }
        # else{
        #   print("Wrong prediction:")
        #   print(test_metadata[[test_idx4[i]]])
        # }
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
    
    cv_LL = c(cv_LL, mean(-NLZ))
    
    write.csv(result,output_file)
    
    # print(paste("Correct predictions: ",correct_predictions))
    # 
    # print(paste("Accuracy: ",correct_predictions/length(test_metadata)))
    # 
    # print(paste("Correlation: ",cor(PMEAN, VOTE)))
    # 
    # print(paste("RSME: ",sqrt(mean((PMEAN- VOTE/100)^2))))
    # 
    # print(paste("Ratio in 95% : ",1-Nout_test/73))
    # 
    # print(paste("Predictive averaged nlZ: ",mean(NLZ)))
    # 
    # print(paste("Mean of predictive std: ",mean(PSTD)))
    # 
    # print(paste("Median of predictive std: ",median(PSTD)))
    # 
    # print(paste("Std of predictive std: ",sd(PSTD)))
    
    }
    LLs[a,b] = mean(cv_LL)
  }
}


for (a in 1:length(input_strs)) {
  print(input_strs[a])
  print(which(LLs[a,]==max(LLs[a,])))
}

# save.image(file = paste('models/LOOCV_',TYPE ,'.RData',sep=''))
