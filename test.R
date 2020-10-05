library(dplyr)
# define horizons in the loyo process
horizons = c('0',
             '7',
             '14',
             '21',
             '28',
             '42',
             '56')


# the optimal index of hyperparameters in the loyo process
# optimal index can be obtained with loocv_nlZs.R
best_gp_idx = read.csv(paste("results/", "GP", "_opthyp.csv", sep=''))
best_gp_idx = best_gp_idx$opt_idx

best_lm_idx = read.csv(paste("results/", "LM", "_opthyp.csv", sep=''))
best_lm_idx = best_lm_idx$opt_idx

PVALUES = c()

gp_avg_nlz = c()
lm_avg_nlz = c()

for (a in 1:length(horizons)) {
  
  # load the prior files
  gp_nlz = paste('results/stan_NLZ', 'GP', '_' , test_year, 'day', horizons[a], '_', best_gp_idx[a] ,'.csv',sep='')
  lm_nlz = paste('results/stan_NLZ', 'LM', '_' , test_year, 'day', horizons[a], '_', best_lm_idx[a] ,'.csv',sep='')
  gp_data <- read.csv(gp_nlz)
  lm_data <- read.csv(lm_nlz)
  
  gp_data = gp_data$nlz
  lm_data = lm_data$nlz
  
  print(strtoi(horizons[a]))
  TEST = t.test(gp_data, lm_data, paired = TRUE)
  PVALUES = c(PVALUES, TEST$p.value)
  print(TEST)
  
  gp_avg_nlz = c(gp_avg_nlz, (gp_data))
  lm_avg_nlz = c(lm_avg_nlz, (lm_data))
}

# t.test(gp_avg_nlz, lm_avg_nlz, paired = TRUE)

results <- data.frame(horizons, PVALUES)
write.csv(results, "results/pvalues.csv")


generate_results <- function(TYPE, best_idx, test_year, horizons){
  RMSE = c() 
  ACCURACY = c()
  ENTROPY = c()
  RATIO = c()
  NLZ = c()
  for (a in length(horizons):1) {
    # load result files
    data = paste('results/stan_LOO', TYPE, '_' , test_year, 'day', horizons[a], '_', best_idx[a] ,'.csv',sep='')
    data <- read.csv(data)
    
    RMSE = c(RMSE, mean(data$rmse))
    
    correct_prediction = 0
    for(state in unique(data$state)){
      tmp = data[data$state==state,c("win","vote")]
      if( which(max(tmp$win)==tmp$win)== which(max(tmp$vote)==tmp$vote)){
        correct_prediction = correct_prediction + 1
      }
    }
    
    ACCURACY = c(ACCURACY, correct_prediction/length(unique(data$state)))
    ENTROPY = c(ENTROPY, mean(log(data$win+1e-10)*data$winners))
    RATIO = c(RATIO, sum((data$vote>data$lower95) & (data$vote<data$upper95))/nrow(data))
    
    nlz = paste('results/stan_NLZ', TYPE, '_' , test_year, 'day', horizons[a], '_', best_idx[a] ,'.csv',sep='')
    nlz <- read.csv(nlz)
    nlz = nlz$nlz
    NLZ = c(NLZ, mean(nlz))
  }
  
  results <- data.frame(RMSE,
                        ACCURACY,
                        ENTROPY,
                        RATIO,
                        NLZ)
  
  results = results %>% mutate_if(is.numeric, round, digits = 4)
  write.csv(t(results), paste("results/", TYPE,test_year,"results.csv",sep=""))
}

generate_results('GP',best_gp_idx, 2018, horizons)

generate_results('LM',best_lm_idx, 2018, horizons)

