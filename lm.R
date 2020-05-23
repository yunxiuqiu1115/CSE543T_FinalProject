library(dplyr)
library(dummies)

horizons = c(120,90,42,28,14,0)
df <- read.csv("data/CNNData1992to2018idx.csv")
c2v = df %>% group_by(Candidateidentifier) %>% summarise_at(vars(Percentage_of_Vote_won_x), funs(mean(., na.rm=TRUE)))

# df$Republican = df$Republican*2-1
df$polls <- df$numberSupport/df$samplesize
colnames(df)[colnames(df) == 'Percentage_of_Vote_won_x'] = 'vote'
df$pollsteridx = as.character(df$pollsteridx)
# df <- dummy.data.frame(data=df, names='pollsteridx')

for (horizon in horizons){
  model <- lm(polls ~ Candidateidentifier 
              + Candidateidentifier:daysLeft - 1, data = df[df$daysLeft<=-horizon,])
  
  result <- summary(model)
  
  output = data.frame(matrix(ncol = 9, nrow = 0))
  colnames(output) = c('cycle','state','candidate','posteriormean','posteriorstd', 'vote','pvi', 'party','experienced')
  
  for (i in 1:nrow(c2v)) {
    c = as.character(c2v[i,1]$Candidateidentifier)
    # a = result$coefficients[paste('Candidateidentifier',c,':daysLeft', sep=''),1]
    skip_to_next <- FALSE
    tryCatch(result$coefficients[paste('Candidateidentifier',c, sep=''),1], error = function(e) { skip_to_next <<- TRUE})
    if(skip_to_next) { next }    
    
    b = result$coefficients[paste('Candidateidentifier',c, sep=''),1]
    std_b = result$coefficients[paste('Candidateidentifier',c, sep=''),2]
    # x = df[df$Candidateidentifier==c,]$daysLeft
    # y = df[df$Candidateidentifier==c,]$polls
    # t = min(x):0
    # z = a*t+b
    # plot(t,z, type="l")
    # points(x,y)
    tmp = df[df$Candidateidentifier==c,c('cycle','state','vote','pvi','experienced','Republican')][1,]
    tmp['posteriormean'] = b
    tmp['posteriorstd'] = std_b
    tmp['Republican'] = 2*tmp['Republican']-1
    tmp['candidate'] = substring(c, 7)
    colnames(tmp)[colnames(tmp) == 'Republican'] = 'party'
    output = rbind(output, tmp)
  }
  
  write.csv(output, paste('results/lm1992-2018all', horizon, '.csv', sep=''), row.names=FALSE)
}


