THIN = 10

for(year in c(2000)){
  YS = c()
  YDOTS = c()
  YREPDOTS = c()
  YREPS = c()
  STATES = c()
  for(i in 1:length(idx2)) {
    cycle = metadata[[idx2[i]]][1]
    state = metadata[[idx2[i]]][2]
    if (year!=cycle){
      next
    }
    vote = data[data$state==state & data$cycle==cycle,c("vote")]
    vote = y[[idx2[i]]]/(sum(y[[idx2[i]]]))
    for(j in 1:2){
      STATES = c(STATES, state)
      tmp = paste('rep_y2[',i,',',j,']',sep='')
      pred = fit_params[[tmp]]
      YDOTS = c(YDOTS, vote[j])
      YREPDOTS = c(YREPDOTS, mean(pred))
      pred = pred[c(rep(FALSE, THIN), TRUE)]
      YS = c(YS, rep(vote[j], length(pred)))
      YREPS = c(YREPS, pred)
    }
  }
  
  for(i in 1:length(idx3)) {
    cycle = metadata[[idx3[i]]][1]
    state = metadata[[idx3[i]]][2]
    if (year!=cycle){
      next
    }
    vote = data[data$state==state & data$cycle==cycle,c("vote")]
    vote = y[[idx3[i]]]/(sum(y[[idx3[i]]]))
    for(j in 1:3){
      STATES = c(STATES, state)
      tmp = paste('rep_y3[',i,',',j,']',sep='')
      pred = fit_params[[tmp]]
      YDOTS = c(YDOTS, vote[j])
      YREPDOTS = c(YREPDOTS, mean(pred))
      pred = pred[c(rep(FALSE, THIN), TRUE)]
      YS = c(YS, rep(vote[j], length(pred)))
      YREPS = c(YREPS, pred)
    }
  }
  
  for(i in 1:length(idx4)) {
    cycle = metadata[[idx4[i]]][1]
    state = metadata[[idx4[i]]][2]
    if (year!=cycle){
      next
    }
    vote = data[data$state==state & data$cycle==cycle,c("vote")]
    vote = y[[idx4[i]]]/(sum(y[[idx4[i]]]))
    for(j in 1:4){
      STATES = c(STATES, state)
      tmp = paste('rep_y4[',i,',',j,']',sep='')
      pred = fit_params[[tmp]]
      YDOTS = c(YDOTS, vote[j])
      YREPDOTS = c(YREPDOTS, mean(pred))
      pred = pred[c(rep(FALSE, THIN), TRUE)]
      YS = c(YS, rep(vote[j], length(pred)))
      YREPS = c(YREPS, pred)
    }
  }
  # Open a pdf file
  TITLE = paste(TYPE, '_' , year, 'day 0', sep='')
  pdf(paste('plots/', TITLE, '.pdf', sep='')) 
  # 2. Create a plot
  plot(YS, YREPS, main=TITLE, type ='p', col = rgb(red = 0.3, green = 0.3, blue = 0.3),
       xlab="y", ylab="y_rep", pch = 20, cex = 0.2)
  lines(seq(from = 0, to = 1, by = 0.01), seq(from = 0, to = 1, by = 0.01), lty=3)
  points(YDOTS, YREPDOTS, col='red', pch=4)
  legend("topleft",
         c("y_rep vs y","y=x",'mean(y_rep) vs y'),
         fill=c(rgb(red = 0.3, green = 0.3, blue = 0.3),"black","red")
  )
  # 3. Close the file
  dev.off() 
}