# define horizons in the loyo process
horizons = c('0',
               '7',
               '14',
               '28',
               '42',
               '90',
               '120')

# number of candidate points for loyo
search_size = 100

# all years in loyo process
cv_years = c(1992,1994,1996,1998,2000,2002,2004,2006,2008,2010,2012,2014,2016)

# model type
TYPES = c('GP' ,'LM')

for (TYPE in TYPES) {
  for (a in 1:length(horizons)) {
    averaged_LL = rep(0, search_size)
    for (b in 1:search_size){
      year_nlZs = c()
      for (cv_year in cv_years) {
        input_file = paste('nlZs/', TYPE, '_' , cv_year, 'day', horizons[a], '_', b,'.csv',sep='')
        nlZ <- read.csv(input_file)
        nlZ = nlZ$x
        year_nlZs = c(year_nlZs, nlZ)
      }
      averaged_LL[b] = -mean(year_nlZs)
    }
    # obtain the maximal index of the averaged ll array
    print(paste(TYPE, ' day ', horizons[a] ,sep=''))
    print(which(averaged_LL==sort(max(averaged_LL))))
  }
}
