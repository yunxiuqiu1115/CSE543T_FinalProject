input_strs = c('0',
               '7',
               '14',
               '28',
               '42',
               '90',
               '120')

search_size = 100

cv_years = c(1992,1994,1996,1998,2000,2002,2004,2006,2008,2010,2012,2014,2016)

TYPES = c('LM')

for (TYPE in TYPES) {
  for (a in 1:length(input_strs)) {
    averaged_LL = rep(0, search_size)
    for (b in 1:search_size){
      year_nlZs = c()
      for (cv_year in cv_years) {
        input_file = paste('nlZs/', TYPE, '_' , cv_year, 'day', input_strs[a], '_', b,'.csv',sep='')
        nlZ <- read.csv(input_file)
        nlZ = nlZ$x
        year_nlZs = c(year_nlZs, nlZ)
      }
      averaged_LL[b] = -mean(year_nlZs)
    }
    print(paste(TYPE, ' day ', input_strs[a] ,sep=''))
    # print(max(averaged_LL))
    # print(sort(averaged_LL,partial=99)[99])
    print(which(averaged_LL==sort(max(averaged_LL))))
  }
}
