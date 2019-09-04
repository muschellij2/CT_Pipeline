rm(list=ls())
library(reshape2)
library(lubridate)
x = read.csv("test.csv", stringsAsFactors = FALSE)
df = x[ , c(grep("time", colnames(x), value=TRUE), "Runner") ]
df$total.time = NULL
head(df)
df$Runner = factor(df$Runner, levels = unique(df$Runner))

long = melt(df, id.vars = "Runner")
long$variable = as.numeric(substr(long$variable, 2, 2))

long$value = hms(long$value)
start = ymd_hms("2015-09-11 9:30:00")
long$tod = start
# vals = long$value[1]
# for (ival in seq(2, nrow(long))){
#     vals = c(vals, long$value[ival]+ vals[ival-1])
# }
hr_to_sec = function(x){
    hr = as.numeric(hour(x)) * 60 * 60
    min = as.numeric(minute(x)) * 60
    sec = as.numeric(second(x))
    hr + min + sec
}
long$v = cumsum(c(0, hr_to_sec(long$value[1:(nrow(long)-1)])))
long$start_time = long$tod + long$v

long$tod = long$v = NULL

df = reshape(long, 
    direction = "wide", 
    idvar = "Runner",
    timevar = "variable")

write.csv(df, 
    file = "Running_Starts.csv", row.names = FALSE)