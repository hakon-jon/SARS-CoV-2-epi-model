# Viral concentration regression 
sim_data <- data.frame(Date=seq(from=1,to=60),
                       CT=-1*seq(from=1,to=60)+rnorm(60))
fit <- lm(CT~Date,data=sim_data)

# Number of Icelandic private mutations 
sim_ice_data <- data.frame(Date=seq(from=1,to=60),
                           Muts=sapply(seq(from=1,to=60),function(x){rpois(lambda=x/14,n=1)}))
fit_private <- glm(Muts~Date,family=poisson(link='identity'),data=sim_ice_data,start=c(Intercept=0,Date=1/14))

