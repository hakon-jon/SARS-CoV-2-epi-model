rm(list=ls())
library('deSolve')
library(ggplot2)
library(reshape2)

alpha_per_indi <- (0.197)/3.5e5 
number_of_possible_generations <- 30

trajectory <- function(kappa_rel,kill_ratio){
    # Growth phase
    parameters <- list(alpha=alpha_per_indi,
                       beta=3.5e5*alpha_per_indi/2, # setting the reproduction number to 2
                       kappa=3.5e5*alpha_per_indi/2*kappa_rel,
                       number_of_states=number_of_possible_generations)

    # Containment phase
    second_parameters <- list(alpha=alpha_per_indi*kill_ratio,
                       beta=3.5e5*alpha_per_indi/2,
                       kappa=3.5e5*alpha_per_indi/2*kappa_rel,
                       number_of_states=parameters$number_of_states)

    #
    times <- seq(0,30,by=1)
    times2 <- seq(0,30,by=1)
    state= c(S=0,I1=0)
    state <- rep(0,2+2*parameters$number_of_states)
    names(state) <- c('S',
                      paste('I',seq(from=1,to=parameters$number_of_states),sep=''),
                      paste('P',seq(from=1,to=parameters$number_of_states),sep=''),
                      'R')
    state['S'] <- 3.5e5
    state['I1'] <- 100

    rate_of_state_change <- function(t,state,parameters){
        with(as.list(c(parameters)),{
             ret <- list()
            I = sum(state[grepl('^I',names(state))])
            P = sum(state[grepl('^P',names(state))])
            ret['dS'] = - alpha *state['S']*(I)

            # The infection rate
            ret['dI1'] = -beta*state['I1']
            for (i in seq(from=2,to=number_of_states)){
                ret[paste('dI',i,sep='')] = alpha*state['S']*state[paste('I',i-1,sep='')]-beta*state[paste('I',i,sep='')]
            }

            # The positive rate 
            for (i in seq(from=1,to=number_of_states)){
                ret[paste('dP',i,sep='')] <- beta*state[paste('I',i,sep='')]-kappa*state[paste('P',i,sep='')]
            }

            ret['dR'] = kappa*P
            return(list(c(ret)))
        })
    }

    first_period <- ode(y=state,times=times,func=rate_of_state_change,parms=parameters)
    # the initial values in the second period are the end time point from the growth epoch. 
    second_period <- ode(y=first_period[nrow(first_period),-1],times=times2,func=rate_of_state_change,parms=second_parameters)
    
    last_time <- (first_period[nrow(first_period),'time'])
    te <- transform(first_period,Epoch='Explosive')
    te2 <- transform(second_period,Epoch='Decaying')
    te2$time <- te2$time+last_time
    combined <- rbind(te,te2)
    combined$Kappa <- kappa
    combined$KillRatio <- kill_ratio
    return(combined)
}

# Calculate epidemilogical trajectories for a range of Kappa.
combined <- data.frame()
for (kappa in c(0.1,0.5,1,2,10)){
    for (kill_ratio in c(10**(log10(0.7)*seq(from=0,to=8)),0.1)){
        combined <- rbind(combined,trajectory(kappa,kill_ratio))
    }
}
combined$NotS <- 3.5e5-combined$S
combined$I <- apply(combined[,grepl('^I',colnames(combined))],1,sum)
combined$P <- apply(combined[,grepl('^P',colnames(combined))],1,sum)
combined$IP <- combined$I+combined$P
combined$InfectedAmongPositive <- (combined$I)/(combined$I+combined$P)
combined$Kappa <- factor(combined$Kappa)
#
combined$AccumulationOfMutations  <- 0
combined$NumberOfGenerations  <- 0
for (i in seq(from=1,to=number_of_possible_generations)){
    combined$AccumulationOfMutations <- combined$AccumulationOfMutations+combined[,paste('I',i,sep='')]/combined$I*(i-1)*0.2+combined[,paste('P',i,sep='')]/combined$P*0.2*(i-1)
    combined$NumberOfGenerations <- combined$NumberOfGenerations+combined[,paste('I',i,sep='')]/combined$I*i+combined[,paste('P',i,sep='')]/combined$P*i
}

#

plot_dat <- subset(melt(combined,
                        id=c('Epoch','time','Kappa','KillRatio')),
                   variable %in% c('S','I','R','P','InfectedAmongPositive','NumberOfGenerations') )
plot_dat$variable <- factor(ifelse(plot_dat$variable=='InfectedAmongPositive','I/(I+P)',as.character(plot_dat$variable)),c('S','I','P','R','I/(I+P)','NumberOfGenerations'))
data_for_main_figure  <- subset(plot_dat,KillRatio==0.1)[,c('time','Kappa','variable','value')]

figure_5B <- ggplot(data_for_main_figure,
                    aes(x=time,
                        y=value,
                        col=Kappa))+
             geom_vline(xintercept=30,
                        linetype='dashed')+
             geom_line()+
             facet_wrap(variable~.,
                        scale='free_y',
                        ncol=2)+
             theme_classic()+
             xlab('Days')+
             ylab('Value')

pdf('figure_5B.pdf')
    plot(figure_5B)
dev.off()

reduction <- subset(plot_dat,Kappa==0.5 & KillRatio!=0.1 & variable %in% c('NumberOfGenerations'))

reduction_relative <- merge(reduction,
                  subset(reduction,KillRatio==1),
                  by=c('Epoch','time','Kappa','variable'))


# Numbers for the main text 
# Amount of reduction
amount_of_reduction_needed <- max(subset(reduction_relative,value.x/value.y<0.7 & time<58)$KillRatio.x)
print(100*(1-amount_of_reduction_needed))


reduction<- transform(reduction,Compatible_Observed=KillRatio==amount_of_reduction_needed)
reduction_relative <- transform(reduction_relative,Compatible_Observed=KillRatio.x==amount_of_reduction_needed)


gg_color <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
colors <- gg_color(length(unique(reduction$KillRatio)))
names(colors) <- unique(signif(reduction$KillRatio,2))
colors[unique(reduction$KillRatio)==amount_of_reduction_needed] <- 'black'

data_for_plot <- rbind(transform(reduction,
                                 variable='Number of generations',
                                 y=value,
                                 col=factor(signif(KillRatio,2)))[,c('time','variable','y','col','Compatible_Observed')],
                  transform(reduction_relative,
                            variable='Relative decrease in generations',
                            y=value.x/value.y,
                            col=factor(signif(KillRatio.x,2)))[,c('time','variable','y','col','Compatible_Observed')])


figure_5C <- ggplot()+
             geom_vline(xintercept=30,
                        linetype='dashed')+
             geom_line(data=data_for_plot,
                       aes(x=time,
                           y=y,
                           col=col))+
             geom_line(data=subset(data_for_plot,
                                   Compatible_Observed),
                       aes(x=time,
                           y=y))+
             facet_wrap(variable~.,
                        scale='free_y')+
             theme_classic()+
             theme(legend.position='bottom')+
             scale_color_manual('Relative viral\nreproduction rate',
                                values=colors)+
             xlab('Days')+
             ylab('Value')
pdf('figure_5C.pdf')
    plot(figure_5C)
dev.off()

