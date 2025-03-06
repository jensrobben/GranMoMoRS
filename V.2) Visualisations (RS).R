##### 0) Preliminary settings  ----

# Clear environment
rm(list = ls())
gc()

# Required packages
packages <- c('dplyr','sf','ncdf4','lubridate','tidyverse','tibble',
              'Polychrome','tidyr','ISOweek', 'ggplot2','RColorBrewer',
              'data.table', 'viridis', 'imputeTS', 'imputeTS',
              'parallel','doParallel', 'spdep', 'tidyr','abind',
              'mvtnorm','randtoolbox','MASS','snow','doSNOW',
              'progress')
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
suppressMessages(sapply(packages, require, character.only = TRUE))

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

# All times and ISO week/year information
t_min <- as.Date('2013-01-01', format = '%Y-%m-%d') # Start of ISO week 1 in 2013
t_max <- as.Date('2024-06-30', format = '%Y-%m-%d') # End of ISO week 26 in 2024

time_frame <- seq(t_min, t_max, by="days")
df_itime   <- data.frame('Date' = time_frame, 'ISOWeek' = isoweek(time_frame),
                         'ISOYear' = isoyear(time_frame))

# Colour scales
ablue  <- '#56B4E9'
agreen <- '#A3D596'
ared   <- '#F8766D'
apurple<- '#B1516D'

##### 1) Parameter estimates ----

# Optimal precision parameter
tau <- 10

# Load output RS model
load(paste0("Results/EMresult_tau", tau, ".Rdata"))

# Load Fisher information matrix
FI <- readRDS(file = 'Results/FI.rds')

# Covariance matrices
B <- 25000
FI.alpha <- Reduce('+',lapply(FI, function(x) x$alpha))/B
FI.beta  <- Reduce('+',lapply(FI, function(x) x$beta))/B

COV.alpha <- ginv(FI.alpha)
COV.beta  <- ginv(FI.beta)

SD.alpha <- sqrt(diag(COV.alpha))
SD.beta  <- sqrt(diag(COV.beta))


# Object
obj[[length(obj)]] <- NULL
iteropt <- length(obj)

# Iteration results
alpha1.iter <- do.call('rbind', lapply(obj, function(x) x$alpha.1))
alpha2.iter <- do.call('rbind', lapply(obj, function(x) x$alpha.2))
beta0.iter  <- do.call('rbind', lapply(obj, function(x) x$beta.0))
beta11.iter <- do.call('rbind', lapply(obj, function(x) x$beta.11))
beta22.iter <- do.call('rbind', lapply(obj, function(x) x$beta.22))
gamma.iter  <- do.call('rbind', lapply(obj, function(x) x$gamma))
U.iter      <- do.call('rbind', lapply(obj, function(x) x$U)) 

theta.iter <- cbind(alpha1.iter, alpha2.iter, beta0.iter, 
                    beta11.iter, beta22.iter, gamma.iter)

# Length vector
d <- ncol(alpha1.iter) + ncol(alpha2.iter) + ncol(beta0.iter) + 
  ncol(beta11.iter) + ncol(beta22.iter) + ncol(gamma.iter) - 2

# Optimal parameter vector
alpha1.opt <- alpha1.iter[iteropt,]
alpha2.opt <- alpha2.iter[iteropt,]
beta0.opt  <- beta0.iter[iteropt,]
beta01.opt <- beta0.opt[1:ncol(xx01)]
beta02.opt <- beta0.opt[-c(1:ncol(xx01))]
beta11.opt <- beta11.iter[iteropt,]
beta22.opt <- beta22.iter[iteropt,]
gamma.opt  <- gamma.iter[iteropt,1] # g2 set at zero (g3 st g1 + g3 = 0)
U.opt      <- U.iter[iteropt,]

alpha.opt  <- c(alpha1.opt, alpha2.opt)
beta.opt   <- c(beta0.opt, beta11.opt, beta22.opt)
theta.opt  <- unname(c(alpha1.opt, alpha2.opt, beta0.opt, beta11.opt,
                       beta22.opt))

##### 2) Visualizations -----
### 2.1) In-sample fit based on BE RS path ----
plot_isf <- function(r, a, agegr){
  df.dtxri <- data.frame('Date' = unique(df$Date),
                         'Deaths' = as.vector(dt[,r]),
                         'bDeaths' = as.vector(bt[,r]),
                         'State' = apply(ms.prob[,c(r,r+R,r+2*R)],1, which.max) - 1,
                         'Age' = rep(unique(df$Age0), each = nT),
                         'Region' = regions[r]) %>%
    dplyr::filter(Age == unique(df$Age0)[a]) %>%
    dplyr::mutate('RSDeaths' = (State == 0) * bDeaths + 
                    (State == 1) * bDeaths * exp(x1[((a-1)*nT+1):(a*nT),,r] %*% alpha1.opt) + 
                    (State == 2) * bDeaths * exp(x2[((a-1)*nT+1):(a*nT),,r] %*% alpha2.opt))
  
  val <- 0.04 * diff(range(df.dtxri$Deaths))
  
  p <- ggplot(df.dtxri) + 
    theme_bw(base_size = 15) + ylab('Death counts') + 
    ggtitle(paste0(shapef$NUTS_NAME[shapef$NUTS_ID == regions[r]], ', Age group ', agegr)) + 
    scale_fill_manual(values = brewer.pal(3,'Dark2'), name = 'State') +
    geom_line(aes(x = Date, y = Deaths, color = 'Obs')) +
    geom_line(aes(x = Date, y = RSDeaths, color = 'RS'), alpha = 0.8, linewidth = 1) + 
    geom_line(aes(x = Date, y = bDeaths, color = 'Base'),linewidth = 1) + 
    geom_rect(aes(xmin = Date - 3.5, xmax = Date + 3.5, ymin = min(Deaths) - val, 
                  ymax = min(Deaths) + val, fill = as.factor(State)), alpha = 0.4) + 
    theme(legend.position = 'bottom') + 
    scale_color_manual(values = c('Obs' = 'gray70', 'RS' = ablue, 'Base' = ared),
                       name = '', 
                       labels = c('Obs' = 'Observed',
                                  'Base' = 'Baseline',
                                  'RS' = 'Regime-switching (BE)'),
                       breaks = c('Obs', 'Base', 'RS')) + 
    coord_cartesian(ylim = c(min(df.dtxri$Deaths, df.dtxri$RSDeaths), 
                             max(c(df.dtxri$Deaths, df.dtxri$RSDeaths))),
                    xlim = c(min(df$Date), max(df$Date) -50)) + 
    theme(plot.title = element_text(size=13.5),
          legend.box = 'vertical')
  
  p
}

# Provence-Alpes-Cote d'azur
p1 <- plot_isf(12, 3, '75-79')
p2 <- plot_isf(12, 6, '90+')
p3 <- plot_isf(21, 3, '75-79')
p4 <- plot_isf(21, 6, '90+')

p2.1 <- ggpubr::ggarrange(plotlist = list(p1,p2,p3, p4),
                          common.legend = TRUE, align = 'hv', 
                          nrow = 2, ncol = 2, legend = 'bottom')
p2.1

### 2.2) Mode spatial effect vector ----
spatU <- data.frame('Region' = colnames(adjmat), 'U' = U.opt)
dfspat <- shapef %>% 
  dplyr::left_join(spatU, by = c('NUTS_ID' = 'Region'))


p2.2 <- ggplot() + 
  theme_bw(base_size = 15) + 
  geom_sf(data = dfspat$geometry, aes(fill = dfspat$U), col = 'gray70') + 
  scale_fill_viridis(name = 'U ', option = 'rocket', direction = -1,
                     na.value = 'white') + 
  ggtitle(bquote('Mode spatial effect U*')) + 
  theme(legend.position = 'bottom',
        legend.title = element_text(margin = margin(t = -1.5, unit = "lines")),
        legend.key.size = unit(1, "cm"))
p2.2

### 2.3) Parameter estimates (alpha) ----
df.alpha1 <- data.frame('a1' = alpha1.opt, 
                        'Age' = factor(c( rep(c('Y65-74','Y75-84','Y_GE85'),
                                              times = 3)),
                                       levels = c('Y65-74','Y75-84','Y_GE85')),
                        'Var' = factor(c(rep(c('HI','HI1','HI2','TA','TA1', 'TA2'), 
                                             each = 3)),
                                       levels = c('TA','TA1','TA2','HI','HI1','HI2')),
                        'SE' = sqrt(diag(COV.alpha))[1:length(alpha1.opt)])

df.alpha2 <- data.frame('a2'  = alpha2.opt, 
                        'Age' = factor(c(rep(c('Y65-74','Y75-84','Y_GE85'), 
                                             times = 6)),
                                       levels = c('Y65-74','Y75-84','Y_GE85')),
                        'Var' = factor(c(rep(c('CI01','CI23','IA01',
                                               'IA23', 'HA01','HA23'),
                                             each = 3)),
                                       levels = c('IA01', 'IA23', 'CI01', 
                                                  'CI23', 'HA01', 'HA23')),
                        'SE' = sqrt(diag(COV.alpha))[-c(1:length(alpha1.opt))])


p.alpha1 <- ggplot(df.alpha1, aes(x = Var, y = a1, group = Age, color = Age)) + 
  theme_bw(base_size = 15) + ylab(bquote(alpha['1'])) + xlab('') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black', size = 0.8) + 
  geom_errorbar(aes(ymin = a1 - 1.96*SE, ymax = a1 + 1.96*SE),
                position = position_dodge(0.55), width = 0.5, 
                size = 1.1) + 
  geom_pointrange(aes(ymin = a1 - 1.96*SE, ymax = a1 + 1.96*SE),
                  position = position_dodge(0.55)) +
  scale_color_manual(name = 'Age', 
                     values = c('Y65-74' = ablue, 'Y75-84' = ared, 'Y_GE85' = agreen),
                     labels = c('Y65-74' = '65-74', 'Y75-84' = '75-84', 
                                'Y_GE85' = '85+')) + 
  scale_x_discrete(labels = c('HI' = bquote('HI'['t']),
                              'HI1' = bquote('HI'['t-1']),
                              'HI2' = bquote('HI'['t-2']),
                              'TA' = bquote('TA'['t']), 
                              'TA1' = bquote('TA'['t-1']), 
                              'TA2' = bquote('TA'['t-2']))) + 
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle = 30, vjust = 0.75)) + 
  ggtitle('Environmental shock state')

p.alpha2 <- ggplot(df.alpha2, aes(x = Var, y = a2, group = Age, color = Age)) + 
  theme_bw(base_size = 15) + ylab(bquote(alpha['2'])) + xlab('') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black', size = 0.8) + 
  geom_errorbar(aes(ymin = a2 - 1.96*SE, ymax = a2 + 1.96*SE),
                position = position_dodge(0.55), width = 0.5, 
                size = 1.1) + 
  geom_pointrange(aes(ymin = a2 - 1.96*SE, ymax = a2 + 1.96*SE),
                  position = position_dodge(0.55)) +
  scale_color_manual(name = 'Age', 
                     values = c('Y65-74' = ablue, 'Y75-84' = ared, 
                                'Y_GE85' = agreen),
                     labels = c('Y65-74' = '65-74', 'Y75-84' = '75-84', 
                                'Y_GE85' = '85+')) + 
  scale_x_discrete(labels = c('HA01'  = bquote('HA'['t,t-1']), 
                              'HA23'  = bquote('HA'['t-2,t-3']), 
                              'CI01' = bquote('CI'['t,t-1']),
                              'CI23'  = bquote('CI'['t-2,t-3']),
                              'IA01'  = bquote('IA'['t,t-1']),
                              'IA23'  = bquote('IA'['t-2,t-3']))) + 
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle = 30, vjust = 0.75)) + 
  ggtitle('Respiratory shock state')


p2.3 <- ggpubr::ggarrange(plotlist = list(p.alpha1, p.alpha2), 
                          widths = c(0.4,0.6), align = 'hv',
                          common.legend = TRUE,
                          legend = 'bottom')
p2.3

### 2.4) Parameter estimates (beta) ----
df.b01 <- data.frame('Val' = beta01.opt, 
                     'Var' = factor(c('INT', rep(c('HI'))),
                                    levels = c('INT','HI')),
                     'SE' = sqrt(diag(COV.beta))[1:length(beta01.opt)])
df.b02 <- data.frame('Val' = beta02.opt, 
                     'Var' = factor(c('INT', c('IA01','HA01')),
                                    levels = c('INT','IA01', 'HA01')),
                     'SE' = sqrt(diag(COV.beta))[(length(beta01.opt)+1):
                                                   (length(beta02.opt) + 
                                                      length(beta01.opt))])
df.b11 <- data.frame('Val' = beta11.opt, 
                     'Var' = factor(c('INT', rep(c('HI','HI1', 'HI2'))),
                                    levels = c('INT','HI','HI1','HI2')),
                     'SE' = sqrt(diag(COV.beta))[(length(beta0.opt)+1):
                                                   (length(beta0.opt)+
                                                      length(beta11.opt))])
df.b22 <- data.frame('Val' = beta22.opt, 
                     'Var' = factor(c('INT', c('IA01','HA01', 'IA23','HA23')),
                                    levels = c('INT','IA01','IA23', 'HA01', 'HA23')),
                     'SE' = sqrt(diag(COV.beta))[-c(1:length(c(beta0.opt,beta11.opt)))])


p.b01 <- ggplot(df.b01, aes(x = Var, y = Val)) + 
  theme_bw(base_size = 15) + ylab(bquote(beta['01'])) + xlab('') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black', size = 0.8) + 
  geom_errorbar(aes(ymin = Val - 1.96*SE, ymax = Val + 1.96*SE),
                position = position_dodge(0.55), width = 0.25, 
                size = 1.1, col = ared) + 
  geom_pointrange(aes(ymin = Val - 1.96*SE, ymax = Val + 1.96*SE),
                  position = position_dodge(0.55), col = ared) +
  scale_x_discrete(labels = c('INT' = 'INT', 
                              'HI' = bquote('HI'['t']), 
                              'HI1' = bquote('HI'['t-1']), 
                              'HI2' = bquote('HI'['t-2']))) +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle = 30, vjust = 0.75)) + 
  ggtitle(expression('Baseline' %->% 'Env. shock state'))

p.b02 <- ggplot(df.b02, aes(x = Var, y = Val)) + 
  theme_bw(base_size = 15) + ylab(bquote(beta['02'])) + xlab('') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black', size = 0.8) + 
  geom_errorbar(aes(ymin = Val - 1.96*SE, ymax = Val + 1.96*SE),
                position = position_dodge(0.55), width = 0.25, 
                size = 1.1, col = ared) + 
  geom_pointrange(aes(ymin = Val - 1.96*SE, ymax = Val + 1.96*SE),
                  position = position_dodge(0.55), col = ared) +
  scale_x_discrete(labels = c('INT' = 'INT', 
                              'HA01' = bquote('HA'['t,t-1']), 
                              'HA23' = bquote('HA'['t-2,t-3']), 
                              'IA01' = bquote('IA'['t,t-1']),
                              'IA23' = bquote('IA'['t-2,t-3']))) + 
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle = 30, vjust = 0.75)) + 
  ggtitle(expression('Baseline' %->% 'Respiratory shock state'))

p.b11 <- ggplot(df.b11, aes(x = Var, y = Val)) + 
  theme_bw(base_size = 15) + ylab(bquote(beta['11'])) + xlab('') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black', size = 0.8) + 
  geom_errorbar(aes(ymin = Val - 1.96*SE, ymax = Val + 1.96*SE),
                position = position_dodge(0.55), width = 0.25, 
                size = 1.1, col = ared) + 
  geom_pointrange(aes(ymin = Val - 1.96*SE, ymax = Val + 1.96*SE),
                  position = position_dodge(0.55), col = ared) +
  scale_x_discrete(labels = c('INT' = 'INT', 
                              'HI' = bquote('HI'['t']), 
                              'HI1' = bquote('HI'['t-1']), 
                              'HI2' = bquote('HI'['t-2']))) +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle = 30, vjust = 0.75)) + 
  ggtitle(expression('Env. shock state' %->% 'Env. shock state'))

p.b22 <- ggplot(df.b22, aes(x = Var, y = Val)) + 
  theme_bw(base_size = 15) + ylab(bquote(beta['22'])) + xlab('') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black', size = 0.8) + 
  geom_errorbar(aes(ymin = Val - 1.96*SE, ymax = Val + 1.96*SE),
                position = position_dodge(0.55), width = 0.25, 
                size = 1.1, col = ared) + 
  geom_pointrange(aes(ymin = Val - 1.96*SE, ymax = Val + 1.96*SE),
                  position = position_dodge(0.55), col = ared) +
  scale_x_discrete(labels = c('INT' = 'INT', 
                              'HA01' = bquote('HA'['t,t-1']), 
                              'HA23' = bquote('HA'['t-2,t-3']), 
                              'IA01' = bquote('IA'['t,t-1']),
                              'IA23' = bquote('IA'['t-2,t-3']))) + 
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle = 30, vjust = 0.75)) + 
  ggtitle(expression('Respiratory shock state' %->% 'Respiratory shock state'))

p2.4 <- ggpubr::ggarrange(plotlist = list(p.b01, p.b02, p.b11, p.b22), 
                          widths = c(0.4,0.6), align = 'hv',
                          nrow = 2, ncol = 2)
p2.4

### 2.5) Uncertainty in in-sample fit ----

# Read results
df.p5 <- readRDS('Results/UncQ1.rds')
df.p6 <- readRDS('Results/UncQ2.rds')
df.p7 <- readRDS('Results/UncQ3.rds')

# Construct data frame to plot
df.p5 <- df.p5 %>%
  dplyr::select(c('Date','Age', 'Region', 'Deaths', '2.5%', '50%', '97.5%'))
df.p6 <- df.p6 %>%
  dplyr::select(c('Date','Age', 'Region', 'Deaths', '2.5%', '50%', '97.5%'))
df.p7 <- df.p7 %>%
  dplyr::select(c('Date','Age', 'Region', 'Deaths', '2.5%', '50%', '97.5%'))
colnames(df.p5) <- c('Date','Age', 'Region', 'Deaths', 'Q1.M1', 'Q2.M1', 'Q3.M1')
colnames(df.p6) <- c('Date','Age', 'Region', 'Deaths', 'Q1.M2', 'Q2.M2', 'Q3.M2')
colnames(df.p7) <- c('Date','Age', 'Region', 'Deaths', 'Q1.M3', 'Q2.M3', 'Q3.M3')


plot_unc <- function(region, age, agegr){
  dfplotf <- df.p5 %>% 
    dplyr::left_join(df.p6, by = c('Date', 'Age', 'Region', 'Deaths')) %>%
    dplyr::left_join(df.p7, by = c('Date', 'Age', 'Region', 'Deaths')) %>%
    dplyr::filter(Age == age, Region == region)
  
  p <- ggplot(dfplotf) + 
    theme_bw(base_size = 15) + ylab('Death counts') + xlab('Date')  +
    ggtitle(paste0(shapef$NUTS_NAME[shapef$NUTS_ID == region], ', Age group ', agegr)) +
    geom_line(aes(x = Date, y = Deaths), col = 'gray50') + 
    geom_ribbon(aes(x = Date, ymin = Q1.M1, ymax = Q3.M1, fill = 'M1'), alpha = 0.35) + 
    geom_ribbon(aes(x = Date, ymin = Q1.M2, ymax = Q1.M1, fill = 'M2'), alpha = 0.5) + 
    geom_ribbon(aes(x = Date, ymin = Q3.M1, ymax = Q3.M2, fill = 'M2'), alpha = 0.5) + 
    geom_ribbon(aes(x = Date, ymin = Q1.M2, ymax = Q1.M3, fill = 'M3'), alpha = 0.25) + 
    geom_ribbon(aes(x = Date, ymin = Q3.M2, ymax = Q3.M3, fill = 'M3'), alpha = 0.25) + 
    geom_line(aes(x = Date, y = Q3.M3), col = ablue) +
    geom_line(aes(x = Date, y = Q1.M3), col = ablue) +
    scale_fill_manual(values = c('M1' = ared, 'M2' = agreen, 'M3' = ablue),
                      name = 'Uncertainty',
                      labels = c('M1' = 'State',
                                 'M2' = 'State + Parameter + Spatial',
                                 'M3' = 'State + Parameter + Spatial + Poisson')) + 
    theme(legend.position = 'bottom',
          plot.title = element_text(size = 14)) + 
    coord_cartesian(xlim = c(min(dfplotf$Date),max(dfplotf$Date)-50))
  
  p
}

q1 <- plot_unc('FRG0', 'Y75-79', '75-79')
q2 <- plot_unc('FRG0', 'Y_GE90', '90+')
q3 <- plot_unc('FRL0', 'Y75-79', '75-79')
q4 <- plot_unc('FRL0', 'Y_GE90', '90+')

p2.5 <- ggpubr::ggarrange(plotlist = list(q1, q2, q3, q4),
                          common.legend = TRUE, nrow = 2,
                          ncol = 2, align = 'hv', legend = 'bottom') 
p2.5
