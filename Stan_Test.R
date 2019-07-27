init <- function(){
  list(CL = exp(rnorm(1, log(10), 0.2)),
       Q = exp(rnorm(1, log(15), 0.2)),
       V1 = exp(rnorm(1, log(35), 0.2)),
       V2 = exp(rnorm(1, log(100), 0.2)),
       ka = exp(rnorm(1, log(1), 0.2)),
       ke0 = exp(rnorm(1,log(1),0.2)),
       sigma = 0.5,
       sigmaResp = 20)
}


parametersToPlot <- c("CL", "Q", "V1", "V2", "ka", "sigma")

## Additional variables to monitor
otherRVs <- c("cObsPred")

parameters <- c(parametersToPlot, otherRVs)
parametersToPlot <- c("lp__", parametersToPlot)

nChains <- 4
nPost <- 1000 ## Number of post-warm-up samples per chain after thinning
nBurn <- 1000 ## Number of warm-up samples per chain after thinning
nThin <- 1
nIter <- (nBurn + nPost) * nThin
nBurnin <- nBurn * nThin

data <- read_rdump("twoCpt.data.R")

library(rstan)
## real[] time, real[] amt, int[] cmt, int[] evid

fit <- stan(file = file.path("./2cmt_model.stan", sep = ""),
            data = data,
            pars = parameters,
            iter = nIter,
            warmup = nBurnin,
            thin = nThin, 
            init = init,
            chains = nChains,
            cores = min(nChains, parallel::detectCores()))


save(fit, file = file.path("./",paste("this.", "Fit.Rsave", sep = "")))

stan_trace(fit, parametersToPlot)

pairs(fit, pars = parametersToPlot)

print(fit, pars = parametersToPlot)


xdata <- data.frame(data$cObs, data$time[data$evid != 1])

xdata <- plyr::rename(xdata, c("data.cObs"="cObs", "data.time.data.evid....1."="time"))

pred <- as.data.frame(fit, pars="cObsPred") %>%
  tidyr::gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols(xdata)

p1 <- ggplot(pred, aes(x=time, y=cObs))

p1 <- p1 + geom_point() + 
  labs(x="time (h)", y="plasma concentration (mg/L)") +
  theme(text = element_text(size=12), axis.text = element_text(size=12),
        legend.position = "none", strip.text = element_text(size=8))
p1 + geom_line(aes(x=time, y=median)) + geom_ribbon(aes(ymin=lb, ymax=ub), alpha=0.25)
