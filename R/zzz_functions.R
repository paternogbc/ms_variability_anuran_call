# Calculates the Coefficient of Variance (NA are removed from each vector)
cv <- function(x,na.rm=TRUE) 100*(sd(x,na.rm=na.rm)/mean(x,na.rm=na.rm))

# Function to gather summary (lmer)
sl <- function(x){
  
  # get estimates
  parameter <- x@call[[2]][[2]]
  predictor <- x@call[[2]][[3]][[2]]
  suppressWarnings({
    estimate  <- tidy(x)$estimate[[2]]
    p.value   <- tidy(Anova(x))$p.value
    Chisq     <- tidy(Anova(x))$statistic
  })
  
  r2_conditional    <- performance(x)$R2_conditional
  r2_marginal       <- performance(x)$R2_marginal
  
  out <- data.frame(parameter = paste(parameter), 
                    predictor = paste(predictor), 
                    estimate,
                    Chisq,
                    p.value, 
                    r2_conditional,
                    r2_marginal)
  return(out)
}

ext_mantel <- function(x, name) {
  data.frame(Source = name, r = x$obs, P_value = x$pvalue) 
  
} 
