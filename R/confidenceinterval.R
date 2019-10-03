# confidence interval for RI
# res is the output of CalculateSens()

ConfidenceInterval <- function(res, niter = 100) {
  dose = sort(unique(df$dose))
  n = length(dose)
  score.simulated = c()
  for(i in 1:niter){
    drug.simulated = data.frame(response = rnorm(n, mean = pred$Prediction, sd = pred$sd), dose = dose)
    res.simulated = CalculateSens(drug.simulated)
    score.simulated[i] = res.simulated$score
  }

  # confidence interval
  lower = quantile(score.simulated, probs = 0.025)
  upper = quantile(score.simulated, probs = 0.975)
  score = list(score = res$score, lower = lower, upper = upper)

  #options(show.error.messages = TRUE)
  return (score)

}

