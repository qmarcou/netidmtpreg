logit <- function(p) log(p / (1.0 - p))
expit <- function(x) 1.0 / (1.0 + exp(-x))