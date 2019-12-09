Cohens.d <- function(Group.1, Group.2) {
  n.1 <- length(Group.1)
  n.2 <- length(Group.2)
  SS1 <- var(Group.1)*(n.1-1)
  SS2 <- var(Group.2)*(n.2-1)
  pooled.sd <- sqrt((SS1 + SS2)/(n.1+n.2-2))
  Result <- (mean(Group.1)-mean(Group.2))/pooled.sd 
  Result
}

Unbiased.d <- function(Group.1, Group.2){
  nu <- length(Group.1) + length(Group.2) - 2
  G.nu <- gamma(nu / 2) / (sqrt(nu / 2) * gamma((nu - 1) / 2)) 
  d <- Cohens.d(Group.1, Group.2)
  Result <- ifelse(nu > 171, d, d * G.nu)
  Result
}

bootstrapDCI <- function(Group.1, Group.2, B = 1000, alpha = 0.05){
  n.1 <- length(Group.1) 
  n.2 <- length(Group.2)
  
  if(is.nan(Cohens.d(Group.1,Group.2)) | abs(Cohens.d(Group.1,Group.2)) == Inf  ){
    data.frame(method = c('percentile', 'bca', 'bootLength'),
               lcl = c(NA, NA, NA),
               ucl = c(NA, NA, NA)
    )
  }else{
  
  Bootstrap.Results <- matrix(NA, B, 1)
  # for(b in 1:B){
  #   Bootstrap.Results[b,1] <- Unbiased.d(sample(Group.1, size=n.1, replace=T),
  #                                        sample(Group.2,size=n.2, replace=T))
  # }
  # 
  # Bootstrap.Results <- Bootstrap.Results %>% na.omit()
  
  while(length(which(!is.na(Bootstrap.Results))) < B){
    ub <- Unbiased.d(sample(Group.1, size=n.1, replace=T),
                     sample(Group.2,size=n.2, replace=T))
    
    if(is.numeric(ub)){
      Bootstrap.Results[length(which(!is.na(Bootstrap.Results)))+1,1] <- ub
    }
    
  }
  
  
  Jackknife.Results <- matrix(NA,n.1+n.2,1)
  Marker.1 <- seq(1, n.1, 1)
  for(sample.1 in 1:n.1){
    Jackknife.Results[sample.1, 1] <- Unbiased.d(Group.1[Marker.1 [-sample.1]],Group.2)
  }
  
  Marker.2 <- seq(1, n.2, 1)
  for(sample.2 in 1:n.2) {
    Jackknife.Results[n.1+sample.2, 1] <- Unbiased.d(Group.1,
                                                     Group.2[Marker.2[-sample.2]]) 
  }
  
  Jackknife.Results <- Jackknife.Results %>% na.omit()
  
  Mean.Jackknife <- mean(Jackknife.Results)
  
  a <- (sum((Mean.Jackknife-Jackknife.Results)^3))/(6*sum((Mean.Jackknife-Jackknife.Results)^2)^(3/2))
  
  propLessD <- ifelse(Bootstrap.Results < Unbiased.d(Group.1, Group.2), 1, 0)
  
  z0 <- qnorm(sum(propLessD)/B)
  
  CI.Low.BCa <- pnorm(z0 + (z0+qnorm(alpha/2))/(1-a*(z0+qnorm(alpha/ 2))))
  CI.Up.BCa <- pnorm(z0 + (z0+qnorm(1-alpha/2))/(1-a*(z0+qnorm(1- alpha/2))))
  
  Percentile.Confidence.Limits <- c(quantile(Bootstrap.Results,alpha/2),
                                    quantile(Bootstrap.Results, 1-alpha/2))
  
  BCa.Confidence.Limits <- c(quantile(Bootstrap.Results, CI.Low.BCa),
                             quantile(Bootstrap.Results, CI.Up.BCa))
  
  data.frame(method = c('percentile', 'bca', 'bootLength'),
             lcl = c(Percentile.Confidence.Limits[1], ifelse(is.nan(BCa.Confidence.Limits[1]), NA, BCa.Confidence.Limits[1]), nrow(Bootstrap.Results)),
             ucl = c(Percentile.Confidence.Limits[2],ifelse(is.nan(BCa.Confidence.Limits[2]), NA, BCa.Confidence.Limits[2]), nrow(Bootstrap.Results))
  )
  }
  
}
