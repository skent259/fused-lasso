
###-------------------------------------------------------###
###---- Setup                                            ----
###-------------------------------------------------------###
{
library(tidyverse)
library(png)
par(pty="s")
}
###-------------------------------------------------------###
###---- Images for testing                               ----
###-------------------------------------------------------###
{
image_noise <- function(image, sd, seed=8) {
  ## inputs:  image, a matrix with entries in [0.1] representing a greyscale image
  ##          sd, a standard deviation which determines the noise added.  
  ##            Scale can be given based a grey scale of 255
  ## output: image + gaussian noise
  
  ## adjust sd to appropriate scale
  if ( sd > 1 ) { 
    sd <- sd/255 
  }
  
  set.seed(seed)
  noisy <- image + sd * matrix( rnorm(prod(dim(image))), nrow = dim(image)[1])
  
  ## fix values between 0 and 1
  noisy[noisy < 0] <- 0
  noisy[noisy > 1] <- 1
  return( noisy )
}

## small (16x16) toy example
cross <- matrix(0.1, nrow=16, ncol=16)
cross[7:10,4:13] <- 0.8
cross[4:13,7:10] <- 0.8
cross_noi_s10 <- image_noise(cross, 10)
cross_noi_s40 <- image_noise(cross, 40)

## read in images from the BM3D paper (standard in image smoothing)
image_names <- c("Cameraman256", "house", "montage", "peppers256")
  
all_true_images <- lapply(image_names, FUN = function(x) {
  image <- readPNG( paste0("./BM3D_images/",x,".PNG") )
})

all_noisy_images <- lapply(image_names, FUN = function(x) {
  image <- readPNG( paste0("./BM3D_images/",x,".PNG") )
  s20 <- image_noise(image, 20)
  s40 <- image_noise(image, 40)
  s100 <- image_noise(image, 100)
  return(list(s20=s20, s40=s40, s100=s100))
})

names(all_noisy_images) <- names(all_true_images) <- image_names
}
###-------------------------------------------------------###
###---- Functions                                        ----
###-------------------------------------------------------###
{
soft_threshold <- function(image, penalty) {
  sign(image) * max( abs(image)-penalty, 0 )
}

find_d1_groups <- function(groups_, group_one_) {
  ## returns groups that are touching the selected group
  
  nrow <- dim(groups_)[1]
  ncol <- dim(groups_)[2]
  ind <- which(groups_ == group_one_)
  
  d1_groups_ind <- c(ind[row(groups_)[ind] != nrow] + 1, # index below all not in last row
                     ind[row(groups_)[ind] != 1]    - 1, # index above all not in first row
                     ind[col(groups_)[ind] != ncol] + nrow, # index right of all not in last column
                     ind[col(groups_)[ind] != 1]    - nrow) %>% # index left of all not in first column
    setdiff(ind) %>%
    unique()
  
  out <- groups_[d1_groups_ind] 
  
  return(out) 
}

fuse_groups <- function(groups_, group_one_, group_two_) {
  groups_[groups_==group_two_] <- group_one_
  return(groups_)
}

w <- function(groups_, group_one_, group_two_) {
  ind_one <- which(groups_==group_one_)
  ind_two <- which(groups_==group_two_)
  nrow <- dim(groups_)[1]
  
  ind_d1_to_one <- c(ind_two+1, ind_two-1, ind_two+nrow, ind_two-nrow)
  return(sum(ind_one %in% ind_d1_to_one))
}

objective_function <- function(image_, estimate_, groups_, lam1_, lam2_) {
  mse <- 1/2  * sum( (image_ - estimate_)^2 ) 
  l1 <- lam1_ * sum( abs(estimate_) )
  
  fused <- lam2_/2 * sum( sapply( unique(as.vector(groups_)), FUN = function(k) {
    ## first, sum over each group k
    d1_groups <- groups_ %>% find_d1_groups(k)
    gamma_k <- estimate_[which(groups_==k)][1]
    
    out_k <- sum( sapply(d1_groups, FUN = function(k_prime) {
      ## second, sum over k_prime that are distance 1 from k
      gamma_k_prime <- estimate_[which(groups_==k_prime)][1]
      
      out_k_prime <- w(groups_, k, k_prime) * abs( gamma_k - gamma_k_prime )
      return(out_k_prime)
    }))
    return(out_k)
  }))
  
  out <- mse + l1 + fused
  return(out)
}

local_objective_function <- function(k, image_, estimate_, groups_, lam1_, lam2_) {
  mse <- 1/2  * sum( (image_[groups_==k] - estimate_[groups_==k])^2 ) 
  l1 <- lam1_ * sum( abs(estimate_[groups_==k]) )
  
  d1_groups <- groups_ %>% find_d1_groups(k)
  gamma_k <- estimate_[which(groups_==k)][1]
  
  fused <- lam2_/2 * sum( sapply(d1_groups, FUN = function(k_prime) {
    ## second, sum over k_prime that are distance 1 from k
    gamma_k_prime <- estimate_[which(groups_==k_prime)][1]
    
    out_k_prime <- w(groups_, k, k_prime) * abs( gamma_k - gamma_k_prime )
    return(out_k_prime)
  }))
  
  out <- mse + l1 + fused
  return(out)
}

gradient_objective_function_k <- function(k_, gamma_k_=NULL, image_, estimate_, groups_, lam1_, lam2_) {
  ## k_ is group of interest, gamma_k_ can be set as any estimate
  ## image_: original image; estimate_: current estimate (gives gamma_k's)
  if(is.null(gamma_k_)) { 
    gamma_k_ <- estimate_[groups_==k_] 
  }
  
  mse <- -sum( (image_[groups_==k_] - gamma_k_) ) # \bar{y_k} - gamma_k
  l1 <- lam1_ * sum( sign(gamma_k_) )
  
  d1_groups <- groups_ %>% find_d1_groups(k_)
  
  fused <-  lam2_/2 * sum( sapply(d1_groups, FUN = function(k_prime) {
    ## sum over k_prime that are distance 1 from k
    gamma_k_prime <- estimate_[which(groups_==k_prime)][1]
    
    out_k_prime <- w(groups_, k_, k_prime) * sign( gamma_k_ - gamma_k_prime )
    return(out_k_prime)
  }))
  
  
  out <- mse + l1 + fused
  return(out)
}

find_gamma_k <- function(k_, image_, estimate_, groups_, lam1_, lam2_) {
  ## minimize objective function over gamma_k, holding all other gamma_k_prime as fixed
  ## This utilizes the fact that the gradient of the objective function is
  ##  piecewise linear with breaks at 0 and gamma estimates of distance 1 neighbors.
  
  d1_groups_ <- groups_ %>% find_d1_groups(k_) %>% unique()
  gamma_k_prime_ <- unique(estimate_[groups_ %in% d1_groups_])
  eps <- 1e-4
  breaks_ <- sort( c( gamma_k_prime_-eps, gamma_k_prime_+eps, 0) )
  
  gradient_at_breaks_ <- sapply(breaks_, FUN = gradient_objective_function_k, 
                                k_=k_, image_=image_, estimate_=estimate_, 
                                groups_=groups_, lam1_=lam1_, lam2_=lam2_)
  
  ## calculate min of objective function
  if(sum(unique(sign(gradient_at_breaks_))) == 0) { # gradient crosses zero
    
    ind <- max( which(gradient_at_breaks_ < 0) )  # index just before graidient crosses zero 
    
    # x = - y_1 * (x_2 - x_1) / (y_2 - y_1) + x_1 gives x at which line hits 0
    find_zero <- function(x,y) {
      slope <- (y[2] - y[1]) / (x[2] - x[1])
      zero <- - y[1] / slope + x[1]
      return(zero)
    }
    
    min_ <- find_zero(breaks_[ind:(ind+1)], gradient_at_breaks_[ind:(ind+1)])
    
  } else { # gradient doesn't cross zero
    
    ## Make more efficient by taking advantage of (39)!!!!!!
    estimate_w_k_prime_ <- lapply(d1_groups_, FUN = function(k_prime_) {
      estimate_[groups_==k_prime_][1] * (groups_ == k_) + estimate_ * (groups_ != k_)
    })
    estimate_w_0 <- 0 * (groups_ == k_) + estimate_ * (groups_ != k_)
    
    ## take gamma_k_prime with smallest objective function
    obj_at_breaks_ <- sapply(estimate_w_k_prime_, FUN = local_objective_function, 
                             k=k_, image_=image_, groups_=groups_, lam1_=lam1_, lam2_=lam2_)
    obj_at_zero <- local_objective_function(k_, image_, estimate_w_0, groups_, lam1_, lam2_)
    
    if ( min(obj_at_breaks_) < obj_at_zero ) {
      k_min <- d1_groups_[obj_at_breaks_ == min(obj_at_breaks_)]
      min_ <- estimate_[groups_==k_min][1]
    } else {
      min_ <- 0
    }
    
  }
  
  return(min_)
  
}

diff_in_objective_function <- function(gamma_m_, k_, k_prime_, image_, estimate_, groups_, lam1_, lam2_) {
  gamma_k_ <- estimate_[groups_==k_][1]
  gamma_k_prime_ <- estimate_[groups_==k_prime_][1]
  n_k_ <- sum(groups_==k_)
  n_k_prime_ <- sum(groups_==k_prime_)
  y_k_ <- mean(image_[groups_==k_])
  y_k_prime_ <- mean(image_[groups_==k_prime_])
  ##w_kl_ <- w(groups_, k_, )
  
  n_m_ <- n_k_ + n_k_prime_
  y_m_ <- (n_k_*y_k_ + n_k_prime_*y_k_prime_) / n_m_
  
  mse_diff <- (-n_k_*(y_k_ - gamma_k_)^2 +
                 -n_k_prime_*(y_k_prime_ - gamma_k_prime_)^2 +
                 n_m_*(y_m_ - gamma_m_)^2 ) 
  
  l1_diff <- lam1_ * (-n_k_*abs(gamma_k_) +
                        -n_k_prime_*abs(gamma_k_prime_) +
                        n_m_*abs(gamma_m_) ) 
  
  fused_k <- function(l,k,gamma_k,groups) {
    ## for given k, and l with d(k,l)=1, finds fused penalty component
    gamma_l <- estimate_[which(groups==l)][1]
    
    out_k_prime <- w(groups, k, l) * abs( gamma_k - gamma_l )
    return(out_k_prime)
  }
  
  d1_groups_k <- find_d1_groups(groups_, k)
  d1_groups_k_prime <- find_d1_groups(groups_, k_prime_)
  d1_groups_m <- union(d1_groups_k, d1_groups_k_prime) %>% setdiff( c(k_,k_prime_) )
  groups_fused_ <- fuse_groups(groups_, k_, k_prime_)
  
  fused_diff <- lam2_/2 * ( -sum( sapply(d1_groups_k, FUN = fused_k, k=k_, gamma_k=gamma_k_, groups=groups_) ) +
                              -sum( sapply(d1_groups_k_prime, FUN = fused_k, k=k_prime_, gamma_k=gamma_k_prime_,groups=groups_) ) +
                              sum( sapply(d1_groups_m, FUN = fused_k, k=k_, gamma_k=gamma_k_, groups=groups_fused_) ) )
  
  return(mse_diff + l1_diff + fused_diff)
}
}
###-------------------------------------------------------###
###---- Main Algorithm                                   ----
###-------------------------------------------------------###

fused_lasso <- function(image, lam1, lam2=0, lam2_max, delta, tol=5e-4, verbose=FALSE) {
  
  ## start groups as all separate
  groups <- matrix(1:prod(dim(image)), nrow=dim(image)[1])
  
  estimate_list <- list()
  i <- 1
  
  while (lam2 < lam2_max) {
    if (lam2 == 0) { 
      ## set initial estimate to soft-thresholded value of image
      estimate <- soft_threshold(image, lam1)
      estimate_list[[i]] <- estimate # record the estimate
      i <- i + 1
      lam2 <- lam2 + delta
    } 
    
    if(verbose) { cat("\n Lambda2:",lam2,"\n") }
    gamma_change <- TRUE
    while (gamma_change) { # keep track of whether we change any gamma
      
      gamma_k_change <- logical(length(unique(as.vector(groups))))
      
      for (k in unique(as.vector(groups))) {
        gamma_k_change[k] <- FALSE # start with no change
        
        if ( !(k %in% as.vector(groups)) ) next # needed if we fuse a group
        
        ###############.
        ## Descent
        gamma_k_new <- find_gamma_k(k_ = k, image, estimate, groups, lam1, lam2)
        gamma_k <- estimate[groups==k][1]
        
        if ( abs(gamma_k_new - gamma_k) > tol ) {
          estimate[groups == k] <- gamma_k_new
          gamma_k_change[k] <- TRUE
          
        } else { # gamma_k_new == gamma_k 
          ###############.
          ## Fusion 
          d1_groups <- find_d1_groups(groups, k)
          #obj <- objective_function(image, estimate, groups, lam1, lam2)
          
          for (k_prime in d1_groups) {
            
            groups_tmp <- fuse_groups(groups, k, k_prime) # provisionally fuse groups k and k_prime
            gamma_m <- find_gamma_k(k_ = k, image, estimate, groups_tmp, lam1, lam2)
            if (is.na(gamma_m)) { 
              break ## not sure why this would happen...
            }
            
            estimate_tmp <- estimate
            estimate_tmp[groups %in% c(k,k_prime)] <- gamma_m
            
            diff <- local_objective_function(k, image, estimate_tmp, groups_tmp, lam1, lam2) -
              local_objective_function(k, image, estimate, groups, lam1, lam2) -
              local_objective_function(k_prime, image, estimate, groups, lam1, lam2)
            improved_criterion <- diff < 0
            
            if (improved_criterion) {
              groups <- groups_tmp
              estimate[groups_tmp %in% c(k,k_prime)] <- gamma_m
              gamma_k_change[k] <- TRUE
              break # ends for k_prime
            } # else keep gamma_k the same
            
          } # end for k_prime
          
        } # end else
        
        
      } # end for k
      
      
      gamma_change <- sum(gamma_k_change, na.rm = TRUE) > 0 # logical for whether ANY gamma changed
      
      ## Display prevalent metrics
      if(verbose) {
        if (prod(dim(image)) < 300) { obj <- objective_function(image, estimate, groups, lam1, lam2) } 
        n_groups <- length(unique(as.vector(groups)))
        cat("  Cycle completed.")
        if (prod(dim(image)) < 300) { cat(" Obj fun:", obj) }
        cat("  Ngroups:", n_groups ,"Nchanged:",sum(gamma_k_change, na.rm = TRUE),"\n")
        image(estimate, col=grey.colors(255)) # view estimate
      }
    } # end while (will break when no gammas change in a cycle of descent & fusion)
    
    ###############.
    ## Smoothing
    estimate_list[[i]] <- estimate # record the estimate
    i <- i + 1
    lam2 <- lam2 + delta
    
    if (length(unique(as.vector(groups))) == 1) { break } ## stop iterating if we have one group
  }
  
  return(estimate_list)
}

###-------------------------------------------------------###
###---- Implementation                                   ----
###-------------------------------------------------------###
system.time({

all_noisy_images_c <- list(cross_noi_s10, cross_noi_s40)
estimates.cross <- lapply(all_noisy_images_c, FUN = fused_lasso,
                          lam1=0, lam2=0, lam2_max=0.1, delta=0.02, verbose=TRUE)

save(estimates.cross, "cross_estimates.RData")

})
###-------------------------------------------------------###
###---- Evaulation                                       ----
###-------------------------------------------------------###

## Peak Signal-to-Noise Ratio
psnr <- function(y_true, y_est) {
  10*log(base = 10, x = ( 1/mean((y_true-y_est)^2) ) )
}

image(estimates.cross[[1]][[5]], col=grey.colors(255))
image(estimates.cross[[2]][[5]], col=grey.colors(255))
image(cross_noi_s10, col=grey.colors(255))

psnr_s10 <- sapply(estimates.cross[[1]], FUN = psnr, y_true=cross)
psnr_s40 <- sapply(estimates.cross[[2]], FUN = psnr, y_true=cross)


###-------------------------------------------------------###
###---- Implementation  (in flsa package)                ----
###-------------------------------------------------------###
{
library(flsa)

lambda2_seq <- seq(0,0.1,length=5)

# tmp <- flsa(all_noisy_images$house$s20, lambda1=0, lambda2=lambda2_seq, verbose=TRUE)

system.time({
  estimates.all <- lapply(all_noisy_images, FUN = flsa,
                          lambda1=0, lambda2=lambda2_seq, verbose=TRUE)
})

save(estimates.all, "all_estimates.RData")

}
###------------------------------------------------------####






