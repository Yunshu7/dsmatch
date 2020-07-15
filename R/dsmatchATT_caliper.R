dsmatchATT_caliper = function(Y, A, X, caliper, replace){
  pg = X[,1]
  ps = X[,2]
  Sinv = solve(cov(X)) # inverse of covariance matrix

  loc.1 = which(A == 1)
  loc.0 = which(A == 0)
  candidate_pool = loc.0 # will change if replace = F

  ndrops = 0
  index.treated = c()
  index.control = c()

  # llist = c()

  for (i in loc.1) {
    if(length(candidate_pool) > 0){
      # satisfy caliper restriction for subject i
      # either score is within caliper
      pg_in = abs(pg[candidate_pool] - pg[i]) <= caliper
      ps_in = abs(ps[candidate_pool] - ps[i]) <= caliper
      candidate_i = candidate_pool[pg_in | ps_in]

      # llist = c(llist, length(candidate_i))

      if(length(candidate_i) > 0){
        # calculate Mahalanobis distance
        if(length(candidate_i) > 1){
          distance_i = t(t(X[candidate_i, ]) - X[i, ]) %*% Sinv %*% (t(X[candidate_i, ]) - X[i, ])
          distance_i = diag(distance_i)

          # choose the closest control subject to match
          match.control.index = candidate_i[which.min(distance_i)]
        }else{
          distance_i = (X[candidate_i, ] - X[i, ]) %*% Sinv %*% (X[candidate_i, ] - X[i, ])
        }
        # add matched pair into the list
        index.control = c(index.control, match.control.index)
        index.treated = c(index.treated, i)

        # if no replacement, remove the matched control subject from candidate pool
        if (replace == F){
          index.remove = which(candidate_pool == match.control.index)
          candidate_pool = candidate_pool[-index.remove]
        }
      }else{
        # if no subject is within the caliper, drop it
        ndrops = ndrops + 1
      }
    }else{
      # if no candidate left, drop it
      ndrops = ndrops + 1
    }
  }

  return(list(index.control = index.control, index.treated = index.treated, ndrops = ndrops))
}
