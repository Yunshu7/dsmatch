is.lp = function(lp, n){
  if(is.vector(lp)){
    if(length(lp) == n){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }else{
    if(dim(lp)[1] == n){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }
}