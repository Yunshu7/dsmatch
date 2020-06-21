is.pg = function(pg, n){
  if(prod(dim(pg) == c(n, 2))){
    return(TRUE)
  }else{
    return(FALSE)
  }
}