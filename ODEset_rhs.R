#builds right hand side of the ODE set
ODEset_rhs = function (t, y, parms){
  with(as.list(y), {
    starting_vector=y
    transition_matrix=parms
    dy=rep(0,length(starting_vector))
    dy=drop(transition_matrix %*% starting_vector)
    out=(dy)
    names(out)=names(y)
    return(list(out))
  })
}


min(result)
hist(result)
