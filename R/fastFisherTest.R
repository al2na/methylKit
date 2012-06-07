dyn.load(paste("fastFisherTest", .Platform$dynlib.ext, sep=""))

#fast.fisher.test <- function(a,b,c,d) {
#    n <- as.integer(length(a))
#    r <- .C("fet",as.integer(a),as.integer(b), as.integer(c), as.integer(d), n, out=double(n))$out
#    return (r);
#}

fast.fisher.test <- function(x) {
    if(is.matrix(x)){
      a=as.integer(x[,1])
      b=as.integer(x[,2])
      c=as.integer(x[,3])
      d=as.integer(x[,4])
      n <- as.integer(length(a))
      r <- .C("fastFisherTest",a,b,c,d,n, pval=double(n))$pval
      return (r);
    } else {
      show("input have to be a matrix with 4 columns")
    }

}
