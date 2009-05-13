print.tclust <-
function (x, ...)
{
    cat("* Results for TCLUST algorithm: *\n")
    cat("trim=", x$par$alpha, ", k=", x$k, "\n")
    cat("Classification (trimmed points are indicated by 0", "):\n")
    print(x$assig)
    cat("Means:\n")
    print(x$center)
    if (x$obj < (-1e+20)) 
        warning("The solution is not reliable. More iterations are probably needed.")
    cat("Trimmed objective function: ", x$obj, "\n")	
}

