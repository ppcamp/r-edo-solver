# # Numerical Calculation"s methods class.
# **NOTE** observes that all methods have _VERBOSE_ parameter,
# this parameter can be used to see steps inside loops && another
# infos. Also they have an _VECTOR_ param, wich gives same info as
# Verbose, returning in Data-frame

# rm(list = ls()) # Clear all local variables

# To do ----------------------------------
# [X] 1. Check exibitions inside methods
#    (verbose flag).
# [±] 2. Create & update documentation
# [?] 3. Make plots from differences
# [#] 4. Do report
# [ ] 5. Fix Secant method

# Formating ------------------------------
hline <- function(value = 2){
    cat(sprintf("------------", sep="", 1:value), sep="")
    cat('\n')
}
hfill <- function(n = 1){
    while (n) {
        cat("\n")
        n=n-1
    }
}


# Constants ------------------------------
# Set default tolerance
self.tol = 10^(-15)
# Set default NUMBER of decimal places
self.DECIMAL = 13


# Others methods -------------------------
about <- function(){
    cat("CalcNum.\nVersion: 1.0\n\nCommands:\n")
    cat("VERBOSE=1\tTo see all infos\n")
}

getFuncValue <- function(kwg, ...){
    # Unpacking parameters
    args = list(...)
    func = args$func

    if ("x" %in% names(args)){func = gsub("x", args$x, func)}
    if ("y" %in% names(args)){func = gsub("y", args$y, func)}
    return(eval( parse(text=func) ))
}


# Solutions of Eqs in One Variable -------
minSteps <- function (kwg,...){
    args = list(...)
    n = (log10(args$b-args$a) - log10(args$E)) / log10(2)
    return (ceiling(n))
}

bisection <- function(kwg, ...){
    # Unpacking parameters
    args = list(...)

    # Tolerance
    tol = if ("tol" %in% names(args)) args$tol else self.tol

    # Decimal places
    DECIMAL = if ("DECIMAL" %in% names(args)) args$DECIMAL else self.DECIMAL

    # Function in point
    fa = round(getFuncValue(func=args$func, x=args$a), digits = DECIMAL)
    pi = fp = 0.0

    if ("VERBOSE" %in% names(args) && args["VERBOSE"]==1){
        cat(sprintf("Function: \t%s\nTol:\t\t%.1E\n", args$func,tol))
        hline(6)
        cat(sprintf("i\tan\t\tbn\t\tpn\t\tf(pn)"), sep = '\n')
        hline(6)
    }

    if ("VECTOR" %in% names(args) && args$VECTOR==1)
        m = matrix(ncol=5)

    for (i in 1:args$n0){
        pi = round((args$a+args$b)/2, digits = DECIMAL)
        fp = round(getFuncValue(func=args$func, x=pi, format=0), digits = DECIMAL)

        if ("VERBOSE" %in% names(args) && args["VERBOSE"]==1)
            cat(sprintf("%d\t%.9f\t%.9f\t%.9f\t%.9f", i,args$a,args$b,pi,fp), step='\n')

        if ("VECTOR" %in% names(args) && args$VECTOR==1) {
            m <- rbind(m, round(c(i, args$a, args$b, pi, fp), digits=9) )
            if (i == args$n0) {
                df = data.frame(An=m[2:nrow(m),2],
                                Bn=m[2:nrow(m),3], Pn=m[2:nrow(m),4],
                                F_Pn=m[2:nrow(m),5])
                m = NULL
                return (df)
            }
        }
        if (fp==0 || ((args$b-args$a)/2 < tol)) return(pi)

        if( fa*fp > 0) {
            args$a = pi
            fa = fp
        }
        else args$b = pi

    }
    return('Fail')
}

newtonIteration <- function(kwg, ...) {
    # Unpacking parameters
    args = list(...)

    # Tolerance
    tol = if ("tol" %in% names(args)) args$tol else self.tol

    # Decimal places
    DECIMAL = if ("DECIMAL" %in% names(args)) args$DECIMAL else self.DECIMAL

    if ("VERBOSE" %in% names(args) && args["VERBOSE"]==1){
        cat(sprintf("Function: \t%s\nTol:\t\t%.1E\n", args$func,tol))
        cat(sprintf("Function: \t%s\n", args$funcDeriv))
        hline()
        cat(sprintf("i\tp"), sep = '\n')
        hline()
        cat(sprintf("%d\t%.9f", 0,args$p0), sep = '\n')

    }

    if ("VECTOR" %in% names(args) && args$VECTOR==1) {
        # Matrix
        m = matrix(ncol = 2)
        m <- rbind(m, round(c(0, args$p0), DECIMAL) )
    }

    for (i in 1:args$n0) {
        # Note:     Package deriv not working
        # Error:    Error in deriv.default(a,'x') :
        #               invalid expression in 'FindSubexprs'
        aux = getFuncValue(func=args$func, x=args$p0) / getFuncValue(func=args$funcDeriv, x=args$p0)
        p = args$p0 - aux

        if ("VERBOSE" %in% names(args) && args["VERBOSE"]==1)
            cat(sprintf("%d\t%.9f", i,p), sep = '\n')

        if ("VECTOR" %in% names(args) && args$VECTOR==1) {
            m <- rbind(m, round(c(i, p), DECIMAL) )
            if (i == args$n0) {
                df = data.frame(Pn=m[2:nrow(m),2])
                m = NULL
                return (df)
            }
        }
        if (abs(p - args$p0) < tol) return(p)
        args$p0 = p
    }
    return('Fail')
}

fixedPointIteration <- function(kwg, ...) {
    # Unpacking parameters
    args = list(...)

    # Tolerance
    tol = if ("tol" %in% names(args)) args$tol else self.tol

    # Decimal places
    DECIMAL = if ("DECIMAL" %in% names(args)) args$DECIMAL else self.DECIMAL

    if ("VERBOSE" %in% names(args) && args["VERBOSE"]==1){
        cat(sprintf("Function: \t%s\nTol:\t\t%.1E\n", args$func,tol))
        hline(2)
        cat(sprintf("i\t p"), sep = '\n')
        hline(2)
        cat(sprintf("%d\t %.9f", 0, args$p0), sep = '\n')
    }

    if ("VECTOR" %in% names(args) && args$VECTOR==1){
        m = matrix(ncol=2)
        m <- rbind(m, round(c(1,args$p0),DECIMAL) )
    }

    for (i in 1:args$n0) {
        p = getFuncValue(func=args$func, x=args$p0)

        if ("VECTOR" %in% names(args) && args$VECTOR==1) {
            m <- rbind(m, round(c(i, p), DECIMAL) )
            if ((i == args$n0) || (abs(p - args$p0) < tol)) {
                df = data.frame(Pn=m[2:nrow(m),2])
                m = NULL
                return (df)
            }
        }
        if ("VERBOSE" %in% names(args) && args["VERBOSE"]==1)
            cat(sprintf("%d\t %.9f", i, p), sep='\n')

        if (abs(p - args$p0) < tol) return(p)
        args$p0 = p
    }
    return('Fail')
}

falsePosition <- function(kwg, ...) {
    # Uncompressing list
    args = list(...)

    # Tolerance
    tol = if ("tol" %in% names(args)) args$tol else self.tol

    # Decimal places
    DECIMAL = if ("DECIMAL" %in% names(args)) args$DECIMAL else self.DECIMAL

    if ("VERBOSE" %in% names(args) && args["VERBOSE"]==1){
        cat(sprintf("Function: \t%s\nTol:\t\t%.1E\nN0\t%d", args$func,tol, args$n0), sep = '\n')
        hline(2)
        cat(sprintf("n\tpn"), sep='\n')
        hline(2)
        cat(sprintf("%d\t%.9f",0, args$p0), sep = '\n')
        cat(sprintf("%d\t%.9f",1, args$p1), sep = '\n')
    }

    if ("VECTOR" %in% names(args) && args$VECTOR==1) {
        m = matrix(ncol=2)
        m <- rbind(m, round(c(1, args$p0), DECIMAL) )
        m <- rbind(m, round(c(2, args$p1), DECIMAL) )
    }

    q0 = round(getFuncValue(func=args$func, x=args$p0), DECIMAL)
    q1 = round(getFuncValue(func=args$func, x=args$p1), DECIMAL)

    for (i in seq(from=2, to=args$n0)) {
        p = round(args$p1 - q1*(args$p1 - args$p0)/(q1 - q0), DECIMAL)

        if ("VERBOSE" %in% names(args) && args["VERBOSE"]==1)
            cat(sprintf("%d\t%.9f",i, p), sep = '\n')

        if ("VECTOR" %in% names(args) && args$VECTOR==1) {
            m <- rbind(m, round(c(i+1, p), DECIMAL) )
            if ((i == args$n0) || (abs(p - args$p1) < tol)) {
                df = data.frame(Pn=m[2:nrow(m),2])
                m = NULL
                return (df)
            }
        }

        if (abs(p - args$p1) < tol) return(p)

        q = round(getFuncValue(func=args$func, x=p), DECIMAL)
        if (q*q1 < 0) {
            args$p0 = round(args$p1, DECIMAL)
            q0 = round(q1, DECIMAL)
        }
        args$p1 = p
        q1 = q
    }
    return('Fail')
}

# Differential Equations -----------------
odeEuler <- function(kwg, ...){
    # Unpacking parameters
    args = list(...)

    h = (args$b-args$a)/args$m
    x = args$a
    y = args$y0

    fxy = getFuncValue(func=args$func, x=x,y=y)

    m = matrix(ncol=3)

    if ("VERBOSE" %in% names(args) && args$VERBOSE==1){
        cat("Values from Euler's method\n")
        cat("i \t x \t\t y \t\t fxy\n")
        cat(sprintf("%i\t%.5f\t\t%.5f\t\t%.5f\n", 0,x,y,fxy))
    }
    if ("VECTOR" %in% names(args) && args$VECTOR==1){
        m = matrix(ncol=4)
        m <- rbind(m, round(c(0, x, y, fxy), digits=9) )
    }

    for (i in 1:args$m){
        x = args$a + i * h
        y = y + h * fxy
        fxy = getFuncValue(func=args$func, x=x,y=y)
        if ("VERBOSE" %in% names(args) && args$VERBOSE==1){
            cat(sprintf("%i\t%.5f\t\t%.5f\t\t%.5f\n", i,x,y,fxy))
        }
        if ("VECTOR" %in% names(args) && args$VECTOR==1) {
            m <- rbind(m, round(c(i, x, y, fxy), digits=9) )
            if (i == args$m) {
                df = data.frame(i=m[2:nrow(m),1], x=m[2:nrow(m),2],
                                y=m[2:nrow(m),3], fxy=m[2:nrow(m),4])
                m = NULL
                return (df)
            }
        }
    }
    c(x,y,fxy)
}

odeRK4 <- function(kwg, ...){
    # Unpacking parameters
    args = list(...)

    h = (args$b-args$a)/args$m
    xt = args$a
    yt = args$y0

    if ("VERBOSE" %in% names(args) && args$VERBOSE==1){
        cat("Values from RK4's method\n")
        cat("i \t x \t\t y\n")
        cat(sprintf("0 \t %.5f \t %.5f\n", xt,yt))
    }

    if ("VECTOR" %in% names(args) && args$VECTOR==1){
        m = matrix(ncol=3)
        m <- rbind(m, round(c(0, xt, yt), digits=9) )
    }

    for (i in 1:args$m){
        x = xt
        y = yt
        k1 = getFuncValue(func=args$func, x=x,y=y)

        x = xt + h/2
        y = yt + (h/2)*k1
        k2 = getFuncValue(func=args$func, x=x,y=y)

        y = yt + (h/2)*k2
        k3 = getFuncValue(func=args$func, x=x,y=y)

        x = xt + h
        y = yt + h*k3
        k4 = getFuncValue(func=args$func, x=x,y=y)

        xt = args$a + i*h
        yt = yt + (h/6)*(k1 + 2*(k2 + k3) + k4)

        if ("VERBOSE" %in% names(args) && args$VERBOSE==1){
            cat(sprintf("%i \t %.5f \t %.5f\n",i,xt,yt))
        }
        if ("VECTOR" %in% names(args) && args$VECTOR==1) {
            m <- rbind(m, round(c(i, xt, yt), digits=9) )
            if (i == args$m) {
                df = data.frame(i=m[2:nrow(m),1], xt=m[2:nrow(m),2],
                                yt=m[2:nrow(m),3])
                m = NULL
                return (df)
            }
        }
    }
    c(xt, yt)
}

# Tests and analysis ---------------------
euler <- odeEuler(func='x-2*y+1',a=0,b=1,y0=1,m=10,VECTOR=1)
rk4 <- odeRK4(func='x-2*y+1',a=0,b=1,y0=1,m=10,VECTOR=1)

x <- euler[1:11, 2]
y_exata <- (1/4)*(3*exp(-2*x)+2*x+1) # Função resolvida analiticamente
y_euler <- euler[1:11, 3]
y_rk4   <- rk4[1:11, 3]

DRACULA_ORCHID= "#2d3436"
WET_ASPHALT = "#34495e"
ORANGE_VILLE = "#e17055"
BRIGHT_YARROW = "#fdcb6e"
MINT_LEAF = "#00b894"
PRUNUS_AVIUM = "#e84393"
ELECTRON_BLUE = "#0984e3"
EXODUS_FRUIT = "#6c5ce7"
cor_exata = ELECTRON_BLUE
cor_rk4   = DRACULA_ORCHID
cor_euler = EXODUS_FRUIT

plot(frame = FALSE,
     ylim=c(0.75,1),
     lwd = 1,
     pch = 19,
     x, y_exata,
     xlab = "x",
     ylab = "Y",
     type = "b",
     col = cor_exata,
     main = "Diferença entre RK4, Euler e o valor Exato"
)
lines(x, y_rk4, col=cor_rk4, lwd = 1)
lines(x, y_euler, col=cor_euler, lwd = 1)
grid(lty="solid")#longdash")
legend(
    "topright",
    c("Exata","RK4","Euler"),
    fill = c(cor_exata,cor_rk4,cor_euler)
)


# https://subscription.packtpub.com/book/big_data_and_business_intelligence/9781849513067/4/ch04lvl1sec06/adding-horizontal-and-vertical-grid-lines
