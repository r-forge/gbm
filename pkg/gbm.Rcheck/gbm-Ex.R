### * <HEADER>
###
attach(NULL, name = "CheckExEnv")
assign("nameEx", 
       local({
	   s <- "__{must remake R-ex/*.R}__"
           function(new) {
               if(!missing(new)) s <<- new else s
           }
       }),
       pos = "CheckExEnv")
## Add some hooks to label plot pages for base and grid graphics
assign("base_plot_hook",
       function() {
           pp <- par(c("mfg","mfcol","oma","mar"))
           if(all(pp$mfg[1:2] == c(1, pp$mfcol[2]))) {
               outer <- (oma4 <- pp$oma[4]) > 0; mar4 <- pp$mar[4]
               mtext(sprintf("help(\"%s\")", nameEx()), side = 4,
                     line = if(outer)max(1, oma4 - 1) else min(1, mar4 - 1),
               outer = outer, adj = 1, cex = .8, col = "orchid", las=3)
           }
       },
       pos = "CheckExEnv")
assign("grid_plot_hook",
       function() {
           grid::pushViewport(grid::viewport(width=grid::unit(1, "npc") - 
                              grid::unit(1, "lines"), x=0, just="left"))
           grid::grid.text(sprintf("help(\"%s\")", nameEx()),
                           x=grid::unit(1, "npc") + grid::unit(0.5, "lines"),
                           y=grid::unit(0.8, "npc"), rot=90,
                           gp=grid::gpar(col="orchid"))
       },
       pos = "CheckExEnv")
setHook("plot.new",     get("base_plot_hook", pos = "CheckExEnv"))
setHook("persp",        get("base_plot_hook", pos = "CheckExEnv"))
setHook("grid.newpage", get("grid_plot_hook", pos = "CheckExEnv"))
assign("cleanEx",
       function(env = .GlobalEnv) {
	   rm(list = ls(envir = env, all.names = TRUE), envir = env)
           RNGkind("default", "default")
	   set.seed(1)
   	   options(warn = 1)
	   .CheckExEnv <- as.environment("CheckExEnv")
	   delayedAssign("T", stop("T used instead of TRUE"),
		  assign.env = .CheckExEnv)
	   delayedAssign("F", stop("F used instead of FALSE"),
		  assign.env = .CheckExEnv)
	   sch <- search()
	   newitems <- sch[! sch %in% .oldSearch]
	   for(item in rev(newitems))
               eval(substitute(detach(item), list(item=item)))
	   missitems <- .oldSearch[! .oldSearch %in% sch]
	   if(length(missitems))
	       warning("items ", paste(missitems, collapse=", "),
		       " have been removed from the search path")
       },
       pos = "CheckExEnv")
assign("ptime", proc.time(), pos = "CheckExEnv")
## at least one package changes these via ps.options(), so do this
## before loading the package.
## Use postscript as incomplete files may be viewable, unlike PDF.
## Choose a size that is close to on-screen devices, fix paper
ps.options(width = 7, height = 7, paper = "a4", reset = TRUE)
grDevices::postscript("gbm-Ex.ps")
		      
assign("par.postscript", graphics::par(no.readonly = TRUE), pos = "CheckExEnv")
options(contrasts = c(unordered = "contr.treatment", ordered = "contr.poly"))
options(warn = 1)    
library('gbm')

assign(".oldSearch", search(), pos = 'CheckExEnv')
assign(".oldNS", loadedNamespaces(), pos = 'CheckExEnv')
cleanEx(); nameEx("calibrate.plot")
### * calibrate.plot

flush(stderr()); flush(stdout())

### Name: calibrate.plot
### Title: Calibration plot
### Aliases: calibrate.plot
### Keywords: hplot

### ** Examples

library(rpart)
data(kyphosis)
y <- as.numeric(kyphosis$Kyphosis)-1
x <- kyphosis$Age
glm1 <- glm(y~poly(x,2),family=binomial)
p <- predict(glm1,type="response")
calibrate.plot(y, p, xlim=c(0,0.6), ylim=c(0,0.6))



cleanEx(); nameEx("gbm")
### * gbm

flush(stderr()); flush(stdout())

### Name: gbm
### Title: Generalized Boosted Regression Modeling
### Aliases: gbm gbm.more gbm.fit
### Keywords: models nonlinear survival nonparametric tree

### ** Examples

# A least squares regression example
# create some data

N <- 1000
X1 <- runif(N)
X2 <- 2*runif(N)
X3 <- ordered(sample(letters[1:4],N,replace=TRUE),levels=letters[4:1])
X4 <- factor(sample(letters[1:6],N,replace=TRUE))
X5 <- factor(sample(letters[1:3],N,replace=TRUE))
X6 <- 3*runif(N)
mu <- c(-1,0,1,2)[as.numeric(X3)]

SNR <- 10 # signal-to-noise ratio
Y <- X1**1.5 + 2 * (X2**.5) + mu
sigma <- sqrt(var(Y)/SNR)
Y <- Y + rnorm(N,0,sigma)

# introduce some missing values
X1[sample(1:N,size=500)] <- NA
X4[sample(1:N,size=300)] <- NA

data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)

# fit initial model
gbm1 <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
    data=data,                   # dataset
    var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease,
                                 # +1: monotone increase,
                                 #  0: no monotone restrictions
    distribution="gaussian",     # bernoulli, adaboost, gaussian,
                                 # poisson, coxph, and quantile available
    n.trees=3000,                # number of trees
    shrinkage=0.005,             # shrinkage or learning rate,
                                 # 0.001 to 0.1 usually work
    interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc.
    bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
    train.fraction = 0.5,        # fraction of data for training,
                                 # first train.fraction*N used for training
    n.minobsinnode = 10,         # minimum total weight needed in each node
    cv.folds = 5,                # do 5-fold cross-validation
    keep.data=TRUE,              # keep a copy of the dataset with the object
    verbose=TRUE)                # print out progress

# check performance using an out-of-bag estimator
# OOB underestimates the optimal number of iterations
best.iter <- gbm.perf(gbm1,method="OOB")
print(best.iter)

# check performance using a 50% heldout test set
best.iter <- gbm.perf(gbm1,method="test")
print(best.iter)

# check performance using 5-fold cross-validation
best.iter <- gbm.perf(gbm1,method="cv")
print(best.iter)

# plot the performance
# plot variable influence
summary(gbm1,n.trees=1)         # based on the first tree
summary(gbm1,n.trees=best.iter) # based on the estimated best number of trees

# compactly print the first and last trees for curiosity
print(pretty.gbm.tree(gbm1,1))
print(pretty.gbm.tree(gbm1,gbm1$n.trees))

# make some new data
N <- 1000
X1 <- runif(N)
X2 <- 2*runif(N)
X3 <- ordered(sample(letters[1:4],N,replace=TRUE))
X4 <- factor(sample(letters[1:6],N,replace=TRUE))
X5 <- factor(sample(letters[1:3],N,replace=TRUE))
X6 <- 3*runif(N)
mu <- c(-1,0,1,2)[as.numeric(X3)]

Y <- X1**1.5 + 2 * (X2**.5) + mu + rnorm(N,0,sigma)

data2 <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)

# predict on the new data using "best" number of trees
# f.predict generally will be on the canonical scale (logit,log,etc.)
f.predict <- predict.gbm(gbm1,data2,best.iter)

# least squares error
print(sum((data2$Y-f.predict)^2))

# create marginal plots
# plot variable X1,X2,X3 after "best" iterations
par(mfrow=c(1,3))
plot.gbm(gbm1,1,best.iter)
plot.gbm(gbm1,2,best.iter)
plot.gbm(gbm1,3,best.iter)
par(mfrow=c(1,1))
# contour plot of variables 1 and 2 after "best" iterations
plot.gbm(gbm1,1:2,best.iter)
# lattice plot of variables 2 and 3
plot.gbm(gbm1,2:3,best.iter)
# lattice plot of variables 3 and 4
plot.gbm(gbm1,3:4,best.iter)

# 3-way plots
plot.gbm(gbm1,c(1,2,6),best.iter,cont=20)
plot.gbm(gbm1,1:3,best.iter)
plot.gbm(gbm1,2:4,best.iter)
plot.gbm(gbm1,3:5,best.iter)

# do another 100 iterations
gbm2 <- gbm.more(gbm1,100,
                 verbose=FALSE) # stop printing detailed progress



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx(); nameEx("print.gbm")
### * print.gbm

flush(stderr()); flush(stdout())

### Name: print.gbm
### Title: ~~function to do ... ~~
### Aliases: print.gbm
### Keywords: models nonlinear survival nonparametric tree

### ** Examples

library( gbm )
data( iris )
iris.mod <- gbm( Species ~ ., distribution="kclass", data=iris,
                 n.trees=2000, shrinkage=.01, cv.folds=5 )
iris.mod
data( lung )
lung.mod <- gbm( Surv(time, status) ~ ., distribution="coxph", data=lung,
                 n.trees=2000, shrinkage=.01, cv.folds=5 )
lung.mod



cleanEx(); nameEx("quantile.rug")
### * quantile.rug

flush(stderr()); flush(stdout())

### Name: quantile.rug
### Title: Quantile rug plot
### Aliases: quantile.rug
### Keywords: aplot

### ** Examples

x <- rnorm(100)
y <- rnorm(100)
plot(x,y)
quantile.rug(x)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
