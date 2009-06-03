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



