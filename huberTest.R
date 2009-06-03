library(MASS)
library(gbm)
x <- mvrnorm( 500, mu=rep(1,5), Sigma=diag(rep(1,5)))
y <- ifelse( x[,1]^2 + x[,2]^2 + x[,3]^3 > 4, 1, 0 )
mc <- sample( 1:500 , size=50 )
#y[mc] <- abs( y-1 )
d <- data.frame(y=y, x=x)

bmod <- gbm( y ~ ., data=d, dist="bernoulli", n.tree=3000,
             shrinkage=.01, cv.folds=3 )
hmod <- gbm( y ~ ., data=d, dist="huberized",
             n.tree=2000, shrinkage=.005, cv.folds=3 )
amod <- gbm( y ~ ., data=d, dist="adaboost",
             n.tree=2000, shrinkage=.01, cv.folds=3 )

par(mfrow=c( 2, 2 ) )
hbest <- gbm.perf(hmod)
bbest <- gbm.perf(bmod)
abest <- gbm.perf(amod)

x <- mvrnorm( 5000, mu=rep(1,5), Sigma=diag(rep(1,5)))
y <- ifelse( x[,1]^2 + x[,2]^2 + x[,3]^3 > 4, 1, 0 )
nd <- data.frame(y=y, x=x)


hpred <- predict( hmod, newdata=nd, type="response" )
bpred <- predict( bmod, newdata=nd, type="response" )

