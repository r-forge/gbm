# Create some output to test the new print function and
# to be comparable with any future changes to gbm.

library( MASS )

data( birthwt )
data( VA )
data( iris ) 
data( fgl )
data( cpus )
data( housing )

set.seed( 20090417 )


bMod <- gbm( low ~ ., data=birthwt,
             n.tree=1000, shrinkage=.01, cv.folds=5 )
bMod

bwt <- birthwt
bwt <- bwt[ sample( 1:nrow( bwt ) ),]
aMod <- gbm( low ~ ., data=bwt, distribution="adaboost",
             n.trees=1000, shrinkage=.01, cv.folds=10,
        train.fraction=.9 )
aMod

cMod <- gbm( Surv( stime, status ) ~ treat + age + Karn + diag.time + cell + prior,
             data = VA, n.tree = 1000, shrinkage=.1, cv.folds = 5 )
cMod

kMod <- gbm( Species ~ . , data=iris , n.tree=1000, shrinkage=.1,
             cv.folds=5, train.fraction=.9 )
kMod

kMod2 <- gbm( type ~ ., data=fgl, n.tree=1000, shrinkage=.01,
              cv.folds=5 )
kMod2

mycpus <- cpus
mycpus <- mycpus[, -1 ]
gMod <- gbm( log( perf ) ~ ., data = mycpus, distribution="gaussian",
             cv.folds=5, n.trees=1000, shrinkage=.01 )
gMod

biMod <- gbm( log(perf) ~ ., data=mycpus,
              cv.folds=5, n.trees=1000, shrinkage=.01 )
biMod

tMod <- gbm( log(perf) ~ ., data=mycpus, distribution="tdist",
             cv.folds=5, n.trees=1000, shrinkage=.01, 
        interaction.depth= 3)
tMod

lMod <- gbm( log(perf) ~ ., data=mycpus, distribution="laplace",
             cv.folds=5, n.trees=1000, shrinkage=.01, 
        interaction.depth= 3)
lMod

qMod <- gbm( log(perf) ~ ., data=mycpus,
             distribution=list(name="quantile", alpha=.7 ),
             cv.folds=5, n.trees=1000, shrinkage=.01, 
        interaction.depth= 3)
qMod

pMod <- gbm( Freq ~ ., data=housing , distribution="poisson",
             n.trees=1000, cv.folds=5 , shrinkage=.01,
        interaction.depth = 3)
pMod
