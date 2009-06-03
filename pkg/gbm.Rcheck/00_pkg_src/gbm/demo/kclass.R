library( gbm, lib.loc = "~/R/library" )
data( iris )
source( "~/R/dev/gbm/R/gbm.R" )

set.seed( 20090415 )



mod <- gbm( Species ~ ., data= iris, distribution="kclass",
            n.tree = 5000, shrinkage=.001 , cv.folds=2,
	    bag.fraction=.8, interaction.depth=3 )
gbm.perf( mod, method="cv" )
mod


riris <- iris
riris$setosa <- iris$Species == "setosa"
riris <- riris[, names( riris ) != "Species ]

junk <- gbm( setosa ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width,
             data=riris, distribution="adaboost" ,
	     n.tree = 5000, shrinkage=.001, cv.folds=5,
	     bag.fraction=.8, interaction.depth=3 )

