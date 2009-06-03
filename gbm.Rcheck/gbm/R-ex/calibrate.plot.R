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



