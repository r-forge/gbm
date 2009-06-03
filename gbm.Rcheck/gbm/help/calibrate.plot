calibrate.plot              package:gbm              R Documentation

_C_a_l_i_b_r_a_t_i_o_n _p_l_o_t

_D_e_s_c_r_i_p_t_i_o_n:

     An experimental diagnostic tool that plots the fitted values
     versus the actual average values. Currently developed for only
     'distribution="bernoulli"'.

_U_s_a_g_e:

     calibrate.plot(y,p,
                    distribution="bernoulli",
                    replace=TRUE,
                    line.par=list(col="black"),
                    shade.col="lightyellow",
                    shade.density=NULL,
                    rug.par=list(side=1),
                    xlab="Predicted value",
                    ylab="Observed average",
                    xlim=NULL,ylim=NULL,
                    knots=NULL,df=6,
                    ...)

_A_r_g_u_m_e_n_t_s:

       y: the outcome 0-1 variable 

       p: the predictions estimating E(y|x) 

distribution: the loss function used in creating 'p'. 'bernoulli' and
          'poisson' are currently the only special options. All others
          default to squared error assuming 'gaussian'

 replace: determines whether this plot will replace or overlay the
          current plot. 'replace=FALSE' is useful for comparing the
          calibration of several methods

line.par: graphics parameters for the line 

shade.col: color for shading the 2 SE region. 'shade.col=NA' implies no
          2 SE region

shade.density: the 'density' parameter for 'polygon'

 rug.par: graphics parameters passed to 'rug'

    xlab: x-axis label corresponding to the predicted values

    ylab: y-axis label corresponding to the observed average

xlim,ylim: x and y-axis limits. If not specified te function will
          select limits

knots,df: these parameters are passed directly to  'ns' for
          constructing a natural spline  smoother for the calibration
          curve

     ...: other graphics parameters passed on to the plot function 

_D_e_t_a_i_l_s:

     Uses natural splines to estimate E(y|p). Well-calibrated
     predictions imply that E(y|p) = p. The plot also includes a
     pointwise 95 band.

_V_a_l_u_e:

     'calibrate.plot' returns no values.

_A_u_t_h_o_r(_s):

     Greg Ridgeway gregr@rand.org

_R_e_f_e_r_e_n_c_e_s:

     J.F. Yates (1982). "External correspondence: decomposition of the
     mean probability score," Organisational Behaviour and Human
     Performance 30:132-156.

     D.J. Spiegelhalter (1986). "Probabilistic Prediction in Patient
     Management and Clinical Trials," Statistics in Medicine 5:421-433.

_E_x_a_m_p_l_e_s:

     library(rpart)
     data(kyphosis)
     y <- as.numeric(kyphosis$Kyphosis)-1
     x <- kyphosis$Age
     glm1 <- glm(y~poly(x,2),family=binomial)
     p <- predict(glm1,type="response")
     calibrate.plot(y, p, xlim=c(0,0.6), ylim=c(0,0.6))

