var documenterSearchIndex = {"docs":
[{"location":"dataAssimilation/4dvar/#D-Var","page":"4d var","title":"4D-Var","text":"","category":"section"},{"location":"dataAssimilation/4dvar/","page":"4d var","title":"4d var","text":"The 3D-Var algorithm attempts to optimize a cost function to obtain the ideal analysis for each point where we have observation data. This can become computationally expensive as we require model evaluations and an optimization routine for every observation point. An alternative approach is to simultaneously optimize accross all observations in order to obtain the ideal initial condition that acheive the best model fit. This approach is similar to sensativity analysis which seeks to fit a model's parameters to data. ","category":"page"},{"location":"dataAssimilation/4dvar/","page":"4d var","title":"4d var","text":"To begin, we construct the 4d-var cost function","category":"page"},{"location":"dataAssimilation/4dvar/","page":"4d var","title":"4d var","text":"beginaligned\n    J(u_0) = frac12left( u_0 - u_0^(b) right)^TB^-1left( u_0 - u_0^(b) right) + frac12sum_kleft(w_k - h(u_k) right)^TR_k^-1left(w_k - h(u_k) right)  \n           = J_b(u_0) + J_m(u_0)\nendaligned","category":"page"},{"location":"dataAssimilation/4dvar/","page":"4d var","title":"4d var","text":"The first term is usefull if we already have an initial guess u_0^(b) for the inital condition in mind. If we do not have one, we may ommit this term. ","category":"page"},{"location":"dataAssimilation/4dvar/","page":"4d var","title":"4d var","text":"As before, we now want to optimize this cost function. To do so, we first observe that ","category":"page"},{"location":"dataAssimilation/4dvar/","page":"4d var","title":"4d var","text":"    u_k = mathcalM^(k)(u_0 theta)","category":"page"},{"location":"dataAssimilation/4dvar/","page":"4d var","title":"4d var","text":"It is easy to obtain the gradient of J_0 so we shall focus on the second term. We find that ","category":"page"},{"location":"dataAssimilation/4dvar/","page":"4d var","title":"4d var","text":"beginaligned\n    nabla_u_0J_m = nabla_u_0Big sum_k frac12  left(w_k - h(u_k) right)^TR_k^-1left(w_k - h(u_k) right) Big\n                    = - sum_k leftdfracpartial partial u_0hleft(mathcalM^(k-1)(u_0)right) right^T R_k^-1left(w_k - h(u_k) right)\n                    = - sum_k leftD_h(u_k)D_M(u_k-1)D_M(u_k-2)cdots D_M(u_0) right^T R_k^-1left(w_k - h(u_k) right)\n                    = - sum_k leftD_M^T(u_0)D_M^T(u_1)cdots D_M^T(u_k-1)D_h^T(u_k) right R_k^-1left(w_k - h(u_k) right)\nendaligned","category":"page"},{"location":"dataAssimilation/4dvar/","page":"4d var","title":"4d var","text":"Given that we can now obtain the gradient of the cost function, the procedure is nearly identical to 3d-var: ","category":"page"},{"location":"dataAssimilation/4dvar/","page":"4d var","title":"4d var","text":"Integrate your model forward to obtain u_k\nEvaluate each of the D_M^T(u_k-10) and D_h(u_k). \nUsing these values, compute nabla J_m(u)\nSet u_0^(new) = u_0^(prev) - eta nabla J(u_0^(prev))\nStop when lvert u_0^(new) - u_0^(prev) rvert converges to your desired tolerance. ","category":"page"},{"location":"dataAssimilation/4dvar/","page":"4d var","title":"4d var","text":"You can of course substitute another optimzation scheme after step 3. ","category":"page"},{"location":"dataAssimilation/sensativity/#Sensativity-Analysis-for-Differential-Equations","page":"Sensativity Analysis","title":"Sensativity Analysis for Differential Equations","text":"","category":"section"},{"location":"dataAssimilation/sensativity/","page":"Sensativity Analysis","title":"Sensativity Analysis","text":"Provided some model for a physical system in the form of a set of differential equations, a natural question is: How can we select the parameters for our model in order to get the best possible fit to some experimental data. Similarly, one may wonder what would happen to the prediction of your model if you were to slightly change the values of some parameters. In other words, how sensative is the output of our model to your choice of parameter values? ","category":"page"},{"location":"dataAssimilation/sensativity/","page":"Sensativity Analysis","title":"Sensativity Analysis","text":"In the most general sense, we may frame the problem as follows. Suppose we have a model of the form ","category":"page"},{"location":"dataAssimilation/sensativity/","page":"Sensativity Analysis","title":"Sensativity Analysis","text":"dfracdudt = f(uttheta)","category":"page"},{"location":"dataAssimilation/sensativity/","page":"Sensativity Analysis","title":"Sensativity Analysis","text":"Given this model, our goal is to optimize a cost function ","category":"page"},{"location":"dataAssimilation/sensativity/","page":"Sensativity Analysis","title":"Sensativity Analysis","text":"J(u theta) = int_0^T g(utheta)dt","category":"page"},{"location":"dataAssimilation/sensativity/","page":"Sensativity Analysis","title":"Sensativity Analysis","text":"where g(utheta) is usually taken to be some quadratic form. ","category":"page"},{"location":"dataAssimilation/sensativity/","page":"Sensativity Analysis","title":"Sensativity Analysis","text":"As an example, we might consider g(u theta) = (u(t)-w(t))^T(u(t)-w(t)) where w(t) denotes the vector of observations at time t. ","category":"page"},{"location":"dataAssimilation/sensativity/","page":"Sensativity Analysis","title":"Sensativity Analysis","text":"Our goal then is to find out how J depends on the parameters theta, in other words, to find partial J  partial theta. To do this, we will use the method of Lagrange multipliers to generate a so called adjoint equation that enables us to find this derivative in a way that minimizes computational cost. As always, this method begins by adding a term that evaluates to 0 into our cost function: ","category":"page"},{"location":"dataAssimilation/sensativity/","page":"Sensativity Analysis","title":"Sensativity Analysis","text":"    mathcalL = int_0^T left g(utheta) + lambda^T(t)left(f-dfracdudtright) right dt","category":"page"},{"location":"dataAssimilation/sensativity/","page":"Sensativity Analysis","title":"Sensativity Analysis","text":"From this, we find ","category":"page"},{"location":"dataAssimilation/sensativity/","page":"Sensativity Analysis","title":"Sensativity Analysis","text":"beginaligned\n    dfracpartial mathcalLpartial theta = int_0^Tleft fracpartial gpartial theta + fracpartial gufracpartial upartial theta + lambda^T(t)left( fracpartial fpartial theta + fracpartial fpartial ufracpartial upartial theta - fracddtfracpartial upartial theta right)rightdt  \n    = int_0^T left fracpartial gpartial theta + lambda^T(t)fracpartial fpartial theta + left( fracpartial gpartial u + lambda^T(t)fracpartial fpartial u - lambda^T(t)fracddt right)fracpartial upartial theta rightdt\nendaligned","category":"page"},{"location":"dataAssimilation/sensativity/","page":"Sensativity Analysis","title":"Sensativity Analysis","text":"This reorganization is nice because the term partial upartial theta is the one thats hard to compute. Therefore, if we can make the terms in the paretheses evaluate to 0, we will be able to remove this pesky term. Let's use integration by parts to further rearrange by moving the ddt.","category":"page"},{"location":"dataAssimilation/sensativity/","page":"Sensativity Analysis","title":"Sensativity Analysis","text":"beginaligned\n    int_0^T-lambda^T(t)fracddtfracpartial upartial theta dt = left-lambda^T(t)fracpartial upartial theta right_0^T + int_0^T fracdlambda^T(t)dtfracpartial upartial thetadt  \n    = lambda^T(0)fracpartial u_0partial theta - lambda^T(T)fracpartial u(T)partial theta + int_0^T left fracdlambdadt right^Tfracpartial upartial thetadt\nendaligned","category":"page"},{"location":"dataAssimilation/sensativity/","page":"Sensativity Analysis","title":"Sensativity Analysis","text":"so that plugging this back into our expression for partial mathcalLpartial theta, we obtain","category":"page"},{"location":"dataAssimilation/sensativity/","page":"Sensativity Analysis","title":"Sensativity Analysis","text":"fracpartial mathcalLpartial theta = int_0^T left fracpartial gpartial theta + lambda^Tfracpartial fpartial theta + left( fracpartial gpartial u + lambda^Tfracpartial fpartial u + leftfracdlambdadtright^T right)fracpartial upartial thetarightdt + lambda^T(0)fracpartial u_0partial theta - lambda^T(T)fracpartial u(T)partial theta","category":"page"},{"location":"dataAssimilation/sensativity/","page":"Sensativity Analysis","title":"Sensativity Analysis","text":"Thus, forcing the nasty terms to dissappear is equivalent find the lambda(t) subject to the differential equations ","category":"page"},{"location":"dataAssimilation/sensativity/","page":"Sensativity Analysis","title":"Sensativity Analysis","text":"beginaligned\n   fracpartial gpartial u + lambda^T(t)fracpartial fpartial u + fracdlambda^T(t)dt = 0  \n   lambda^T(T) = 0\nendaligned","category":"page"},{"location":"dataAssimilation/sensativity/","page":"Sensativity Analysis","title":"Sensativity Analysis","text":"or by taking the transpose: ","category":"page"},{"location":"dataAssimilation/sensativity/","page":"Sensativity Analysis","title":"Sensativity Analysis","text":"beginaligned\n    fracddtlambda = - left fracpartial gpartial u right^T - left fracpartial fpartial u right^Tlambda   \n    lambda(T) = 0\nendaligned","category":"page"},{"location":"dataAssimilation/sensativity/#Summary","page":"Sensativity Analysis","title":"Summary","text":"","category":"section"},{"location":"dataAssimilation/sensativity/","page":"Sensativity Analysis","title":"Sensativity Analysis","text":"To find the sensativities partial Jpartial theta, we perform the following: ","category":"page"},{"location":"dataAssimilation/sensativity/","page":"Sensativity Analysis","title":"Sensativity Analysis","text":"Integrate the model dudt = f(uttheta) forward to obtain u(t).\nIntegrate the adjoint model dlambdadt = -(partial f partial u)^Tlambda - (partial g  partial u)^T backwards in time from T to 0 to obtain lambda(t).\nEvaluate partial J  partial theta = int_0^Tleft( partial g partial theta + lambda^T partial fpartial thetaright)dt + lambda^T(0)partial u_0partial theta","category":"page"},{"location":"dataAssimilation/kalman/#Kalman-Filtering","page":"Kalman Filter","title":"Kalman Filtering","text":"","category":"section"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"Given some model for the error covariance matrices Q_k and R_k, we would like a method that propagates both our model and the errors forward. This way we may guarantee that the accuracy of our analysis doesn't come at the cost of higher uncertainty. ","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"The original implementation of the Kalman filter was for strictly linear systems. We will first develop the analysis for this simplified case adn then will generalize to the Extended Kalman Filter (EKF) that can handle fully nonlinear situations.","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"In the linear case, our system may be written as ","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"beginaligned\n    u_k+1^(t) = M_ku_k^(t) + xi_k+1^(p)  \n    w_k = H_ku_k^(t) + xi_k^(m)\nendaligned","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"where M_k and H_k are now matrices defining the linear problem. ","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"The goal of the Kalman filter is to derive the analysis u^(a) which optimizes the trace of the analysis error covariance matrix (i.e. sum of squared errors): ","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"mathrmTrleft( P_kright) = E(u_k^(t)-u_k^(a))^T(u_k^(t)-u_k^(a))","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"Finding the analysis consists of two steps: the forecast step and the assimilation step.","category":"page"},{"location":"dataAssimilation/kalman/#Forecast-Step","page":"Kalman Filter","title":"Forecast Step","text":"","category":"section"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"Assume we have the analysis at time t_k denoted u_k^(a). Then the forecast for time t_k+1 is","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"    u_k+1^(b) = M_ku_k^(a)","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"The background error is therefore ","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"beginaligned\n    xi_k+1^(b) = u_k+1^(t) - u_k+1^(b)  \n    = M_ku_k^(t)+xi_k+1^(p) - M_ku_k^(a)  \n    = M_kleft(u_k^(t)-u_k^(a) right) + xi_k+1^(p)  \n    = M_kxi_k^(a) + xi_k+1^(p)\nendaligned","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"We may now evaluate the covariance matrix of our background estimate as: ","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"beginaligned\n    B_k+1 = Exi_k+1^(b)(xi_k+1^(b))^T  \n    = Eleftleft(M_kxi_k^(a) + xi_k+1^p right) left(M_kxi_k^(a) + xi_k+1^p right)^T right  \nendaligned","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"If we presume that Exi_k^(b)(xi_k+1^(p))^T = 0, then the cross terms vanish and we are left with ","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"boxedB_k+1 = M_kP_kM_k^T + Q_k+1","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"Thus we now have the background (i.e forecast) estimate of the state at t_k+1 and its covariance matrix. Given a measurement w_k+1 at the same time with covariance matrix R_k+1, then we may now perform the assimilation step where we fuse the two sources of information to obtain u_k+1^(a) and P_k+1.","category":"page"},{"location":"dataAssimilation/kalman/#Data-Assimilation-Step","page":"Kalman Filter","title":"Data Assimilation Step","text":"","category":"section"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"Let's suppose that the analysis has the form ","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"u_k+1^(a) = nu + K_k+1w_k+1","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"for some vector nuinR^n and matrix K_k+1inR^mtimes n. In a perfect world, we would have Eu_k^(t)-u_k^(a) = 0. Therefore, ","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"beginaligned\n    0 = Eu_k^(t) - u_k^(a)  \n    = E(u_k^(b) + xi_k^(b)) - (nu + K_kw_k)  \n    = E(u_k^(b) + xi_k^(b)) - (nu + K_kH_ku_k^(t) + K_kxi_k^(m))  \n    = Eu_k^(b) + Exi_k^(b) - Enu -K_kH_kEu_k^(t) - K_kExi_k^(m)\n    = u_k^(b) + 0 - nu - K_kH_ku_k^(b) - 0  \n    = u_k^(b) - nu - K_kH_ku_k^(b)  \n    Rightarrow nu = u_k^(b) - K_kH_ku_k^(b)\nendaligned","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"which we now substitute to obtain ","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"boxedu_k^(a) = u_k^(b) + K_k(w_k - H_ku_k^(b))","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"Now that we know the form for the analysis we may derive the optimal matrix K_k by optimization of P_k. We have","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"beginaligned\n\txi_k^(a) = u_k^(t) - u_k^(a)  \n                = M_k-1u_k-1^(t) + xi_k^(p) - u_k^(b) - K_kleft(w_k - H_ku_k^(b) right)  \n                = M_k-1u_k-1^(t) + xi_k^(p) - M_k-1u_k-1^(a) - K_kleft(H_ku_k^(t) + xi_k^(m) - H_ku_k^(b) right)  \n                = M_k-1u_k-1^(t) + xi_k^(p) - M_k-1u_k-1^(a) - K_kH_ku_k^(t) - K_kxi_k^(m) + K_kH_ku_k^(b)  \n                = M_k-1u_k-1^(t) + xi_k^(p) - M_k-1u_k-1^(a) - K_kH_ku_k^(t) - K_kxi_k^(m) + K_kH_ku_k^(b)  \n                = Big M_k-1(xi_k-1^(a)+u_k-1^(a)) + xi_k^(p) - M_k-1u_k-1^(a) Big - K_kH_ku_k^(t) - K_kxi_k^(m) + K_kH_ku_k^(b)  \n                = Big M_k-1xi_k-1^(a) + xi_k^(p) Big - K_kH_ku_k^(t) + K_kH_ku_k^(b) - K_kxi_k^(m) \n                = M_k-1xi_k-1^(a) + xi_k^(p) - K_kH_k(M_k-1u_k-1^(t) + xi_k^(b)) + K_kH_ku_k^(b) - K_kxi_k^(m) \n                = M_k-1xi_k-1^(a) + xi_k^(p) - K_kH_kM_k-1(xi_k-1^(a) + u_k-1^a) - K_kH_kxi_k^(b) + K_kH_ku_k^(b) - K_kxi_k^(m) \n                = M_k-1xi_k-1^(a) + xi_k^(p) - K_kH_kM_k-1(xi_k-1^(a) + u_k-1^a) - K_kH_kxi_k^(b) + K_kH_kM_k-1u_k-1^(a) - K_kxi_k^(m) \n                = M_k-1xi_k-1^(a) + xi_k^(p) - K_kH_kM_k-1xi_k-1^(a) - K_kH_kxi_k^(b) - K_kxi_k^(m) \n                = big(I-K_kH_k big)(M_k-1xi_k-1^(a) - xi_k^p) - K_kxi_k^(m)\nendaligned","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"and therefore the covariance matrix is ","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"beginaligned\n    P_k = Exi_k^(a)(xi_k^(a))^T  \n        = EBigleft(big(I-K_kH_k big)(M_k-1xi_k-1^(a) - xi_k^p) - K_kxi_k^(m) right) left(big(I-K_kH_k big)(M_k-1xi_k-1^(a) - xi_k^p) - K_kxi_k^(m) right)^T Big  \n        = big(I-K_kH_k big)M_k-1Exi_k-1^(a)(xi_k-1^(a))^TM_k-1^Tbig(I-K_kH_k big)^T + big(I-K_kH_k big)Exi_k^(p)(xi_k^(p))^Tbig(I-K_kH_k big)^T - K_kExi_k^(m)(xi_k^(m))^TK_k^T  \n        = big(I-K_kH_k big)B_kbig(I-K_kH_k big)^T - K_kR_kK_k^T\nendaligned","category":"page"},{"location":"dataAssimilation/kalman/#Deriving-K_k","page":"Kalman Filter","title":"Deriving K_k","text":"","category":"section"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"The Kalman filter is defined at that K_k which which minimizes the sum of squared analysis errors, i.e. the trace of the analysis error covariance matrix. The following identies will be useful: ","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"beginaligned\n    mathopnabla_Atexttr(AB) = B^T  \n    mathopnabla_Atexttr(BA^T) = B  \n    mathopnabla_Atexttr(ABA^T) = AB^T + AB   \nendaligned","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"from which we obtain ","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"beginaligned\n    0 = mathopnabla_K_ktexttr(P_k)  \n      = mathopnabla_K_kBig B_k -B_kH_k^TK_k^T - K_kH_kB_k  + K_kH_kB_kH_k^TB_k^T - K_kR_kK_k^T Big  \n      = -B_kH_k^T - (H_kB_k)^T + K_kleftH_kB_kH_k^T + (H_kB_kH_k^T)^T - R_k+R_k^Tright \n      = -2B_kH_k^T + 2K_kleft(H_kB_kH_k^2 - R_k right)  \n  Rightarrow K_k = B_kH_k^TBig H_kB_kH_k^T - R_k Big^-1\nendaligned","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"we now substitute this result to obtain a simplified form for P_k. ","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"beginaligned\n    P_k = left(I - K_kH_k right)B_kleft(I - K_kH_k right)^T + K_kR_kK_k^T  \n        = left(I - K_kH_k right)B_k - left(I - K_kH_k right)B_kleft(K_kH_k right)^T + K_kR_kK_k^T \n        = left(I - K_kH_k right)B_k -left left(I - K_kH_k right)B_kleft(K_kH_k right)^T + K_kR_kK_k^T right \n        = left(I - K_kH_k right)B_k -left left(I - K_kH_k right)B_kH_k^TK_k^T + K_kR_kK_k^T right \n        = left(I - K_kH_k right)B_k -left left(I - K_kH_k right)B_kH_k^T + K_kR_k rightK_k^T \n        = left(I - K_kH_k right)B_k -left B_kH_k^T - K_kleft( H_kB_kH_k^T + R_k right)  rightK_k^T \n        = left(I - K_kH_k right)B_k -left B_kH_k^T - B_kH_k^T rightK_k^T \n        = left(I - K_kH_k right)B_k \nendaligned","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"NOTE: we have used the fact that covariance matrices are symmetric. ","category":"page"},{"location":"dataAssimilation/kalman/#Summary","page":"Kalman Filter","title":"Summary","text":"","category":"section"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"Let's summarize the whole process. We have ","category":"page"},{"location":"dataAssimilation/kalman/#.-Initialization","page":"Kalman Filter","title":"0. Initialization","text":"","category":"section"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"We must set the system to some initial condition. This means we must define u_0^a and P_0. We must also come up with a model for the process noise covariance Q_k and measurement error covariance R_k.","category":"page"},{"location":"dataAssimilation/kalman/#.-Forecast-Step","page":"Kalman Filter","title":"1. Forecast Step","text":"","category":"section"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"beginaligned\n    u_k^(b) = M_k-1u_k-1^(a)  \n    B_k = M_k-1P_k-1M_k-1^T + Q_k\nendaligned","category":"page"},{"location":"dataAssimilation/kalman/#.-Assimilation-Step","page":"Kalman Filter","title":"2. Assimilation Step","text":"","category":"section"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"beginaligned\n    K_k = B_kH_k^TBig H_kB_kH_k^T - R_k Big^-1   \n    u_k^(a) = u_k^(b) + K_k(w_k - H_ku_k^(b)) \n    P_k = left(I - K_kH_k right)B_k \nendaligned","category":"page"},{"location":"dataAssimilation/ekf/#Extended-Kalman-Filter","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"","category":"section"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"Given the nonlinear nature of many scientific models it is desirable to extend the Kalman Filter to be able to handle nonlinear models f(cdot) (and by extension, their update function mathcalM(cdot)), and nonlinear observation functions h(cdot). This can be accomplished so long as these functions are sufficiently smooth (C^1 to be precises) so as to admit valid Taylor approximations to first order. That is, ","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"beginaligned\n    mathcalM(u_k) approx mathcalM(u_k^(a)) + D_M(u_k^(a))xi_k^(a)  h(u_k) approx h(u_k^(b)) + D_h(u_k^(b))xi_k^(b)  \n    D_M = leftdfracpartial mathcalM_ipartial u_j right  D_h = left dfracpartial h_ipartial u_jright\nendaligned","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"where mathcalM_i and h_i  denote the ith component functions of mathcalM and h. ","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"Using these substitutions for the previously linear functions M_k and H_k, we may follow the same derivation to obtain the following procedure. ","category":"page"},{"location":"dataAssimilation/ekf/#.-Initialization","page":"Extended Kalman Filter","title":"0. Initialization","text":"","category":"section"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"To begin we must choose values for u_0^(a) and P_0. We must also provide models for Q_k and R_k. ","category":"page"},{"location":"dataAssimilation/ekf/#.-Forecast-Step","page":"Extended Kalman Filter","title":"1. Forecast Step","text":"","category":"section"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"beginaligned\n    u_k^(b) = mathcalM(u_k-1^(a))  \n    B_k = D_M(u_k-1^(a))P_k-1D_M^T(u_k-1^(a)) + Q_k\nendaligned","category":"page"},{"location":"dataAssimilation/ekf/#.-Assimilation-Step","page":"Extended Kalman Filter","title":"2. Assimilation Step","text":"","category":"section"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"beginaligned\n    K_k = B_kD_M^T(u_k^(b))left D_h(u_k^(b))B_kD_h^T(u_k^(b)) + R_k right^-1\n    u_k^(a) = u_k^(b) + K_k(w_k - h(u_k^(b)))  \n    P_k = left( I - K_kD_h(u_k^(b)) right)B_k\nendaligned","category":"page"},{"location":"fdocs/#Function-Documentation","page":"Function docs","title":"Function Documentation","text":"","category":"section"},{"location":"fdocs/","page":"Function docs","title":"Function docs","text":"Modules = [DataAssimilation]","category":"page"},{"location":"dataAssimilation/dataassim/#Data-Assimilation-Overview","page":"Data Assimilation Oveview","title":"Data Assimilation Overview","text":"","category":"section"},{"location":"dataAssimilation/dataassim/#Overview","page":"Data Assimilation Oveview","title":"Overview","text":"","category":"section"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"The proper application of scientific models to make real-world predictions requires that we commit ourselves to a full accounting of all possible sources of uncertainty when reporting results. Further, the explosion of big data across scientific fields provieds a plethora observational data that our models are typically unequipped to incorporate when making predictions. The field of Data Assimilation addresses this problem by providing a family of techniques engineered to combine model output together with observational data whilst enabling a complete accounting the sources of uncertainty. For chaotic systems in particular, data assimilation enables integration on long time scales that would be impossible via models alone. ","category":"page"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"In this overview, we will follow the examples from this nice paper. ","category":"page"},{"location":"dataAssimilation/dataassim/#Framing-the-Problem","page":"Data Assimilation Oveview","title":"Framing the Problem","text":"","category":"section"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"Data assimilation can be understood most generally in terms of dyscrete dynamical systems. This enables us to apply the methods to most mathematical models from gridded PDE solvers to systems of ordinary differential equations. Our goal is to find the best prediction for the system state vector u that combines our model predictions, also known as forecasts, with observational data. Model predictions are summarized via the discrete update equation: ","category":"page"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"u_k+1 = mathcalM(u_k theta)","category":"page"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"For ODE systems, mathcalM represents the time integration scheme for a system of ODEs like ","category":"page"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"dfracdudt = f(u t theta)","category":"page"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"To measure the performance of our assimilation scheme, we denote the true value of the state vector asd u^(t). The output of our model is denoted u^(b) (b subscript for background). The discrepancy between the true value and our forecast is denoted xi^(b) = u^(t) - u^(b) characterizing the extent to which our model prediction is imperfect. ","category":"page"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"The observations of our system are denoted by w_k = w(t_k). These observations do not neccessarily need to be components of the state vector u, but rather, are related to it via the observation function h. For example, one may attempt to predict sea surface temperature using data assimilation with data from satellite observations. The function h would then be the Stefan-Boltzmann law. However, real world data is noisy, which we must take into account. We write ","category":"page"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"w_k = h(u_k) + xi_k^(m)","category":"page"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"where xi_k^(m) denotes this measurement noise. ","category":"page"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"Given our model predictions u_k^(b) and observations w_k, we seek to obtain the optimal or best-possible prediction called the analysis, u^(a). This analysis will still not be perfect, so we further specify the analysis error via ","category":"page"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"xi^(a) = u^(t) - u^(a)","category":"page"},{"location":"dataAssimilation/dataassim/#Summary","page":"Data Assimilation Oveview","title":"Summary","text":"","category":"section"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"beginaligned\n    u_k^(t) in R^n textthe true state vector  \n    u_k^(b) in R^n textthe kth model forecast  \n    u_k^(a) in R^n textthe analysis  \n    w_k in R^m textthe kth observation vector  \n    xi^(b) in R^n textthe model forecast error\n    xi^(m) in R^m textthe observation noise vector \n    xi^(a) in R^n textthe analysis error\n    xi^(p) in R^n textthe process noise if we used our model on the true state\n    mathcalMR^ntoR^n textthe model update function\n    fR^ntoR^n textdifferential equation model \n    hR^ntoR^m  textobservation function\nendaligned","category":"page"},{"location":"dataAssimilation/dataassim/#Assumptions","page":"Data Assimilation Oveview","title":"Assumptions","text":"","category":"section"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"To make possible the derivation of a unique analysis u^(a), the following assumptions are in order.","category":"page"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"beginaligned\n    Exi_k^(b) = 0  Exi_k^(b)(xi_j^(b))^T = 0 text for  kneq j\n    Exi_k^(m) = 0  Exi_k^(m)(xi_j^(m))^T = 0 text for  kneq j\n    Exi_k^(b)(u_0)^T = 0  Exi_k^(m)(u_0)^T = 0\n    Exi_k^(b)xi_j^(m) = 0     \n    Eu_k^(t) = u_k^(b)  \nendaligned","category":"page"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"We also define the error covariance matrices","category":"page"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"beginaligned\n    Q_k = Exi_k^(p)(xi_k^(p))^T  \n    R_k = Exi_k^(m)(xi_k^(m))^T  \n    B_k = Exi_k^(b)(xi_k^(b))^T \nendaligned","category":"page"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"which we will use in our consideration of the final error of our analysis.","category":"page"},{"location":"dataAssimilation/3dvar/#D-Var","page":"3d var","title":"3D-Var","text":"","category":"section"},{"location":"dataAssimilation/3dvar/","page":"3d var","title":"3d var","text":"For the Kalman Filter and the EKF, we derived the optimal way to combine observation with simulation so as to minimize the trace of the analysis error covariance matrix, P_k. An alternative approach is to recast the problem as a pure optimzation problem where rather than finding a filter K_k that will add an innovation to u_k^(b) to obtain the analysis u_k^(a), we obtain the analysis by optimizing the following cost function","category":"page"},{"location":"dataAssimilation/3dvar/","page":"3d var","title":"3d var","text":"J(u) = frac12left(w - h(u) right)^TR^-1left(w - h(u) right) + frac12left(u - u^(b) right)^TB^-1frac12left(u - u^(b) right)","category":"page"},{"location":"dataAssimilation/3dvar/","page":"3d var","title":"3d var","text":"which we can justify as coming from the joint probability distribution assuming Gaussian errors","category":"page"},{"location":"dataAssimilation/3dvar/","page":"3d var","title":"3d var","text":"mathcalP(uw) = Cexpleft(- frac12left(u - u^(b) right)^TB^-1frac12left(u - u^(b) right) right)cdotexpleft(-  frac12left(w - h(u) right)^TR^-1left(w - h(u) right) right)","category":"page"},{"location":"dataAssimilation/3dvar/","page":"3d var","title":"3d var","text":"with model error covariance B and measurement error covariance R as before. This is clearly a very strong assumption.","category":"page"},{"location":"dataAssimilation/3dvar/","page":"3d var","title":"3d var","text":"To optimize J(u), we begin by taking it's gradient.","category":"page"},{"location":"dataAssimilation/3dvar/","page":"3d var","title":"3d var","text":"nabla_uJ(u) = -D_h^TR^-1(w-h(u)) + B^-1(u-u^(b))","category":"page"},{"location":"dataAssimilation/3dvar/","page":"3d var","title":"3d var","text":"Thus, finding the analysis u^(a) ammounts to solving the system ","category":"page"},{"location":"dataAssimilation/3dvar/","page":"3d var","title":"3d var","text":"a-D_h^TR^-1(w-h(u^(a))) + B^-1(u^(a)-u^(b)) = 0","category":"page"},{"location":"dataAssimilation/3dvar/","page":"3d var","title":"3d var","text":"As for Kalman filtering, let's begin with the assumption that our model and observation function are linear. ","category":"page"},{"location":"dataAssimilation/3dvar/#Linear-Case","page":"3d var","title":"Linear Case","text":"","category":"section"},{"location":"dataAssimilation/3dvar/","page":"3d var","title":"3d var","text":"Suppose that we have h(u) = Hu so that D_h(u) = H. Then, we have ","category":"page"},{"location":"dataAssimilation/3dvar/","page":"3d var","title":"3d var","text":"beginaligned\n    D_h^TR^-1(w-Hu^(a)) = B^-1(u^(a)-u^(b))  \n    D_h^TR^-1w - D_h^TR^-1Hu^(a) = B^-1u^(a) - B^-1u^(b)  \n    left(D_h^TR^-1H + B^-1 right)u^(a) = D_h^TR^-1 + B^-1u^(b) \n    left(H^TR^-1H + B^-1 right)u^(a) = H^TR^-1 + B^-1u^(b)\nendaligned","category":"page"},{"location":"dataAssimilation/3dvar/","page":"3d var","title":"3d var","text":"Thus we see that the analysis is given by ","category":"page"},{"location":"dataAssimilation/3dvar/","page":"3d var","title":"3d var","text":"    u^(a) = u^(b) + BH^Tleft( R + HB^TH right)^-1(w-Hu^(b))","category":"page"},{"location":"dataAssimilation/3dvar/","page":"3d var","title":"3d var","text":"which agrees with what we found for the Linear Kalman Filter. ","category":"page"},{"location":"dataAssimilation/3dvar/#Nonlinear-Case","page":"3d var","title":"Nonlinear Case","text":"","category":"section"},{"location":"dataAssimilation/3dvar/","page":"3d var","title":"3d var","text":"To deal with the nonlinearity, we can expand h about an initial guess u^(c) which we will later choose to be u^(b) for convenience. ","category":"page"},{"location":"dataAssimilation/3dvar/","page":"3d var","title":"3d var","text":"    h(u^(a)) approx h(u^(c)) + D_h(u^(c))Delta u","category":"page"},{"location":"dataAssimilation/3dvar/","page":"3d var","title":"3d var","text":"Using this, we have ","category":"page"},{"location":"dataAssimilation/3dvar/","page":"3d var","title":"3d var","text":"beginaligned\n    D_h^T(u^(a))R^-1(w-h(u^(a))) = B^-1(u^(a) - u^(b))  \n    D_h^T(u^(a))R^-1(w-h(u^(c))-D_h(u^(c))Delta u) approx B^-1(u^(c) + Delta u - u^(b))  \n    D_h^T(u^(c))R^-1(w-h(u^(c))-D_h(u^(c))Delta u) approx B^-1(u^(c) + Delta u - u^(b)) \nendaligned","category":"page"},{"location":"dataAssimilation/3dvar/","page":"3d var","title":"3d var","text":"which we now solve for the update Delta u to obtain the linear system ","category":"page"},{"location":"dataAssimilation/3dvar/","page":"3d var","title":"3d var","text":"left(B^-1 + D_h^T(u^(c))R^-1D_h(u^(c)) right)Delta u = B^-1(u^(b)-u^(c)) + D_h^T(u^(c))R^-1(w-h(u^(c)))","category":"page"},{"location":"dataAssimilation/3dvar/","page":"3d var","title":"3d var","text":"Thus we have the following prescription ","category":"page"},{"location":"dataAssimilation/3dvar/","page":"3d var","title":"3d var","text":"To begin, take u^(c) == u^(b). \nSolve the system ","category":"page"},{"location":"dataAssimilation/3dvar/","page":"3d var","title":"3d var","text":"left(B^-1 + D_h^T(u^(c))R^-1D_h(u^(c)) right)Delta u = B^-1(u^(b)-u^(c)) + D_h^T(u^(c))R^-1(w-h(u^(c)))","category":"page"},{"location":"dataAssimilation/3dvar/","page":"3d var","title":"3d var","text":"to obtain Delta u","category":"page"},{"location":"dataAssimilation/3dvar/","page":"3d var","title":"3d var","text":"Update your guess using your favorite optimization algorithm. For example, in steppest descent, choose a learning rate eta and set ","category":"page"},{"location":"dataAssimilation/3dvar/","page":"3d var","title":"3d var","text":"u_textnew^(c)  = u_textprev^(c) + etaDelta u","category":"page"},{"location":"dataAssimilation/3dvar/","page":"3d var","title":"3d var","text":"Repeat the procedure until lvert u_textnew^(c) - u_textprev^(c) rvert converges to a desired tolerance.","category":"page"},{"location":"dataAssimilation/3dvar/","page":"3d var","title":"3d var","text":"In both the linear and nonlinear case, it should be noted that we have not added time indices to our state vectors. This is an indication that the 3d-var procedure is performed at every time where you have observation data.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = DataAssimilation","category":"page"},{"location":"#DataAssimilation","page":"Home","title":"DataAssimilation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for DataAssimilation.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Pages = [\n    \"index.md\",\n    \"dataAssimilation/\"\n    \"fdocs.md\",\n]\nDepth = 2","category":"page"}]
}
