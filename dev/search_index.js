var documenterSearchIndex = {"docs":
[{"location":"dataAssimilation/kalman/#Kalman-Filtering","page":"Kalman Filter","title":"Kalman Filtering","text":"","category":"section"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"Given some model for the error covariance matrices Q_k and R_k, we would like a method that propagates both our model and the errors forward. This way we may guarantee that the accuracy of our analysis doesn't come at the cost of higher uncertainty. ","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"The original implementation of the Kalman filter was for strictly linear systems. We will first develop the analysis for this simplified case adn then will generalize to the Extended Kalman Filter (EKF) that can handle fully nonlinear situations.","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"In the linear case, our system may be written as ","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"beginaligned\n    u_k+1^(t) = M_ku_k^(t) + xi_k+1^(p)  \n    w_k = H_ku_k^(t) + xi_k^(m)\nendaligned","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"where M_k and H_k are now matrices defining the linear problem. ","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"The goal of the Kalman filter is to derive the analysis u^(a) which optimizes the trace of the analysis error covariance matrix (i.e. sum of squared errors): ","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"mathrmTrleft( P_kright) = E(u_k^(t)-u_k^(a))^T(u_k^(t)-u_k^(a))","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"Finding the analysis consists of two steps: the forecast step and the assimilation step.","category":"page"},{"location":"dataAssimilation/kalman/#Forecast-Step","page":"Kalman Filter","title":"Forecast Step","text":"","category":"section"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"Assume we have the analysis at time t_k denoted u_k^(a). Then the forecast for time t_k+1 is","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"    u_k+1^(b) = M_ku_k^(a)","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"The background error is therefore ","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"beginaligned\n    xi_k+1^(b) = u_k+1^(t) - u_k+1^(b)  \n    = M_ku_k^(t)+xi_k+1^(p) - M_ku_k^(a)  \n    = M_kleft(u_k^(t)-u_k^(a) right) + xi_k+1^(p)  \n    = M_kxi_k^(a) + xi_k+1^(p)\nendaligned","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"We may now evaluate the covariance matrix of our background estimate as: ","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"beginaligned\n    B_k+1 = Exi_k+1^(b)(xi_k+1^(b))^T  \n    = Eleftleft(M_kxi_k^(a) + xi_k+1^p right) left(M_kxi_k^(a) + xi_k+1^p right)^T right  \nendaligned","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"If we presume that Exi_k^(b)(xi_k+1^(p))^T = 0, then the cross terms vanish and we are left with ","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"boxedB_k+1 = M_kP_kM_k^T + Q_k+1","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"Thus we now have the background (i.e forecast) estimate of the state at t_k+1 and its covariance matrix. Given a measurement w_k+1 at the same time with covariance matrix R_k+1, then we may now perform the assimilation step where we fuse the two sources of information to obtain u_k+1^(a) and P_k+1.","category":"page"},{"location":"dataAssimilation/kalman/#Data-Assimilation-Step","page":"Kalman Filter","title":"Data Assimilation Step","text":"","category":"section"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"Let's suppose that the analysis has the form ","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"u_k+1^(a) = nu + K_k+1w_k+1","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"for some vector nuinR^n and matrix K_k+1inR^mtimes n. In a perfect world, we would have Eu_k^(t)-u_k^(a) = 0. Therefore, ","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"beginaligned\n    0 = Eu_k^(t) - u_k^(a)  \n    = E(u_k^(b) + xi_k^(b)) - (nu + K_kw_k)  \n    = E(u_k^(b) + xi_k^(b)) - (nu + K_kH_ku_k^(t) + K_kxi_k^(m))  \n    = Eu_k^(b) + Exi_k^(b) - Enu -K_kH_kEu_k^(t) - K_kExi_k^(m)\n    = u_k^(b) + 0 - nu - K_kH_ku_k^(b) - 0  \n    = u_k^(b) - nu - K_kH_ku_k^(b)  \n    Rightarrow nu = u_k^(b) - K_kH_ku_k^(b)\nendaligned","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"which we now substitute to obtain ","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"boxedu_k^(a) = u_k^(b) + K_k(w_k - H_ku_k^(b))","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"Now that we know the form for the analysis we may derive the optimal matrix K_k by optimization of P_k. We have","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"beginaligned\n    xi_k^(a) = u_k^(t) - u_k^(a)  \n                = M_k-1u_k-1^(t) + xi_k^(p) - u_k^(b) - K_kleft(w_k - H_ku_k^(b) right)  \n                = M_k-1u_k-1^(t) + xi_k^(p) - M_k-1u_k-1^(a) - K_kleft(H_ku_k^(t) + xi_k^(m) - H_ku_k^(b) right)  \n                = M_k-1u_k-1^(t) + xi_k^(p) - M_k-1u_k-1^(a) - K_kH_ku_k^(t) - K_kxi_k^(m) + K_kH_ku_k^(b)  \n                = M_k-1u_k-1^(t) + xi_k^(p) - M_k-1u_k-1^(a) - K_kH_ku_k^(t) - K_kxi_k^(m) + K_kH_ku_k^(b)  \n                = Bigleft M_k-1(xi_k-1^(a)+u_k-1^(a)) + xi_k^(p) - M_k-1u_k-1^(a) Bigright - K_kH_ku_k^(t) - K_kxi_k^(m) + K_kH_ku_k^(b)  \n                = Bigleft M_k-1xi_k-1^(a) + xi_k^(p) Bigright - K_kH_ku_k^(t) + K_kH_ku_k^(b) - K_kxi_k^(m) \n                = M_k-1xi_k-1^(a) + xi_k^(p) - K_kH_k(M_k-1u_k-1^(t) + xi_k^(b)) + K_kH_ku_k^(b) - K_kxi_k^(m) \n                = M_k-1xi_k-1^(a) + xi_k^(p) - K_kH_kM_k-1(xi_k-1^(a) + u_k-1^a) - K_kH_kxi_k^(b) + K_kH_ku_k^(b) - K_kxi_k^(m) \n                = M_k-1xi_k-1^(a) + xi_k^(p) - K_kH_kM_k-1(xi_k-1^(a) + u_k-1^a) - K_kH_kxi_k^(b) + K_kH_kM_k-1u_k-1^(a) - K_kxi_k^(m) \n                = M_k-1xi_k-1^(a) + xi_k^(p) - K_kH_kM_k-1xi_k-1^(a) - K_kH_kxi_k^(b) - K_kxi_k^(m) \n                = big(I-K_kH_k big)(M_k-1xi_k-1^(a) - xi_k^p) - K_kxi_k^(m)\nendaligned","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"and therefore the covariance matrix is ","category":"page"},{"location":"dataAssimilation/kalman/","page":"Kalman Filter","title":"Kalman Filter","text":"beginaligned\n    P_k = Exi_k^(a)(xi_k^(a))^T  \n        = EBigleft(big(I-K_kH_k big)(M_k-1xi_k-1^(a) - xi_k^p) - K_kxi_k^(m) right) left(big(I-K_kH_k big)(M_k-1xi_k-1^(a) - xi_k^p) - K_kxi_k^(m) right)^T Big  \n        = big(I-K_kH_k big)M_k-1Exi_k-1^(a)(xi_k-1^(a))^TM_k-1^Tbig(I-K_kH_k big)^T + big(I-K_kH_k big)Exi_k^(p)(xi_k^(p))^Tbig(I-K_kH_k big)^T - K_kExi_k^(m)(xi_k^(m))^TK_k^T  \n        = big(I-K_kH_k big)B_kbig(I-K_kH_k big)^T - K_kR_kK_k^T\nendaligned","category":"page"},{"location":"dataAssimilation/ekf/#Extended-Kalman-Filter","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"","category":"section"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"For the Extended Kalman Filter we seek to find the x_k^a that minimizes xi_k^(a) = u_k^(t) - x_k^(t) in the least-squares sense. ","category":"page"},{"location":"dataAssimilation/ekf/#Initialization","page":"Extended Kalman Filter","title":"Initialization","text":"","category":"section"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"To begin, we apply the use-what-you-have strategy and set ","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"beginaligned\n    x_0^a = mu_0  \n    P_0 = E(x_0-x_0^a)(x_0-x_0^a)^T\nendaligned","category":"page"},{"location":"dataAssimilation/ekf/#Deriving-the-Forecast-Covariance","page":"Extended Kalman Filter","title":"Deriving the Forecast Covariance","text":"","category":"section"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"Consider now the forecast error produced by our model f at time index k. We have ","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"beginaligned\n    e_k^f = x_k - x_k^f  \n        = f(x_k-1) + w_k-1 - f(x_k-1^a)  \nendaligned","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"If our functions f and h are smooth enough (in this case, C^1), we may expand in a Taylor series and determine their linearization via the Jacobian about the analysis vector x_k-1^a. That is, ","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"f(x_k-1) approx f(x_k-1^a) + J_f(x_k-1^a)(x_k-1-x_k-1^a) + text Higher Order Terms","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"where ","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"J_f = left dfracpartial f_ipartial x_jright ","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"Therefore we can acheive the approximation ","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"beginaligned\n    e_k^f approx f(x_k-1^a) + J_f(x_k-1^a)e_k-1 + w_k-1 - f(x_k-1^a)  \n        = J_f(x_k-1^a)e_k-1 + w_k-1\nendaligned","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"We can now use this error vector to form the forecast error covariance matrix","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"beginaligned\n    P_k^f = Ee_k^f(e_k^f)^T  \n        approx E(J_f e_k-1 + w_k-1)(J_f e_k-1 + w_k-1)^T  \n        approx E(J_f e_k-1 + w_k-1)(e_k-1^T J_f^T  + w_k-1^T)  \n        approx J_fEe_k-1e_k-1^TJ_f^T + J_fEe_k-1w_k-1^T + Ew_k-1e_k-1^TJ_f^T + Ew_k-1w_k-1^T  \n        = J_f(x_k-1^a)P_k-1J_f^T(x_k-1^a) + Q_k-1\nendaligned","category":"page"},{"location":"dataAssimilation/ekf/#Deriving-Data-Assimilation-Step","page":"Extended Kalman Filter","title":"Deriving Data Assimilation Step","text":"","category":"section"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"At time index k we have x_k^f, P_k^f, z_k, and R_k. We will now use these to derive the optimal analysis x_k^a. ","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"Let's assume that the result is of the form","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"x_k^a = a + K_kz_k  in R^n","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"with ainR^n and K_kinR^mtimes n. We desire that Ex_k - x_k^a = 0. Thus, ","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"beginaligned\n    0 = Ex_k - x_k^a  \n    = E(x_k^f + e_k^f) - (a+K_k z_k)  \n    = E(x_k^f + e_k^f) - (a+K_k h(x_k)  +K_kv_k)  \n    = E_x_kx_k^f + E_x_ke_k^f - E_x_ka - K_kE_x_kh(x_k) - K_kE_x_kv_k  \n    = x_k^f + 0 - a - K_kE_x_kh(x_k)  - 0  \n    = x_k^f - a - K_kEh(x_k)  \n    Rightarrow a = x_k^f - K_kEh(x_k)\nendaligned","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"We now substitute to obtain","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"x_k^a - x_k^f - K_kEh(x_k) + K_kz_k","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"Since h is smooth enough, we may expand it about x_k^f to find that ","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"h(x_k) approx h(x_k^f) + J_h(x_k^f)(x_k-x_k^f)","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"where ","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"J_h = left dfracpartial h_ipartial x_j right","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"This leads to the approximation ","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"beginaligned\n    Eh(x_k) approx Eh(x_k^f) + J_h(x_k^f)e_k^f  \n        = Eh(x_k^f) + J_hEe_k^f  \n        = Eh(x_k^f) = h(x_k^f)\nendaligned","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"The analysis then becomes ","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"beginaligned\n    x_k^a = x_k^f - K_kEh(x_k) + K_k z_k  \n        approx x_k^f - K_kh(x_k^f) + K_kz_k  \n        = x_k^f + K_k(z_k - h(x_k^f))\nendaligned","category":"page"},{"location":"dataAssimilation/ekf/#Determining-the-Filter-Matrix-K_k","page":"Extended Kalman Filter","title":"Determining the Filter Matrix K_k","text":"","category":"section"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"To determine the K_k that will give us the optimal x_k^a, we will first derive the analysis error covariance matrix and then optimize its trace (i.e. the sum of squared errors) with respect to K_k. Following the same procedure as before, we begin by infestigating the error e_k=x_k-x_k^a.","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"beginaligned\n    e_k = x_k - x_k^a  \n    = f(x_k-1) + w_k-1 - x_k^f - K_k(z_k-h(x_k^f))  \n    approx f(x_k-1) + w_k-1 - f(x_k-1^a) - K_k(h(x_k) - h(x_k^f) + v_k)  \n    approx J_f(x_k-1^a)e_k-1 + w_k-1 - K_kJ_h(x_k^f)(J_f(x_k-1^a)e_k-1 + w_k-1) - K_kv_k  \n    = J_f(x_k-1^a) e_k-1 - K_k J_h(x_k^f)J_f(x_k-1^a)e_k-1 + w_k-1 - K_kJ_h(x_k^f) w_k-1 - K_kv_k  \n    = left(I-K_kJ_h(x_k^f) right)J_f(x_k-1^a)e_k-1  + left(I - K_kJ_h(x_k^f)right)w_k-1 - K_kv_k\nendaligned","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"We now form the covariance matrix","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"beginaligned\n    P_k = Ee_ke_k^T  \n    = left( I - K_k J_h(x_k^f)right)J_f(x_k-1^a)P_k-1J_f^T(x_k-1^a)left( I - K_k J_h(x_k^f)right)^T  \n     + left( I - K_k J_h(x_k^f)right)Q_k-1left( I - K_k J_h(x_k^f)right)^T + K_kR_kK_k^T  \n    = left( I - K_k J_h(x_k^f)right)left J_f(x_k-1^a)P_k-1 J_f^T(x_k-1^a) + Q_k-1  rightleft( I - K_k J_h(x_k^f)right)^T + K_kR_kK_k^T  \n    = left( I - K_k J_h(x_k^f)right)P_k^f left( I - K_k J_h(x_k^f)right)^T + K_kR_kK_k^T\nendaligned","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"We now seek to minimize texttr(P_k) with respect to K_k. The following identities will be helpful: ","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"beginaligned\n    mathopnabla_Atexttr(AB) = B^T  \n    mathopnabla_Atexttr(BA) = B  \n    mathopnabla_Atexttr(ABA) = AB^T + AB   \nendaligned","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"From these, we obtain ","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"beginaligned\n    0 = mathopnabla_K_ktexttr(P_k)  \n        = -left(J_h(x_k^f) P_k^f right)^T - P_k^fleft( J_h(x_k^f) right)^T + K_kleft( left J_h(x_k^f)P_k^f J_h^T(x_k^f)right^T  + J_h(x_k^f P_k^f J_h^T(x_k^f))right)  \n        + K_k(R_k^T + R_k)  \n        = -2P_k^f J_h^T(x_k^f) + 2K_kleft J_h(x_k^f) P_k^f J_h^T(x_k^f) + R_k right  \n    K_k = P_k^f J_h^T(x_k^f)left J_h(x_k^f)P_k^f J_h^T(x_k^f) + R_k right^-1\nendaligned","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"Lastly, we substitute this result back into P_k to derive a simpler expression. ","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"beginaligned\n    P_k = left(I - K_k J_h(x_k^f) right)P_k^fleft(I - K_k J_h(x_k^f)right)^T + K_kR_kK_k^T \n    = left(I-K_kJ_h(x_k^f)right)P_k^f - left( I - K_kJ_h(x_k^f) right)P_k^fleft( K_kJ_h(x_k^f)right)^T + K_kR_kK_k^T  \n    = left(I-K_kJ_h(x_k^f)right)P_k^f - left P_k^fJ_h^T(x_k^f) - K_k J_h(x_k^f)P_k^fJ_h^T(x_k^f) - K_kR_krightK_k^T  \n    = left(I-K_kJ_h(x_k^f)right)P_k^f - left P_k^fJ_h^T(x_k^f) - K_k left( J_h(x_k^f)P_k^fJ_h^T(x_k^f) - R_kright) rightK_k^T  \n    = left(I-K_kJ_h(x_k^f)right)P_k^f - leftP_k^fJ_h^T(x_k^f) - P_k^fJ_h^T(x_k^f)rightK_k^T  \n    P_k = left(I-K_kJ_h(x_k^f)right)P_k^f   \nendaligned","category":"page"},{"location":"dataAssimilation/ekf/#Summary","page":"Extended Kalman Filter","title":"Summary","text":"","category":"section"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"Let's summarize the whole process:  We have a dynamical system of the form ","category":"page"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"beginaligned\n    x_k = f(x_-1) + w_k-1  \n    z_k = h(x_k) + v_k\nendaligned","category":"page"},{"location":"dataAssimilation/ekf/#.-Initialization","page":"Extended Kalman Filter","title":"0. Initialization","text":"","category":"section"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"beginaligned\n    x_0^a = mu_0  \n    P_0 = E(x_0-x_0^a)(x_0-x_0^a)^T \nendaligned","category":"page"},{"location":"dataAssimilation/ekf/#.-Forecast-Step","page":"Extended Kalman Filter","title":"1. Forecast Step","text":"","category":"section"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"beginaligned\n    x_k^f = f(x_k-1^a)  \n    P_k^f = J_f(x_k-1^a)P_k-1J_f^T(x_k-1) + Q_k-1\nendaligned","category":"page"},{"location":"dataAssimilation/ekf/#.-Assimilation","page":"Extended Kalman Filter","title":"2. Assimilation","text":"","category":"section"},{"location":"dataAssimilation/ekf/","page":"Extended Kalman Filter","title":"Extended Kalman Filter","text":"beginaligned\n    x_k^a = x_k^f + K_k(z_k - h(x_k^f))  \n    K_k = P_k^fJ_h^T(x_k^f)leftJ_h(x_k^f)P_k^f(J_h^T(x_k^f) + R_k right^-1  \n    P_k = left(I - K_kJ_h(x_k^f) right)P_k^f\nendaligned","category":"page"},{"location":"fdocs/#Function-Documentation","page":"Function docs","title":"Function Documentation","text":"","category":"section"},{"location":"fdocs/","page":"Function docs","title":"Function docs","text":"Modules = [DataAssimilation]","category":"page"},{"location":"dataAssimilation/dataassim/#Data-Assimilation-Overview","page":"Data Assimilation Oveview","title":"Data Assimilation Overview","text":"","category":"section"},{"location":"dataAssimilation/dataassim/#Overview","page":"Data Assimilation Oveview","title":"Overview","text":"","category":"section"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"The proper application of scientific models to make real-world predictions requires that we commit ourselves to a full accounting of all possible sources of uncertainty when reporting results. Further, the explosion of big data across scientific fields provieds a plethora observational data that our models are typically unequipped to incorporate when making predictions. The field of Data Assimilation addresses this problem by providing a family of techniques engineered to combine model output together with observational data whilst enabling a complete accounting the sources of uncertainty. For chaotic systems in particular, data assimilation enables integration on long time scales that would be impossible via models alone. ","category":"page"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"In this overview, we will follow the examples from this nice paper. ","category":"page"},{"location":"dataAssimilation/dataassim/#Framing-the-Problem","page":"Data Assimilation Oveview","title":"Framing the Problem","text":"","category":"section"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"Data assimilation can be understood most generally in terms of dyscrete dynamical systems. This enables us to apply the methods to most mathematical models from gridded PDE solvers to systems of ordinary differential equations. Our goal is to find the best prediction for the system state vector u that combines our model predictions, also known as forecasts, with observational data. Model predictions are summarized via the discrete update equation: ","category":"page"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"u_k+1 = mathcalM(u_k theta)","category":"page"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"For ODE systems, mathcalM represents the time integration scheme for a system of ODEs like ","category":"page"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"dfracdudt = f(u t theta)","category":"page"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"To measure the performance of our assimilation scheme, we denote the true value of the state vector asd u^(t). The output of our model is denoted u^(b) (b subscript for background). The discrepancy between the true value and our forecast is denoted xi^(b) = u^(t) - u^(b) characterizing the extent to which our model prediction is imperfect. ","category":"page"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"The observations of our system are denoted by w_k = w(t_k). These observations do not neccessarily need to be components of the state vector u, but rather, are related to it via the observation function h. For example, one may attempt to predict sea surface temperature using data assimilation with data from satellite observations. The function h would then be the Stefan-Boltzmann law. However, real world data is noisy, which we must take into account. We write ","category":"page"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"w_k = h(u_k) + xi_k^(m)","category":"page"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"where xi_k^(m) denotes this measurement noise. ","category":"page"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"Given our model predictions u_k^(b) and observations w_k, we seek to obtain the optimal or best-possible prediction called the analysis, u^(a). This analysis will still not be perfect, so we further specify the analysis error via ","category":"page"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"xi^(a) = u^(t) - u^(a)","category":"page"},{"location":"dataAssimilation/dataassim/#Summary","page":"Data Assimilation Oveview","title":"Summary","text":"","category":"section"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"beginaligned\n    u_k^(t) in R^n textthe true state vector  \n    u_k^(b) in R^n textthe kth model forecast  \n    u_k^(a) in R^n textthe analysis  \n    w_k in R^m textthe kth observation vector  \n    xi^(b) in R^n textthe model forecast error\n    xi^(m) in R^m textthe observation noise vector \n    xi^(a) in R^n textthe analysis error\n    xi^(p) in R^n textthe process noise if we used our model on the true state\n    mathcalMR^ntoR^n textthe model update function\n    fR^ntoR^n textdifferential equation model \n    hR^ntoR^m  textobservation function\nendaligned","category":"page"},{"location":"dataAssimilation/dataassim/#Assumptions","page":"Data Assimilation Oveview","title":"Assumptions","text":"","category":"section"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"To make possible the derivation of a unique analysis u^(a), the following assumptions are in order.","category":"page"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"beginaligned\n    Exi_k^(b) = 0  Exi_k^(b)(xi_j^(b))^T = 0 text for  kneq j\n    Exi_k^(m) = 0  Exi_k^(m)(xi_j^(m))^T = 0 text for  kneq j\n    Exi_k^(b)(u_0)^T = 0  Exi_k^(m)(u_0)^T = 0\n    Exi_k^(b)xi_j^(m) = 0     \n    Eu_k^(t) = u_k^(b)  \nendaligned","category":"page"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"We also define the error covariance matrices","category":"page"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"beginaligned\n    Q_k = Exi_k^(p)(xi_k^(p))^T  \n    R_k = Exi_k^(m)(xi_k^(m))^T  \n    B_k = Exi_k^(b)(xi_k^(b))^T \nendaligned","category":"page"},{"location":"dataAssimilation/dataassim/","page":"Data Assimilation Oveview","title":"Data Assimilation Oveview","text":"which we will use in our consideration of the final error of our analysis.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = DataAssimilation","category":"page"},{"location":"#DataAssimilation","page":"Home","title":"DataAssimilation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for DataAssimilation.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Pages = [\n    \"index.md\",\n    \"dataAssimilation/\"\n    \"fdocs.md\",\n]\nDepth = 2","category":"page"}]
}
