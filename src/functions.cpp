#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>

using namespace std;
using namespace Rcpp;


// [[Rcpp::export]]
List make_initial(arma::vec Y, arma::vec S, arma::mat Z, arma::mat X){
  int n = X.n_rows;
  
  arma::mat Px = X * (X.t() * X).i() * X.t();
  arma::mat W = arma::join_rows(S, Z);
  arma::vec betagamma = (W.t() * Px * W).i() * W.t() * Px * Y;
  
  arma::vec alpha_init = (X.t() * X).i() * X.t() * S;
  arma::vec beta_init = betagamma(arma::span(1, betagamma.n_elem - 1));
  double gamma_init = betagamma(0);
  
  arma::vec res1 = S - X * alpha_init;
  double sigmaSsquare_init = arma::as_scalar(res1.t() * res1) / n;
  
  arma::vec res2 = Y - gamma_init * S - Z * beta_init;
  double sigmaYsquare_init = arma::as_scalar(res2.t() * res2) / n;
  
  double rho_init = arma::as_scalar(res1.t() * res2) / (n * sqrt(sigmaSsquare_init) * sqrt(sigmaYsquare_init));
  
  return List::create(alpha_init, beta_init, gamma_init, sigmaSsquare_init, sigmaYsquare_init, rho_init);
}






// No Missing Data on Exposure
// MLE of alpha, beta, gamma
arma::vec alpha_beta_gamma_compute(arma::mat X, arma::vec Y, arma::mat Z, arma::vec S,
                                   arma::vec alpha, arma::vec beta, double gamma){
  int n = Y.n_elem;
  int npara = alpha.n_elem + beta.n_elem + 1;
  
  arma::vec parameter_old = arma::zeros(npara);
  arma::vec parameter_new = arma::zeros(npara);
  arma::vec parameter_change = arma::zeros(npara);
  arma::vec parameter_new_temp = arma::zeros(npara);
  
  arma::vec alpha_new = arma::zeros(alpha.n_elem);
  arma::vec alpha_new_temp = arma::zeros(alpha.n_elem);
  arma::vec beta_new = arma::zeros(beta.n_elem);
  arma::vec beta_new_temp = arma::zeros(beta.n_elem);
  double gamma_new = 0;
  double gamma_new_temp = 0;
  
  arma::mat HessianMat = arma::zeros(npara, npara);
  arma::vec grad = arma::zeros(npara);
  
  parameter_new(arma::span(0, alpha.n_elem - 1)) = alpha;
  parameter_new(arma::span(alpha.n_elem, alpha.n_elem + beta.n_elem - 1)) = beta;
  parameter_new(alpha.n_elem + beta.n_elem) = gamma;
  
  arma::vec eS = arma::zeros(n);
  arma::vec eY = arma::zeros(n);
  arma::vec eS_temp = arma::zeros(n);
  arma::vec eY_temp = arma::zeros(n);
  arma::vec XteS = arma::zeros(alpha.n_elem);
  arma::vec XteY = arma::zeros(alpha.n_elem);
  arma::vec ZteS = arma::zeros(beta.n_elem);
  arma::vec ZteY = arma::zeros(beta.n_elem);
  double A1, A2, A3, SteS, SteY, A1_temp, A2_temp, A3_temp;
  int m;
  
  double dist = 1;
  double epsilon = pow(10, -5);
  
  while(dist > epsilon){
    parameter_old = parameter_new;
    
    alpha_new = parameter_new(arma::span(0, alpha.n_elem-1));
    beta_new = parameter_new(arma::span(alpha.n_elem, alpha.n_elem + beta.n_elem - 1));
    gamma_new = parameter_new(alpha.n_elem + beta.n_elem);
    
    eS = S - X * alpha_new;
    eY = Y - gamma_new * S - Z * beta_new;
    XteS = X.t() * eS;
    XteY = X.t() * eY;
    ZteS = Z.t() * eS;
    ZteY = Z.t() * eY;
    SteS = arma::as_scalar(S.t() * eS);
    SteY = arma::as_scalar(S.t() * eY);
    A1 = arma::as_scalar(eS.t() * eS);
    A2 = arma::as_scalar(eY.t() * eY);
    A3 = arma::as_scalar(eY.t() * eS);
    
    HessianMat(arma::span(0, alpha.n_elem-1), arma::span(0, alpha.n_elem-1))
      = 2*A2*X.t()*X - 2*XteY*XteY.t()
      - (2*A3*XteY - 2*A2*XteS) * (2*A3*XteY - 2*A2*XteS).t() / (A1*A2-A3*A3);
    HessianMat(arma::span(0, alpha.n_elem-1), arma::span(alpha.n_elem, alpha.n_elem + beta.n_elem - 1))
      = 4*XteS*ZteY.t() - 2*XteY*ZteS.t() - 2*A3*X.t()*Z
      - (2*A3*XteY - 2*A2*XteS) * (2*A3*ZteS - 2*A1*ZteY).t() / (A1*A2-A3*A3);
    HessianMat(arma::span(0, alpha.n_elem-1), alpha.n_elem + beta.n_elem)
      = 4*XteS*SteY - 2*XteY*SteS - 2*A3*X.t()*S
      - (2*A3*SteS - 2*A1*SteY) * (2*A3*XteY - 2*A2*XteS) / (A1*A2-A3*A3);
          
    HessianMat(arma::span(alpha.n_elem, alpha.n_elem + beta.n_elem - 1), arma::span(0, alpha.n_elem-1))
      = 4*ZteY*XteS.t() - 2*ZteS*XteY.t() - 2*A3*Z.t()*X
      - (2*A3*ZteS - 2*A1*ZteY) * (2*A3*XteY - 2*A2*XteS).t() / (A1*A2-A3*A3);
    HessianMat(arma::span(alpha.n_elem, alpha.n_elem + beta.n_elem - 1), arma::span(alpha.n_elem, alpha.n_elem + beta.n_elem - 1))
      = 2*A1*Z.t()*Z - 2*ZteS*ZteS.t()
      - (2*A3*ZteS - 2*A1*ZteY) * (2*A3*ZteS - 2*A1*ZteY).t() / (A1*A2-A3*A3);
    HessianMat(arma::span(alpha.n_elem, alpha.n_elem + beta.n_elem - 1), alpha.n_elem + beta.n_elem)
      = 2*A1*Z.t()*S - 2*ZteS*SteS
      - (2*A3*SteS - 2*A1*SteY) * (2*A3*ZteS - 2*A1*ZteY) / (A1*A2-A3*A3);
                
    HessianMat(alpha.n_elem + beta.n_elem, arma::span(0, alpha.n_elem-1))
      = 4*SteY*XteS.t() - 2*SteS*XteY.t() - 2*A3*S.t()*X
      - (2*A3*SteS - 2*A1*SteY) * (2*A3*XteY - 2*A2*XteS).t() / (A1*A2-A3*A3);
    HessianMat(alpha.n_elem + beta.n_elem, arma::span(alpha.n_elem, alpha.n_elem + beta.n_elem - 1))
      = 2*A1*S.t()*Z - 2*SteS*ZteS.t()
      - (2*A3*SteS - 2*A1*SteY) * (2*A3*ZteS - 2*A1*ZteY).t() / (A1*A2-A3*A3);
    HessianMat(alpha.n_elem + beta.n_elem, alpha.n_elem + beta.n_elem)
      = 2*A1*arma::as_scalar(S.t()*S) - 2*pow(SteS, 2)
      - (2*A3*SteS - 2*A1*SteY) * (2*A3*SteS - 2*A1*SteY) / (A1*A2-A3*A3);
                      
    grad(arma::span(0, alpha.n_elem-1)) = 2*A3*XteY - 2*A2*XteS;
    grad(arma::span(alpha.n_elem, alpha.n_elem + beta.n_elem - 1)) = 2*A3*ZteS - 2*A1*ZteY;
    grad(alpha.n_elem + beta.n_elem) = 2*A3*SteS - 2*A1*SteY;
                      
    parameter_change = HessianMat.i() * grad;
                      
    parameter_new_temp = parameter_old - parameter_change;
    alpha_new_temp = parameter_new_temp(arma::span(0, alpha.n_elem-1));
    beta_new_temp = parameter_new_temp(arma::span(alpha.n_elem, alpha.n_elem + beta.n_elem - 1));
    gamma_new_temp = parameter_new_temp(alpha.n_elem + beta.n_elem);
                      
    eS_temp = S - X * alpha_new_temp;
    eY_temp = Y - gamma_new_temp * S - Z * beta_new_temp;
    A1_temp = arma::as_scalar(eS_temp.t() * eS_temp);
    A2_temp = arma::as_scalar(eY_temp.t() * eY_temp);
    A3_temp = arma::as_scalar(eY_temp.t() * eS_temp);
                      
    m = 0;
    while( (A1*A2-A3*A3) < (A1_temp*A2_temp-A3_temp*A3_temp) ){
      m = m + 1;
      
      parameter_new_temp = parameter_old - parameter_change / pow(2, m);
      alpha_new_temp = parameter_new_temp(arma::span(0, alpha.n_elem-1));
      beta_new_temp = parameter_new_temp(arma::span(alpha.n_elem, alpha.n_elem + beta.n_elem - 1));
      gamma_new_temp = parameter_new_temp(alpha.n_elem + beta.n_elem);
      
      eS_temp = S - X * alpha_new_temp;
      eY_temp = Y - gamma_new_temp * S - Z * beta_new_temp;
      A1_temp = arma::as_scalar(eS_temp.t() * eS_temp);
      A2_temp = arma::as_scalar(eY_temp.t() * eY_temp);
      A3_temp = arma::as_scalar(eY_temp.t() * eS_temp);
    }
    
    parameter_new = parameter_new_temp;
    
    dist = pow(arma::as_scalar(parameter_change.t() * parameter_change), 0.5) / pow(2, m);
    //dist = pow(arma::as_scalar(grad.t() * grad), 0.5);
    //cout << parameter_new << endl;
                      
    //for(int i=0; i<temp2_new.n_elem; i++)
    //  cout << temp2_new(i) << " ";
    //cout << endl;
  }
  //cout << "A" << " ";
  return parameter_new;
}



// MLE of sigmaSsquare
double sigmaSsquare_compute(arma::vec alpha, arma::mat X, arma::vec S){
  arma::mat mat_temp1 = S - X * alpha;
  //cout << "B" << " ";
  return arma::as_scalar(mat_temp1.t() * mat_temp1) / X.n_rows;
}



// MLE of sigmaYsquare
double sigmaYsquare_compute(arma::vec Y, arma::mat Z, arma::vec S, arma::vec beta, double gamma){
  arma::mat mat_temp1 = Y - gamma * S - Z * beta;
  //cout << "C" << " ";
  return arma::as_scalar(mat_temp1.t() * mat_temp1) / Z.n_rows;
}



// MLE of rho
double rho_compute(arma::mat X, arma::vec Y, arma::mat Z, arma::vec S,
                   arma::vec alpha, arma::vec beta, double gamma, double sigmaSsquare, double sigmaYsquare){
  arma::mat mat_temp1 = Y - gamma * S - Z * beta;
  arma::mat mat_temp2 = S - X * alpha;
  //cout << "D" << " ";
  return arma::as_scalar(mat_temp1.t() * mat_temp2) / ( Z.n_rows * pow(sigmaSsquare, 0.5) * pow(sigmaYsquare, 0.5) );
}



// Estimated Covariance Matrix
arma::mat est_cov_compute(arma::mat X, arma::vec Y, arma::mat Z, arma::vec S,
                          arma::vec alpha, arma::vec beta, double gamma,
                          double sigmaSsquare, double sigmaYsquare, double rho){
  int n = X.n_rows;
  arma::mat Ic = arma::zeros(X.n_cols + Z.n_cols + 4, X.n_cols + Z.n_cols + 4);
  
  arma::vec eS = S - X * alpha;
  arma::vec eY = Y - gamma * S - Z * beta;
  double sigmaS = pow(sigmaSsquare, 0.5);
  double sigmaY = pow(sigmaYsquare, 0.5);
  arma::vec XteS = X.t() * eS;
  arma::vec XteY = X.t() * eY;
  arma::vec ZteS = Z.t() * eS;
  arma::vec ZteY = Z.t() * eY;
  double SteS = arma::as_scalar(S.t() * eS);
  double SteY = arma::as_scalar(S.t() * eY);
  //double eSteS = arma::as_scalar(eS.t() * eS);
  //double eSteY = arma::as_scalar(eS.t() * eY);
  //double eYteY = arma::as_scalar(eY.t() * eY);
  
  Ic(arma::span(0, X.n_cols-1), arma::span(0, X.n_cols-1)) = X.t()*X / ( (1-pow(rho,2)) * sigmaSsquare );
  Ic(arma::span(0, X.n_cols-1), arma::span(X.n_cols, X.n_cols+Z.n_cols-1)) = -rho*X.t()*Z / ( (1-pow(rho,2))*sigmaS*sigmaY );
  Ic(arma::span(0, X.n_cols-1), X.n_cols+Z.n_cols) = -rho*X.t()*S / ( (1-pow(rho,2))*sigmaS*sigmaY );
  Ic(arma::span(0, X.n_cols-1), X.n_cols+Z.n_cols+1) = XteS / ( (1-rho*rho) * pow(sigmaSsquare, 2) ) - rho*XteY / (2 * (1-rho*rho) * pow(sigmaS, 3) * sigmaY);
  Ic(arma::span(0, X.n_cols-1), X.n_cols+Z.n_cols+2) = -rho*XteY / (2 * (1-rho*rho) * sigmaS * pow(sigmaY, 3));
  Ic(arma::span(0, X.n_cols-1), X.n_cols+Z.n_cols+3) = -2*rho*XteS / ((1-rho*rho) * (1-rho*rho) * sigmaSsquare) + (1+rho*rho)*XteY / ((1-rho*rho) * (1-rho*rho) * sigmaS * sigmaY);
  
  Ic(arma::span(X.n_cols, X.n_cols + Z.n_cols - 1), arma::span(0, X.n_cols-1)) = -rho*Z.t()*X / ( (1-rho*rho)*sigmaS*sigmaY );
  Ic(arma::span(X.n_cols, X.n_cols + Z.n_cols - 1), arma::span(X.n_cols, X.n_cols+Z.n_cols-1)) = Z.t()*Z / ( (1-rho*rho)*sigmaYsquare );
  Ic(arma::span(X.n_cols, X.n_cols + Z.n_cols - 1), X.n_cols+Z.n_cols) = Z.t()*S / ( (1-rho*rho)*sigmaYsquare );
  Ic(arma::span(X.n_cols, X.n_cols + Z.n_cols - 1), X.n_cols+Z.n_cols+1) = -rho*ZteS / (2 * (1-rho*rho) * pow(sigmaS, 3) * sigmaY);
  Ic(arma::span(X.n_cols, X.n_cols + Z.n_cols - 1), X.n_cols+Z.n_cols+2) = ZteY / ( (1-rho*rho) * pow(sigmaYsquare, 2) ) - rho*ZteS / (2 * (1-rho*rho) * sigmaS * pow(sigmaY, 3));
  Ic(arma::span(X.n_cols, X.n_cols + Z.n_cols - 1), X.n_cols+Z.n_cols+3) = -2*rho*ZteY / ((1-rho*rho) * (1-rho*rho) * sigmaYsquare) + (1+rho*rho)*ZteS / ((1-rho*rho) * (1-rho*rho) * sigmaS * sigmaY);
  
  Ic(X.n_cols + Z.n_cols, arma::span(0, X.n_cols-1)) = -rho*S.t()*X / ( (1-rho*rho)*sigmaS*sigmaY );
  Ic(X.n_cols + Z.n_cols, arma::span(X.n_cols, X.n_cols+Z.n_cols-1)) = S.t()*Z / ( (1-rho*rho)*sigmaYsquare );
  Ic(X.n_cols + Z.n_cols, X.n_cols+Z.n_cols) = arma::as_scalar(S.t()*S) / ( (1-rho*rho)*sigmaYsquare );
  Ic(X.n_cols + Z.n_cols, X.n_cols+Z.n_cols+1) = -rho*SteS / (2 * (1-rho*rho) * pow(sigmaS, 3) * sigmaY);
  Ic(X.n_cols + Z.n_cols, X.n_cols+Z.n_cols+2) = SteY / ( (1-rho*rho) * pow(sigmaYsquare, 2) ) - rho*SteS / (2 * (1-rho*rho) * sigmaS * pow(sigmaY, 3));
  Ic(X.n_cols + Z.n_cols, X.n_cols+Z.n_cols+3) = -2*rho*SteY / ((1-rho*rho) * (1-rho*rho) * sigmaYsquare) + (1+rho*rho)*SteS / ((1-rho*rho) * (1-rho*rho) * sigmaS * sigmaY);
  
  Ic(X.n_cols + Z.n_cols+1, arma::span(0, X.n_cols-1)) = XteS.t() / ( (1-rho*rho) * pow(sigmaSsquare, 2) ) - rho*XteY.t() / (2 * (1-rho*rho) * pow(sigmaS, 3) * sigmaY);
  Ic(X.n_cols + Z.n_cols+1, arma::span(X.n_cols, X.n_cols+Z.n_cols-1)) = -rho*ZteS.t() / (2 * (1-rho*rho) * pow(sigmaS, 3) * sigmaY);
  Ic(X.n_cols + Z.n_cols+1, X.n_cols+Z.n_cols) = Ic(X.n_cols + Z.n_cols, X.n_cols+Z.n_cols+1);
  Ic(X.n_cols + Z.n_cols+1, X.n_cols+Z.n_cols+1) = n*(2.0-rho*rho) / ( 4*(1-rho*rho)*pow(sigmaSsquare, 2) );
  Ic(X.n_cols + Z.n_cols+1, X.n_cols+Z.n_cols+2) = -n*rho*rho / ( 4*(1-rho*rho)*sigmaSsquare*sigmaYsquare );
  Ic(X.n_cols + Z.n_cols+1, X.n_cols+Z.n_cols+3) = -n*rho / ( 2*(1-rho*rho)*sigmaSsquare );
  
  Ic(X.n_cols + Z.n_cols+2, arma::span(0, X.n_cols-1)) = -rho*XteY.t() / (2 * (1-rho*rho) * sigmaS * pow(sigmaY, 3));
  Ic(X.n_cols + Z.n_cols+2, arma::span(X.n_cols, X.n_cols+Z.n_cols-1)) = ZteY.t() / ( (1-rho*rho) * pow(sigmaYsquare, 2) ) - rho*ZteS.t() / (2 * (1-rho*rho) * sigmaS * pow(sigmaY, 3));
  Ic(X.n_cols + Z.n_cols+2, X.n_cols+Z.n_cols) = Ic(X.n_cols + Z.n_cols, X.n_cols+Z.n_cols+2);
  Ic(X.n_cols + Z.n_cols+2, X.n_cols+Z.n_cols+1) = Ic(X.n_cols + Z.n_cols+1, X.n_cols+Z.n_cols+2);
  Ic(X.n_cols + Z.n_cols+2, X.n_cols+Z.n_cols+2) = n*(2.0-rho*rho) / ( 4*(1-rho*rho)*pow(sigmaYsquare, 2) );
  Ic(X.n_cols + Z.n_cols+2, X.n_cols+Z.n_cols+3) = -n*rho / ( 2*(1-rho*rho)*sigmaYsquare );
  
  Ic(X.n_cols + Z.n_cols+3, arma::span(0, X.n_cols-1)) = -2*rho*XteS.t() / ((1-rho*rho) * (1-rho*rho) * sigmaSsquare) + (1+rho*rho)*XteY.t() / ((1-rho*rho) * (1-rho*rho) * sigmaS * sigmaY);
  Ic(X.n_cols + Z.n_cols+3, arma::span(X.n_cols, X.n_cols+Z.n_cols-1)) = -2*rho*ZteY.t() / ((1-rho*rho) * (1-rho*rho) * sigmaYsquare) + (1+rho*rho)*ZteS.t() / ((1-rho*rho) * (1-rho*rho) * sigmaS * sigmaY);
  Ic(X.n_cols + Z.n_cols+3, X.n_cols+Z.n_cols) = Ic(X.n_cols + Z.n_cols, X.n_cols+Z.n_cols+3);
  Ic(X.n_cols + Z.n_cols+3, X.n_cols+Z.n_cols+1) = Ic(X.n_cols + Z.n_cols+1, X.n_cols+Z.n_cols+3);
  Ic(X.n_cols + Z.n_cols+3, X.n_cols+Z.n_cols+2) = Ic(X.n_cols + Z.n_cols+2, X.n_cols+Z.n_cols+3);
  Ic(X.n_cols + Z.n_cols+3, X.n_cols+Z.n_cols+3) = n*(1.0+rho*rho) / ( (1-rho*rho)*(1-rho*rho) );
  
  //cout << "E" << "  ";
  return Ic.i();
}



struct para_est_nomissing{
  arma::vec alpha;
  arma::vec beta;
  double gamma;
  double sigmaSsquare;
  double sigmaYsquare;
  double rho;
  arma::mat cov;
};



struct para_est_nomissing MR_MLE_NoMissing(arma::mat X, arma::vec Y, arma::mat Z, arma::vec S,
                                           arma::vec alpha_init, arma::vec beta_init, double gamma_init,
                                           double sigmaSsquare_init, double sigmaYsquare_init, double rho_init, 
                                           double epsilon){
  struct para_est_nomissing mle;
  int npara = alpha_init.n_elem + beta_init.n_elem + 4;
  
  //cout << "A" << "  ";
  arma::vec alpha_new = alpha_init;
  arma::vec beta_new = beta_init;
  double gamma_new = gamma_init;
  double sigmaSsquare_new = sigmaSsquare_init;
  double sigmaYsquare_new = sigmaYsquare_init;
  double rho_new = rho_init;
  arma::mat est_covariance_mat = arma::zeros(npara, npara);
  arma::vec alpha_beta_gamma = arma::zeros(npara - 3);
  
  alpha_beta_gamma = alpha_beta_gamma_compute(X, Y, Z, S, alpha_init, beta_init, gamma_init);
  alpha_new = alpha_beta_gamma.subvec(0, X.n_cols - 1);
  beta_new = alpha_beta_gamma.subvec(X.n_cols, X.n_cols + Z.n_cols - 1);
  gamma_new = alpha_beta_gamma(X.n_cols + Z.n_cols);
  //for(int i=0; i<alpha_beta_gamma.n_elem; i++)
  //  cout << alpha_beta_gamma(i) << "  ";
  //cout << endl;
  //cout << "1" << " ";
  
  sigmaSsquare_new = sigmaSsquare_compute(alpha_new, X, S);
  //cout << "2" << "  ";
  
  sigmaYsquare_new = sigmaYsquare_compute(Y, Z, S, beta_new, gamma_new);
  //cout << "3" << "  ";
  
  rho_new = rho_compute(X, Y, Z, S, alpha_new, beta_new, gamma_new, sigmaSsquare_new, sigmaYsquare_new);
  //cout << "4" << endl;
  
  est_covariance_mat = est_cov_compute(X, Y, Z, S, alpha_new, beta_new, gamma_new, sigmaSsquare_new, sigmaYsquare_new, rho_new);
  
  mle.alpha = alpha_new;
  mle.beta = beta_new;
  mle.gamma = gamma_new;
  mle.sigmaSsquare = sigmaSsquare_new;
  mle.sigmaYsquare = sigmaYsquare_new;
  mle.rho = rho_new;
  mle.cov = est_covariance_mat;
  
  return mle;
}




// [[Rcpp::export]]
List MR_estimation_NoMissing_CORR(arma::mat X, arma::vec Y, arma::mat Z, arma::vec S,
                                  arma::vec alpha_init, arma::vec beta_init, double gamma_init,
                                  double sigmaSsquare_init, double sigmaYsquare_init, double rho_init,
                                  double epsilon){
  struct para_est_nomissing mle = MR_MLE_NoMissing(X, Y, Z, S, alpha_init, beta_init, gamma_init, sigmaSsquare_init, sigmaYsquare_init, rho_init, epsilon);
  return List::create(mle.alpha, mle.beta, mle.gamma, mle.sigmaSsquare, mle.sigmaYsquare, mle.rho, mle.cov);
}






// Missing Data on Exposure
// a update
double a_update(double gamma, double sigmaSsquare, double sigmaYsquare, double rho){
  return (1.0 - pow(rho, 2)) / (gamma*gamma/sigmaYsquare + 1.0/sigmaSsquare + 2.0*rho*gamma/sqrt(sigmaSsquare*sigmaYsquare));
}



// b update
arma::vec b_update(double a, arma::vec alpha, arma::vec beta, double gamma, double sigmaSsquare, double sigmaYsquare, double rho,
                   arma::mat X, arma::mat Z, arma::vec Y, arma::vec R){
  int n = Y.n_elem;
  arma::vec b=arma::zeros(n);
  
  for(int i=0; i<n; i++){
    if(R(i) != 0){
      b(i) = a / (1 - pow(rho, 2)) * (gamma * ( Y(i) - arma::as_scalar(Z.row(i) * beta) ) / sigmaYsquare +
        arma::as_scalar(X.row(i) * alpha) / sigmaSsquare +
        rho * ( Y(i) - arma::as_scalar(Z.row(i) * beta) + gamma * arma::as_scalar(X.row(i) * alpha) ) / sqrt(sigmaSsquare*sigmaYsquare) );
    }
  }
  
  return b;
}



// E(Si) and E(Si^2) update
arma::mat ES12_update(double a, arma::vec b, arma::vec R, arma::vec S){
  int n = S.n_elem;
  double L0i, U0i;
  arma::mat ES12=arma::zeros(n, 2);
  
  for(int i=0; i<n; i++){
    if (R(i) == 0) {
      ES12(i, 0) = S(i);
      ES12(i, 1) = pow(S(i), 2) ;
    }
    else if (R(i) == 3) {
      ES12(i, 0) = b(i);
      ES12(i, 1) = pow(b(i), 2) + a ;
    }
    else if (R(i) == 1) {
      L0i = (S(i) - b(i)) / pow(a, 0.5);
      ES12(i, 0) = b(i) - pow(a, 0.5) * arma::normpdf(L0i) / arma::normcdf(L0i);
      ES12(i, 1) = pow(b(i), 2) + a - arma::normpdf(L0i) / arma::normcdf(L0i) * (2*b(i)*pow(a, 0.5) + a*L0i);
    }
    else{
      U0i = (-S(i) + b(i)) / pow(a, 0.5);
      ES12(i, 0) = b(i) + pow(a, 0.5) * arma::normpdf(U0i) / arma::normcdf(U0i);
      ES12(i, 1) = pow(b(i), 2) + a + arma::normpdf(U0i) / arma::normcdf(U0i) * (2*b(i)*pow(a, 0.5) - a*U0i);
    }
  }
  
  return ES12;
}



// alpha, beta, gamma update
arma::vec alpha_beta_gamma_update(arma::mat X, arma::vec Y, arma::mat Z, arma::mat ES12,
                                  arma::vec alpha, arma::vec beta, double gamma){
  int n = Y.n_elem;
  int npara = alpha.n_elem + beta.n_elem + 1;
  
  arma::vec parameter_old = arma::zeros(npara);
  arma::vec parameter_new = arma::zeros(npara);
  arma::vec parameter_new_temp = arma::zeros(npara);
  arma::vec parameter_change = arma::zeros(npara);
  
  arma::mat HessianMat = arma::zeros(npara, npara);
  arma::vec grad = arma::zeros(npara);
  
  arma::vec alpha_new_temp = arma::zeros(alpha.n_elem);
  arma::vec beta_new_temp = arma::zeros(beta.n_elem);
  double gamma_new_temp = 0;
  
  arma::vec eS = arma::zeros(n);
  arma::vec eY = arma::zeros(n);
  arma::vec eS_temp = arma::zeros(n);
  arma::vec eY_temp = arma::zeros(n);
  double SS, YY, SY, SteS, SteY, SS_temp, YY_temp, SY_temp;
  arma::vec XteS = arma::zeros(alpha.n_elem);
  arma::vec XteY = arma::zeros(alpha.n_elem);
  arma::vec ZteS = arma::zeros(beta.n_elem);
  arma::vec ZteY = arma::zeros(beta.n_elem);
  int m;
  
  parameter_new(arma::span(0, alpha.n_elem-1)) = alpha;
  parameter_new(arma::span(alpha.n_elem, alpha.n_elem + beta.n_elem - 1)) = beta;
  parameter_new(alpha.n_elem + beta.n_elem) = gamma;
  
  arma::vec E1 = ES12.col(0);
  double A6 = sum(ES12.col(1));
  double temp0 = A6 - arma::as_scalar(E1.t() * E1);
  
  parameter_old = parameter_new;
  
  eS = E1 - X * alpha;
  eY = Y - gamma * E1 - Z * beta;
  SteS = A6 - arma::as_scalar(E1.t() * X * alpha);
  SteY = arma::as_scalar(E1.t() * (Y - Z*beta)) - gamma * A6;
  XteS = X.t() * eS;
  XteY = X.t() * eY;
  ZteS = Z.t() * eS;
  ZteY = Z.t() * eY;
  SS = arma::as_scalar(eS.t() * eS) + temp0;
  YY = arma::as_scalar(eY.t() * eY) + temp0 * pow(gamma, 2);
  SY = arma::as_scalar(eS.t() * eY) - temp0 * gamma;
  
  HessianMat(arma::span(0, alpha.n_elem-1), arma::span(0, alpha.n_elem-1))
    = 2*YY*X.t()*X - 2*XteY*XteY.t()
    - (2*SY*XteY - 2*YY*XteS) * (2*SY*XteY - 2*YY*XteS).t() / (SS*YY-SY*SY);
  HessianMat(arma::span(0, alpha.n_elem-1), arma::span(alpha.n_elem, alpha.n_elem + beta.n_elem - 1))
    = 4*XteS*ZteY.t() - 2*XteY*ZteS.t() - 2*SY*X.t()*Z
    - (2*SY*XteY - 2*YY*XteS) * (2*SY*ZteS - 2*SS*ZteY).t() / (SS*YY-SY*SY);
  HessianMat(arma::span(0, alpha.n_elem-1), alpha.n_elem + beta.n_elem)
    = 4*SteY*XteS - 2*SteS*XteY - 2*SY*X.t()*E1
    - (2*SY*SteS - 2*SS*SteY) * (2*SY*XteY - 2*YY*XteS) / (SS*YY-SY*SY);
        
  HessianMat(arma::span(alpha.n_elem, alpha.n_elem + beta.n_elem - 1), arma::span(0, alpha.n_elem-1))
    = 4*ZteY*XteS.t() - 2*ZteS*XteY.t() - 2*SY*Z.t()*X
    - (2*SY*ZteS - 2*SS*ZteY) * (2*SY*XteY - 2*YY*XteS).t() / (SS*YY-SY*SY);
  HessianMat(arma::span(alpha.n_elem, alpha.n_elem + beta.n_elem - 1), arma::span(alpha.n_elem, alpha.n_elem + beta.n_elem - 1))
    = 2*SS*Z.t()*Z - 2*ZteS*ZteS.t()
    - (2*SY*ZteS - 2*SS*ZteY) * (2*SY*ZteS - 2*SS*ZteY).t() / (SS*YY-SY*SY);
  HessianMat(arma::span(alpha.n_elem, alpha.n_elem + beta.n_elem - 1), alpha.n_elem + beta.n_elem)
    = 2*SS*Z.t()*E1 - 2*SteS*ZteS
    - (2*SY*SteS - 2*SS*SteY) * (2*SY*ZteS - 2*SS*ZteY) / (SS*YY-SY*SY);
              
  HessianMat(alpha.n_elem + beta.n_elem, arma::span(0, alpha.n_elem-1))
    = 4*SteY*XteS.t() - 2*SteS*XteY.t() - 2*SY*E1.t()*X
    - (2*SY*SteS - 2*SS*SteY) * (2*SY*XteY - 2*YY*XteS).t() / (SS*YY-SY*SY);
  HessianMat(alpha.n_elem + beta.n_elem, arma::span(alpha.n_elem, alpha.n_elem + beta.n_elem - 1))
    = 2*SS*E1.t()*Z - 2*SteS*ZteS.t()
    - (2*SY*SteS - 2*SS*SteY) * (2*SY*ZteS - 2*SS*ZteY).t() / (SS*YY-SY*SY);
  HessianMat(alpha.n_elem + beta.n_elem, alpha.n_elem + beta.n_elem)
    = 2*SS*A6 - 2*pow(SteS, 2)
    - (2*SY*SteS - 2*SS*SteY) * (2*SY*SteS - 2*SS*SteY) / (SS*YY-SY*SY);
                    
  grad(arma::span(0, alpha.n_elem-1)) = 2*SY*XteY - 2*YY*XteS;
  grad(arma::span(alpha.n_elem, alpha.n_elem + beta.n_elem - 1)) = 2*SY*ZteS - 2*SS*ZteY;
  grad(alpha.n_elem + beta.n_elem) = 2*SY*SteS - 2*SS*SteY;
                    
  parameter_change = HessianMat.i() * grad;
                    
  parameter_new_temp = parameter_old - parameter_change;
  alpha_new_temp = parameter_new_temp(arma::span(0, alpha.n_elem-1));
  beta_new_temp = parameter_new_temp(arma::span(alpha.n_elem, alpha.n_elem + beta.n_elem - 1));
  gamma_new_temp = parameter_new_temp(alpha.n_elem + beta.n_elem);
  
  eS_temp = E1 - X * alpha_new_temp;
  eY_temp = Y - gamma_new_temp * E1 - Z * beta_new_temp;
  SS_temp = arma::as_scalar(eS_temp.t() * eS_temp) + temp0;
  YY_temp = arma::as_scalar(eY_temp.t() * eY_temp) + temp0 * pow(gamma_new_temp, 2);
  SY_temp = arma::as_scalar(eY_temp.t() * eS_temp) - temp0 * gamma_new_temp;
                    
  m = 0;
  while( (SS*YY-SY*SY) < (SS_temp*YY_temp-SY_temp*SY_temp) ){
    m = m + 1;
                      
    parameter_new_temp = parameter_old - parameter_change / pow(2, m);
    alpha_new_temp = parameter_new_temp(arma::span(0, alpha.n_elem-1));
    beta_new_temp = parameter_new_temp(arma::span(alpha.n_elem, alpha.n_elem + beta.n_elem - 1));
    gamma_new_temp = parameter_new_temp(alpha.n_elem + beta.n_elem);
    
    eS_temp = E1 - X * alpha_new_temp;
    eY_temp = Y - gamma_new_temp * E1 - Z * beta_new_temp;
    SS_temp = arma::as_scalar(eS_temp.t() * eS_temp) + temp0;
    YY_temp = arma::as_scalar(eY_temp.t() * eY_temp) + temp0 * pow(gamma_new_temp, 2);
    SY_temp = arma::as_scalar(eY_temp.t() * eS_temp) - temp0 * gamma_new_temp;
  }
  
  //cout << "m=" << m << endl;
  return parameter_new_temp;
}



// sigmaSsquare update
double sigmaSsquare_update(arma::vec alpha, arma::mat X, arma::mat ES12){
  arma::mat eS = ES12.col(0) - X * alpha;
  double temp1 = arma::as_scalar(eS.t() * eS);
  double temp2 = sum(ES12.col(1)) - arma::as_scalar(ES12.col(0).t() * ES12.col(0));
  return (temp1 + temp2) / X.n_rows;
}



// sigmaYsquare update
double sigmaYsquare_update(arma::vec Y, arma::mat Z, arma::mat ES12, arma::vec beta, double gamma){
  arma::mat eY = Y - gamma * ES12.col(0) - Z * beta;
  double temp1 = arma::as_scalar(eY.t() * eY);
  double temp2 = sum(ES12.col(1)) - arma::as_scalar(ES12.col(0).t() * ES12.col(0));
  return (temp1 + gamma * gamma * temp2) / Z.n_rows;
}



// rho update
double rho_update(arma::mat X, arma::vec Y, arma::mat Z, arma::mat ES12,
                  arma::vec alpha, arma::vec beta, double gamma, double sigmaSsquare, double sigmaYsquare){
  arma::mat eY = Y - gamma * ES12.col(0) - Z * beta;
  arma::mat eS = ES12.col(0) - X * alpha;
  double temp1 = arma::as_scalar(eS.t() * eY);
  double temp2 = sum(ES12.col(1)) - arma::as_scalar(ES12.col(0).t() * ES12.col(0));
  return (temp1 - gamma * temp2) / ( Z.n_rows * sqrt(sigmaSsquare) * sqrt(sigmaYsquare) );
}



// Louis Formula: Computaion of Ic
arma::mat Info_complete(arma::mat X, arma::vec Y, arma::mat Z, arma::mat ES12,
                        arma::vec alpha, arma::vec beta, double gamma, double sigmaSsquare, double sigmaYsquare, double rho){
  int n = X.n_rows;
  arma::mat Ic = arma::zeros(X.n_cols + Z.n_cols + 4, X.n_cols + Z.n_cols + 4);
  
  arma::vec E1 = ES12.col(0);
  arma::vec eS = E1 - X * alpha;
  arma::vec eY = Y - gamma * E1 - Z * beta;
  //double temp0 = sum(ES12.col(1)) - arma::as_scalar(E1.t() * E1);
  double sigmaS = pow(sigmaSsquare, 0.5);
  double sigmaY = pow(sigmaYsquare, 0.5);
  arma::vec XteS = X.t() * eS;
  arma::vec XteY = X.t() * eY;
  arma::vec ZteS = Z.t() * eS;
  arma::vec ZteY = Z.t() * eY;
  double SteS = sum(ES12.col(1)) - arma::as_scalar(E1.t() * X * alpha);
  double SteY = arma::as_scalar(E1.t() * (Y - Z*beta)) - gamma * sum(ES12.col(1));
  //double SS = arma::as_scalar(eS.t() * eS) + temp0;
  //double YY = arma::as_scalar(eY.t() * eY) + pow(gamma, 2) * temp0;
  //double SY = arma::as_scalar(eS.t() * eY) - gamma * temp0;
  
  Ic(arma::span(0, X.n_cols-1), arma::span(0, X.n_cols-1)) = X.t()*X / ( (1-pow(rho,2)) * sigmaSsquare );
  Ic(arma::span(0, X.n_cols-1), arma::span(X.n_cols, X.n_cols+Z.n_cols-1)) = -rho*X.t()*Z / ( (1-pow(rho,2))*sigmaS*sigmaY );
  Ic(arma::span(0, X.n_cols-1), X.n_cols+Z.n_cols) = -rho*X.t()*E1 / ( (1-pow(rho,2))*sigmaS*sigmaY );
  Ic(arma::span(0, X.n_cols-1), X.n_cols+Z.n_cols+1) = XteS / ( (1-rho*rho) * pow(sigmaSsquare, 2) ) - rho*XteY / (2 * (1-rho*rho) * pow(sigmaS, 3) * sigmaY);
  Ic(arma::span(0, X.n_cols-1), X.n_cols+Z.n_cols+2) = -rho*XteY / (2 * (1-rho*rho) * sigmaS * pow(sigmaY, 3));
  Ic(arma::span(0, X.n_cols-1), X.n_cols+Z.n_cols+3) = -2*rho*XteS / ((1-rho*rho) * (1-rho*rho) * sigmaSsquare) + (1+rho*rho)*XteY / ((1-rho*rho) * (1-rho*rho) * sigmaS * sigmaY);
  
  Ic(arma::span(X.n_cols, X.n_cols + Z.n_cols - 1), arma::span(0, X.n_cols-1)) = -rho*Z.t()*X / ( (1-rho*rho)*sigmaS*sigmaY );
  Ic(arma::span(X.n_cols, X.n_cols + Z.n_cols - 1), arma::span(X.n_cols, X.n_cols+Z.n_cols-1)) = Z.t()*Z / ( (1-rho*rho)*sigmaYsquare );
  Ic(arma::span(X.n_cols, X.n_cols + Z.n_cols - 1), X.n_cols+Z.n_cols) = Z.t()*E1 / ( (1-rho*rho)*sigmaYsquare );
  Ic(arma::span(X.n_cols, X.n_cols + Z.n_cols - 1), X.n_cols+Z.n_cols+1) = -rho*ZteS / (2 * (1-rho*rho) * pow(sigmaS, 3) * sigmaY);
  Ic(arma::span(X.n_cols, X.n_cols + Z.n_cols - 1), X.n_cols+Z.n_cols+2) = ZteY / ( (1-rho*rho) * pow(sigmaYsquare, 2) ) - rho*ZteS / (2 * (1-rho*rho) * sigmaS * pow(sigmaY, 3));
  Ic(arma::span(X.n_cols, X.n_cols + Z.n_cols - 1), X.n_cols+Z.n_cols+3) = -2*rho*ZteY / ((1-rho*rho) * (1-rho*rho) * sigmaYsquare) + (1+rho*rho)*ZteS / ((1-rho*rho) * (1-rho*rho) * sigmaS * sigmaY);
  
  Ic(X.n_cols + Z.n_cols, arma::span(0, X.n_cols-1)) = -rho*E1.t()*X / ( (1-rho*rho)*sigmaS*sigmaY );
  Ic(X.n_cols + Z.n_cols, arma::span(X.n_cols, X.n_cols+Z.n_cols-1)) = E1.t()*Z / ( (1-rho*rho)*sigmaYsquare );
  Ic(X.n_cols + Z.n_cols, X.n_cols+Z.n_cols) = sum(ES12.col(1)) / ( (1-rho*rho)*sigmaYsquare );
  Ic(X.n_cols + Z.n_cols, X.n_cols+Z.n_cols+1) = -rho*SteS / (2 * (1-rho*rho) * pow(sigmaS, 3) * sigmaY);
  Ic(X.n_cols + Z.n_cols, X.n_cols+Z.n_cols+2) = SteY / ( (1-rho*rho) * pow(sigmaYsquare, 2) ) - rho*SteS / (2 * (1-rho*rho) * sigmaS * pow(sigmaY, 3));
  Ic(X.n_cols + Z.n_cols, X.n_cols+Z.n_cols+3) = -2*rho*SteY / ((1-rho*rho) * (1-rho*rho) * sigmaYsquare) + (1+rho*rho)*SteS / ((1-rho*rho) * (1-rho*rho) * sigmaS * sigmaY);
  
  Ic(X.n_cols + Z.n_cols+1, arma::span(0, X.n_cols-1)) = XteS.t() / ( (1-rho*rho) * pow(sigmaSsquare, 2) ) - rho*XteY.t() / (2 * (1-rho*rho) * pow(sigmaS, 3) * sigmaY);
  Ic(X.n_cols + Z.n_cols+1, arma::span(X.n_cols, X.n_cols+Z.n_cols-1)) = -rho*ZteS.t() / (2 * (1-rho*rho) * pow(sigmaS, 3) * sigmaY);
  Ic(X.n_cols + Z.n_cols+1, X.n_cols+Z.n_cols) = Ic(X.n_cols + Z.n_cols, X.n_cols+Z.n_cols+1);
  Ic(X.n_cols + Z.n_cols+1, X.n_cols+Z.n_cols+1) = n*(2.0-rho*rho) / ( 4*(1-rho*rho)*pow(sigmaSsquare, 2) );
  Ic(X.n_cols + Z.n_cols+1, X.n_cols+Z.n_cols+2) = -n*rho*rho / ( 4*(1-rho*rho)*sigmaSsquare*sigmaYsquare );
  Ic(X.n_cols + Z.n_cols+1, X.n_cols+Z.n_cols+3) = -n*rho / ( 2*(1-rho*rho)*sigmaSsquare );
  
  Ic(X.n_cols + Z.n_cols+2, arma::span(0, X.n_cols-1)) = -rho*XteY.t() / (2 * (1-rho*rho) * sigmaS * pow(sigmaY, 3));
  Ic(X.n_cols + Z.n_cols+2, arma::span(X.n_cols, X.n_cols+Z.n_cols-1)) = ZteY.t() / ( (1-rho*rho) * pow(sigmaYsquare, 2) ) - rho*ZteS.t() / (2 * (1-rho*rho) * sigmaS * pow(sigmaY, 3));
  Ic(X.n_cols + Z.n_cols+2, X.n_cols+Z.n_cols) = Ic(X.n_cols + Z.n_cols, X.n_cols+Z.n_cols+2);
  Ic(X.n_cols + Z.n_cols+2, X.n_cols+Z.n_cols+1) = Ic(X.n_cols + Z.n_cols+1, X.n_cols+Z.n_cols+2);
  Ic(X.n_cols + Z.n_cols+2, X.n_cols+Z.n_cols+2) = n*(2.0-rho*rho) / ( 4*(1-rho*rho)*pow(sigmaYsquare, 2) );
  Ic(X.n_cols + Z.n_cols+2, X.n_cols+Z.n_cols+3) = -n*rho / ( 2*(1-rho*rho)*sigmaYsquare );
  
  Ic(X.n_cols + Z.n_cols+3, arma::span(0, X.n_cols-1)) = -2*rho*XteS.t() / ((1-rho*rho) * (1-rho*rho) * sigmaSsquare) + (1+rho*rho)*XteY.t() / ((1-rho*rho) * (1-rho*rho) * sigmaS * sigmaY);
  Ic(X.n_cols + Z.n_cols+3, arma::span(X.n_cols, X.n_cols+Z.n_cols-1)) = -2*rho*ZteY.t() / ((1-rho*rho) * (1-rho*rho) * sigmaYsquare) + (1+rho*rho)*ZteS.t() / ((1-rho*rho) * (1-rho*rho) * sigmaS * sigmaY);
  Ic(X.n_cols + Z.n_cols+3, X.n_cols+Z.n_cols) = Ic(X.n_cols + Z.n_cols, X.n_cols+Z.n_cols+3);
  Ic(X.n_cols + Z.n_cols+3, X.n_cols+Z.n_cols+1) = Ic(X.n_cols + Z.n_cols+1, X.n_cols+Z.n_cols+3);
  Ic(X.n_cols + Z.n_cols+3, X.n_cols+Z.n_cols+2) = Ic(X.n_cols + Z.n_cols+2, X.n_cols+Z.n_cols+3);
  Ic(X.n_cols + Z.n_cols+3, X.n_cols+Z.n_cols+3) = n*(1.0+rho*rho) / pow((1-rho*rho), 2);
  
  //cout << "C" << "  ";
  return Ic;
}



// Louis Formula: Computaion of E[I(mis|obs)]
arma::mat Info_mis_given_obs(arma::vec Y, arma::vec S, arma::vec R, arma::mat Z, arma::mat X,
                             arma::vec alpha, arma::vec beta, double gamma, double sigmaSsquare, double sigmaYsquare, double rho,
                             double a, arma::vec b, arma::mat ES12){
  int n = R.n_elem;
  int npara = alpha.n_elem + beta.n_elem + 4;
  
  arma::mat temp = arma::zeros(npara, npara);
  arma::mat temp33 = arma::zeros(3, 3);
  arma::vec temp31 = arma::zeros(3);
  
  arma::mat Vi = arma::zeros(npara, 3);
  double L0i, U0i;
  
  double sigmaS = pow(sigmaSsquare, 0.5);
  double sigmaY = pow(sigmaYsquare, 0.5);
  double ealpha = 0;
  double ebeta = 0;
  
  for(int i=0; i<n; i++){
    if(R(i) != 0){
      temp31(0) = 1;
      temp31(1) = ES12(i, 0);
      temp31(2) = ES12(i, 1);
      
      temp33(0, 0) = 1;
      temp33(0, 1) = ES12(i, 0);
      temp33(1, 0) = ES12(i, 0);
      temp33(0, 2) = ES12(i, 1);
      temp33(1, 1) = ES12(i, 1);
      temp33(2, 0) = ES12(i, 1);
      if (R(i) == 3) {
        temp33(1, 2) = pow(b(i), 3) + 3 * a * b(i);
        temp33(2, 1) = temp33(1, 2);
        temp33(2, 2) = pow(b(i), 4) + 6 * pow(b(i), 2) * a + 3 * pow(a, 2);
      }
      else if (R(i) == 1){
        L0i = (S(i) - b(i)) / pow(a, 0.5);
        temp33(1, 2) = pow(b(i), 3) + 3*a*b(i) - arma::normpdf(L0i) / arma::normcdf(L0i) * (3*pow(b(i),2)*pow(a, 0.5) + 2*pow(a, 1.5) + 3*a*b(i)*L0i + pow(a, 1.5)*pow(L0i,2));
        temp33(2, 1) = temp33(1, 2);
        temp33(2, 2) = pow(b(i), 4) + 6*pow(b(i), 2)*a + 3*pow(a, 2) - arma::normpdf(L0i) / arma::normcdf(L0i) * (4*pow(b(i),3)*pow(a, 0.5) + 8*b(i)*pow(a, 1.5) + (3*pow(a, 2) + 6*a*pow(b(i), 2))*L0i + 4*b(i)*pow(a, 1.5)*pow(L0i,2) + pow(a,2)*pow(L0i,3));
      }
      else{
        U0i = (-S(i) + b(i)) / pow(a, 0.5);
        temp33(1, 2) = pow(b(i), 3) + 3*a*b(i) + arma::normpdf(U0i) / arma::normcdf(U0i) * (3*pow(b(i),2)*pow(a, 0.5) + 2*pow(a, 1.5) - 3*a*b(i)*U0i + pow(a, 1.5)*pow(U0i,2));
        temp33(2, 1) = temp33(1, 2);
        temp33(2, 2) = pow(b(i), 4) + 6*pow(b(i), 2)*a + 3*pow(a, 2) + arma::normpdf(U0i) / arma::normcdf(U0i) * (4*pow(b(i),3)*pow(a, 0.5) + 8*b(i)*pow(a, 1.5) - (3*pow(a, 2) + 6*a*pow(b(i), 2))*U0i + 4*b(i)*pow(a, 1.5)*pow(U0i,2) - pow(a,2)*pow(U0i,3));
      }
      
      // compute Vi
      ealpha = arma::as_scalar(X.row(i) * alpha);
      ebeta = Y(i) - arma::as_scalar(Z.row(i) * beta);
      
      Vi(arma::span(0, X.n_cols-1), 0) = -ealpha*X.row(i).t() / ((1-rho*rho)*sigmaSsquare)
        - rho*ebeta*X.row(i).t() / ((1-rho*rho)*sigmaS*sigmaY);
      Vi(arma::span(0, X.n_cols-1), 1) = X.row(i).t() / ((1-rho*rho)*sigmaSsquare)
        + rho*gamma*X.row(i).t() / ((1-rho*rho)*sigmaS*sigmaY);
      
      Vi(arma::span(X.n_cols, X.n_cols+Z.n_cols-1), 0) = ebeta * Z.row(i).t() / ((1-rho*rho)*sigmaYsquare)
        + rho*ealpha*Z.row(i).t() / ((1-rho*rho)*sigmaS*sigmaY);
      Vi(arma::span(X.n_cols, X.n_cols+Z.n_cols-1), 1) = -gamma * Z.row(i).t() / ((1-rho*rho)*sigmaYsquare)
        - rho*Z.row(i).t() / ((1-rho*rho)*sigmaS*sigmaY);
      
      Vi(X.n_cols+Z.n_cols, 1) = ebeta / ((1-rho*rho)*sigmaYsquare) + rho * ealpha / ((1-rho*rho)*sigmaS*sigmaY);
      Vi(X.n_cols+Z.n_cols, 2) = -gamma / ((1-rho*rho)*sigmaYsquare) - rho / ((1-rho*rho)*sigmaS*sigmaY);
      
      Vi(X.n_cols+Z.n_cols+1, 0)
        = -1.0 / (2*sigmaSsquare)
        + pow(ealpha, 2) / (2 * pow(sigmaSsquare, 2) * (1-rho*rho))
        + rho * ealpha * ebeta / (2*(1-rho*rho)*pow(sigmaSsquare, 1.5)*sigmaY);
      Vi(X.n_cols+Z.n_cols+1, 1)
        = -ealpha / ( (1-rho*rho) * pow(sigmaSsquare, 2) )
        - rho * (ebeta + gamma*ealpha) / (2*(1-rho*rho)*pow(sigmaSsquare, 1.5)*sigmaY);
      Vi(X.n_cols+Z.n_cols+1, 2)
        = 1.0 / ( 2*(1-rho*rho)*pow(sigmaSsquare, 2) )
        + rho*gamma / ( 2*(1-rho*rho)*pow(sigmaSsquare, 1.5)*sigmaY );
            
      Vi(X.n_cols+Z.n_cols+2, 0)
        = -1.0 / (2*sigmaYsquare)
        + pow(ebeta, 2) / (2*pow(sigmaYsquare, 2)*(1-rho*rho))
        + rho * ealpha * ebeta / (2*(1-rho*rho)*pow(sigmaYsquare, 1.5)*sigmaS);
      Vi(X.n_cols+Z.n_cols+2, 1)
        = -gamma * ebeta / ( (1-rho*rho) * pow(sigmaYsquare, 2) )
        - rho * (ebeta + gamma*ealpha) / (2*(1-rho*rho)*pow(sigmaYsquare, 1.5)*sigmaS);
      Vi(X.n_cols+Z.n_cols+2, 2)
        = pow(gamma, 2) / (2*pow(sigmaYsquare, 2)*(1-rho*rho))
        + rho*gamma / ( 2*(1-rho*rho)*pow(sigmaYsquare, 1.5)*sigmaS );
                  
      Vi(X.n_cols+Z.n_cols+3, 0)
        = rho / (1-rho*rho)
        - rho*pow(ealpha, 2) / (sigmaSsquare*(1-rho*rho)*(1-rho*rho))
        - rho*pow(ebeta, 2) / (sigmaYsquare*(1-rho*rho)*(1-rho*rho))
        - (1+rho*rho)*ealpha*ebeta / (sigmaS*sigmaY*(1-rho*rho)*(1-rho*rho));
      Vi(X.n_cols+Z.n_cols+3, 1)
        = 2*rho*ealpha / (sigmaSsquare*(1-rho*rho)*(1-rho*rho))
        + 2*rho*gamma*ebeta / (sigmaYsquare*(1-rho*rho)*(1-rho*rho))
        + (1+rho*rho)*(ebeta + gamma*ealpha) / (sigmaS*sigmaY*(1-rho*rho)*(1-rho*rho));
      Vi(X.n_cols+Z.n_cols+3, 2)
        = -rho / (sigmaSsquare*(1-rho*rho)*(1-rho*rho))
        - rho*gamma*gamma / (sigmaYsquare*(1-rho*rho)*(1-rho*rho))
        - (1+rho*rho)*gamma / (sigmaS*sigmaY*(1-rho*rho)*(1-rho*rho));
                        
      temp = temp + Vi * (temp33 - temp31 * temp31.t()) * Vi.t();
    }
  }
  //cout << "D" << "  ";
  return temp;
}



// Computation of estimated covariance matrix
arma::mat est_cov_Louis(arma::vec Y, arma::vec S, arma::vec R, arma::mat Z, arma::mat X,
                        arma::vec alpha, arma::vec beta, double gamma, double sigmaSsquare, double sigmaYsquare, double rho,
                        double a, arma::vec b, arma::mat ES12){
  arma::mat mat1= Info_complete(X, Y, Z, ES12, alpha, beta, gamma, sigmaSsquare, sigmaYsquare, rho);
  arma::mat mat2= Info_mis_given_obs(Y, S, R, Z, X, alpha, beta, gamma, sigmaSsquare, sigmaYsquare, rho, a, b, ES12);
  arma::mat mat3= (mat1 - mat2).i();
  //cout << "E" << "  ";
  return mat3;
}



struct para_est_EM{
  arma::vec alpha;
  arma::vec beta;
  double gamma;
  double sigmaSsquare;
  double sigmaYsquare;
  double rho;
  arma::mat cov;
  int num_iter;
};



struct para_est_EM MR_MLE_rho(arma::vec Y, arma::vec S, arma::vec R, arma::mat Z, arma::mat X,
                                                                                        arma::vec alpha_init, arma::vec beta_init, double gamma_init, double sigmaSsquare_init, double sigmaYsquare_init, double rho_init,
                                                                                        double epsilon){
  struct para_est_EM mle;
  int n = Z.n_rows;             // number of observations
  int npara = alpha_init.n_elem + beta_init.n_elem + 4;
  
  //cout << "A" << "  ";
  arma::vec alpha_old = alpha_init;
  arma::vec alpha_new = alpha_init;
  arma::vec alpha_change = arma::zeros(alpha_init.n_elem);
  arma::vec beta_old = beta_init;
  arma::vec beta_new = beta_init;
  arma::vec beta_change = arma::zeros(beta_init.n_elem);
  double gamma_old = gamma_init;
  double gamma_new = gamma_init;
  double gamma_change = 0;
  double sigmaSsquare_old = sigmaSsquare_init;
  double sigmaSsquare_new = sigmaSsquare_init;
  double sigmaSsquare_change = 0;
  double sigmaYsquare_old = sigmaYsquare_init;
  double sigmaYsquare_new = sigmaYsquare_init;
  double sigmaYsquare_change = 0;
  double rho_old = rho_init;
  double rho_new = rho_init;
  double rho_change = 0;
  
  arma::mat est_covariance_mat = arma::zeros(npara, npara);
  arma::vec alpha_beta_gamma = arma::zeros(npara - 3);
  
  //cout << "B" << "  ";
  double a_new = 0;
  arma::vec b_new = arma::zeros(n);
  arma::mat ES12_new = arma::zeros(n, 2);
  
  int count = 0;
  double dist = 1;
  
  while((dist > epsilon) & (count < 1000)){
    count++;
    //cout << "iteration:" << count << "  ";
    
    // E-Step
    a_new = a_update(gamma_new, sigmaSsquare_new, sigmaYsquare_new, rho_new);
    //cout << "1" << "  ";
    b_new = b_update(a_new, alpha_new, beta_new, gamma_new, sigmaSsquare_new, sigmaYsquare_new, rho_new, X, Z, Y, R);
    //cout << "2" << "  ";
    ES12_new = ES12_update(a_new, b_new, R, S);
    //cout << "3" << "  ";
    
    // M-step
    alpha_old = alpha_new;
    beta_old = beta_new;
    gamma_old = gamma_new;
    sigmaSsquare_old = sigmaSsquare_new;
    sigmaYsquare_old = sigmaYsquare_new;
    rho_old = rho_new;
    //cout << "4" << "  ";
    
    alpha_beta_gamma = alpha_beta_gamma_update(X, Y, Z, ES12_new, alpha_old, beta_old, gamma_old);
    alpha_new = alpha_beta_gamma.subvec(0, X.n_cols - 1);
    beta_new = alpha_beta_gamma.subvec(X.n_cols, X.n_cols + Z.n_cols - 1);
    gamma_new = alpha_beta_gamma(X.n_cols + Z.n_cols);
    
    alpha_change = alpha_new - alpha_old;
    beta_change = beta_new - beta_old;
    gamma_change = gamma_new - gamma_old;
    //cout << alpha_new << "  " << beta_new << "  " << gamma_new << "  ";
    
    sigmaSsquare_new = sigmaSsquare_update(alpha_new, X, ES12_new);
    //cout << sigmaSsquare_new << "  ";
    sigmaSsquare_change = sigmaSsquare_new - sigmaSsquare_old;
    
    sigmaYsquare_new = sigmaYsquare_update(Y, Z, ES12_new, beta_new, gamma_new);
    //cout << sigmaYsquare_new << "  ";
    sigmaYsquare_change = sigmaYsquare_new - sigmaYsquare_old;
    
    rho_new = rho_update(X, Y, Z, ES12_new, alpha_new, beta_new, gamma_new, sigmaSsquare_new, sigmaYsquare_new);
    //cout << rho_new << endl ;
    rho_change = rho_new - rho_old;
    
    dist = sqrt(arma::as_scalar(alpha_change.t() * alpha_change) +
      arma::as_scalar(beta_change.t() * beta_change) + pow(gamma_change, 2) +
      pow(sigmaSsquare_change, 2) + pow(sigmaYsquare_change, 2) + pow(rho_change, 2));
    
    //cout << alpha_new(0) << " "<< gamma_new << " " << dist << endl ;
  }
  
  // value of ES12 at MLE
  a_new = a_update(gamma_new, sigmaSsquare_new, sigmaYsquare_new, rho_new);
  b_new = b_update(a_new, alpha_new, beta_new, gamma_new, sigmaSsquare_new, sigmaYsquare_new, rho_new, X, Z, Y, R);
  ES12_new = ES12_update(a_new, b_new, R, S);
  
  est_covariance_mat = est_cov_Louis(Y, S, R, Z, X, alpha_new, beta_new, gamma_new, sigmaSsquare_new, sigmaYsquare_new, rho_new, a_new, b_new, ES12_new);
  
  mle.alpha = alpha_new;
  mle.beta = beta_new;
  mle.gamma = gamma_new;
  mle.sigmaSsquare = sigmaSsquare_new;
  mle.sigmaYsquare = sigmaYsquare_new;
  mle.rho = rho_new;
  mle.cov = est_covariance_mat;
  mle.num_iter = count;
  
  return mle;
}




// [[Rcpp::export]]
List MR_estimation_rho(arma::vec Y, arma::vec S, arma::vec R, arma::mat Z, arma::mat X,
                       arma::vec alpha_init, arma::vec beta_init, double gamma_init, double sigmaSsquare_init, double sigmaYsquare_init, double rho_init,
                       double epsilon){
  struct para_est_EM mle = MR_MLE_rho(Y, S, R, Z, X, alpha_init, beta_init, gamma_init, sigmaSsquare_init, sigmaYsquare_init, rho_init, epsilon);
  return List::create(mle.alpha, mle.beta, mle.gamma, mle.sigmaSsquare, mle.sigmaYsquare, mle.rho, mle.cov, mle.num_iter);
}

