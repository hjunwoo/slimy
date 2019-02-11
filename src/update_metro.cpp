#include <RcppEigen.h>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

double zprob(const Eigen::VectorXd &yx, double a,double b, double hv,
             const Eigen::MatrixXd &xpa, const Eigen::MatrixXd &xtx,
             const Eigen::VectorXd &cw, double m,double sg){

  int nsample = yx.size();
  double v = 2*b + yx.squaredNorm();
  int npa = xpa.cols();
  if(npa > 0){
    Eigen::VectorXd xy = xpa.transpose()*yx;
    v -= hv*xy.dot(xtx*xy);
  }
  double z = -(nsample/2 + a)*log(v);
  if(npa > 0){
    Eigen::FullPivLU<Eigen::MatrixXd> lu(xtx);
    z += 0.5*log(lu.determinant());
  }
  Eigen::VectorXd e = sg*yx +
    Eigen::VectorXd::Constant(yx.size(),m);
  double lkh = e.dot(cw);
  for(int j=0; j<nsample; j++) lkh -= exp(e(j));
  z += lkh;

  return(z);
}
//
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd update_metro(const Eigen::MatrixXd &ci,
                             const Eigen::MatrixXd &xi,
      const Rcpp::IntegerVector &W, const Rcpp::List &hyper,
      const Rcpp::List &po,
      const Eigen::MatrixXd &A,
      const Rcpp::NumericVector &dy,
      const Rcpp::IntegerVector &updaten,
      const Rcpp::IntegerVector &seed){

   const gsl_rng_type *T;
   gsl_rng *r;
   gsl_rng_env_setup();
   T = gsl_rng_default;
   r = gsl_rng_alloc(T);
   unsigned long Seed = seed[0];
   gsl_rng_set(r, Seed);
   int upn = updaten[0];

   double Dy = dy[0];
   int nsample = ci.rows();
   int p = ci.cols();
   int q = W.size();
   double b = hyper["b"];
   double a = hyper["a"];
   double hv = hyper["v"];
   Rcpp::NumericVector mu = po["mu"];
   Rcpp::NumericVector sigma = po["sigma"];

   Eigen::MatrixXd xi2 = xi;

   int *ksamp = new int[upn];
   int *sid = new int[nsample];
   for(int k=0; k<nsample; k++) sid[k] = k;
   for(int i=0; i<q; i++){
     int idx =  W[i];
     std::vector<int> pa;
     for(int j=0; j<p; j++) if(A(j,idx)==1) pa.push_back(j);
     int npa = pa.size();
     Eigen::VectorXd yx = xi.block(0,idx,nsample,1);  // xi[,idx]
     Eigen::MatrixXd xpa(nsample,npa);
     for(int i=0; i<npa; i++) for(int z=0; z<nsample; z++)
       xpa(z,i) = xi(z,pa[i]);
     Eigen::MatrixXd xtx;
     if(npa > 0){
       Eigen::MatrixXd I;
       I.setIdentity(npa,npa);
       xtx = I + hv*xpa.transpose()*xpa;
       xtx = xtx.ldlt().solve(I);
     }
     double m = mu[idx];
     double sg = sigma[idx];
     gsl_ran_choose(r, ksamp, upn, sid, nsample, sizeof(int));   // choose upn sample from sample id's
     for(int u=0; u<upn; u++){
       int k = ksamp[u];
       double x0 = zprob(yx, a, b, hv, xpa, xtx,ci.col(idx), m, sg);
       double delta = gsl_ran_gaussian(r, Dy);
       yx[k] += delta;
       double x1 = zprob(yx, a, b, hv, xpa, xtx, ci.col(idx), m, sg);
       bool accept = false;
       if(x1 > x0) accept = true;
       else accept = (exp(x1-x0) > gsl_rng_uniform(r));
       if(accept)
         xi2(k, idx) = yx[k];
       else
         yx[k] -= delta;
     }
  }
  delete[] ksamp;
  delete[] sid;
  gsl_rng_free(r);

  return xi2;
}
