#include <RcppEigen.h>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

double sample(const std::vector<double> &gridy, const std::vector<double> &prob,
           gsl_rng *r){

   double x = gsl_rng_uniform(r);
   double sum = 0;
   int i = 0;
   for(i=0; i<int(gridy.size()); i++){
     sum += prob[i];
     if(sum > x) break;
   }
   return(gridy[i]);
}

//
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd update_field(const Eigen::MatrixXd &ci,
                             const Eigen::MatrixXd &xi,
      const Rcpp::IntegerVector &W, const Rcpp::List &hyper,
      const Rcpp::List &po,
      const Eigen::MatrixXd &A, const Rcpp::NumericVector &dmax,
      const Rcpp::NumericVector &dy, const Rcpp::IntegerVector &updaten,
      const Rcpp::IntegerVector &seed){

   const gsl_rng_type *T;
   gsl_rng *r;
   gsl_rng_env_setup();
   T = gsl_rng_default;
   r = gsl_rng_alloc(T);
   unsigned long Seed = seed[0];
   gsl_rng_set(r, Seed);
   int upn = updaten[0];

   double Dmax = dmax[0];
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

// std::cout << ci(0,0) << std::endl;
// std::cout << ci(999,4) << std::endl;
// exit(1);

   std::vector<double> gridy;
   for(double x=-Dmax; x <=Dmax; x+= Dy) gridy.push_back(x);

   int *ksamp = new int[upn];
   int *sid = new int[nsample];
   for(int k=0; k<nsample; k++) sid[k] = k;
   for(int i=0; i<q; i++){
     int idx =  W[i];
     std::vector<int> pa;
     for(int j=0; j<p; j++) if(A(j,idx)==1) pa.push_back(j);
     int npa = pa.size();
     Eigen::VectorXd yx = xi.block(0,idx,nsample,1);  // xi[,idx]
     gsl_ran_choose(r, ksamp, upn, sid, nsample, sizeof(int));   // choose upn sample from sample id's
     for(int u=0; u<upn; u++){
       int k = ksamp[u];
//k = 445;
       std::vector<double> prob;
       double zmax=0;
       for(int iy=0; iy<int(gridy.size()); iy++){
         yx[k] = gridy[iy];
//std::cout << iy << " " << yx[k] << std::endl;
         double v = 2*b + yx.squaredNorm();
//if(iy==15) std::cout << iy << " " << v << std::endl;
         Eigen::MatrixXd xtx;
         if(npa > 0){
           Eigen::MatrixXd xpa(nsample,npa);
           for(int i=0; i<npa; i++) for(int z=0; z<nsample; z++)
             xpa(z,i) = xi(z,pa[i]);
//if(iy==15) std::cout << iy << " " << xpa(20,1) << std::endl;
           Eigen::MatrixXd I;
           I.setIdentity(npa,npa);
           xtx = I + hv*xpa.transpose()*xpa;
           xtx = xtx.ldlt().solve(I);
//if(iy==15) std::cout << iy << " " << xtx(2,1) << std::endl;
           Eigen::VectorXd xy = xpa.transpose()*yx;
           v -= hv*xy.dot(xtx*xy);
//if(iy==15) std::cout << iy << " " << v << std::endl;
//if(iy==3) exit(1);
         }
         double z = -(nsample/2 + a)*log(v);
//if(iy==15) std::cout << iy << " " << z << std::endl;
         if(npa > 0){
           Eigen::FullPivLU<Eigen::MatrixXd> lu(xtx);
           z += 0.5*log(lu.determinant());
         }
//std::cout << iy << " " << z << std::endl;
         double m = mu[idx];
         double sg = sigma[idx];
         Eigen::VectorXd e = sg*yx +
                   Eigen::VectorXd::Constant(yx.size(),m);
         double lkh = e.dot(ci.col(idx));
         for(int j=0; j<nsample; j++) lkh -= exp(e(j));
         z += lkh;
//std::cout << iy << " " << z << std::endl;
         prob.push_back(z);
         if(iy==0 || z > zmax) zmax = z;
       }
//std::cout << zmax << std::endl;
       double sum = 0;
       for(int iy=0; iy<int(prob.size()); iy++){
         double f = exp(prob[iy] - zmax);
         sum += f;
         prob[iy] = f;
       }
       for(int iy=0; iy<int(prob.size()); iy++) prob[iy] /= sum;
//std::cout << prob[50] << std::endl;
//std::cout << prob[61] << std::endl;
//std::cout << prob[72] << std::endl;
       double ys = sample(gridy, prob, r);
       xi2(k, idx) = ys;
     }
  }
  delete[] ksamp;
  delete[] sid;
  gsl_rng_free(r);

  return xi2;
}
