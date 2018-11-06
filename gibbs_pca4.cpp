#include <Rcpp.h>
using namespace Rcpp;

// a_kt is determined by eigenvectors, not drawn from posterior

// [[Rcpp::export]]

List gibbs_pca4(int M, int K, NumericMatrix centeredData, NumericVector adjMat, NumericVector a_0) {
  int D = centeredData.nrow();
  int N = centeredData.ncol();
  // set up initials
  NumericVector xi(M*K*D);
  NumericMatrix sigma2(M,N);
  NumericMatrix a(K,D);
  // Initialize loop variables.
  for (int k = 0; k<K;k++){
    a(k,0) = a_0(k);
    for (int t = 0; t<D;t++){
      xi(k*D+t) = R::rnorm(0,1);
    }
  }
      
  for (int t=1;t<D;t++){
    for (int k=0;k<K;k++){
      float tempa1 = 0;
      float tempa2 = 0;
      float tempa3 = 0;
      for (int j=0;j<N;j++){
        tempa1 = tempa1 + adjMat(K*D*j+k*D+t-1)*adjMat(K*D*j+k*D+t);
        tempa2 = tempa2 + adjMat(K*D*j+k*D+t);
        tempa3 = tempa3 + adjMat(K*D*j+k*D+t)*adjMat(K*D*j+k*D+t);
      }
      a(k,t) = (2*a(k,t-1)*tempa1 - tempa2)/(2*tempa3);
    }
  }
  
  for (int i = 1; i < M; i++){
    NumericMatrix tempxi(D,N);
    NumericMatrix temp(D,N);
    for (int j = 0; j < N; j++){
      float rateIG = 0;
      for (int t = 0; t < D; t++){
        for (int k = 0; k < K; k++){
          temp(t,j) = temp(t,j) + a(k,t)*(adjMat(K*D*j+k*D+t))*xi(K*D*(i-1)+k*D+t);
        }
        rateIG = rateIG + (centeredData(t,j) - temp(t,j))*(centeredData(t,j) - temp(t,j))/2;
      }
      sigma2(i,j) = 1/R::rgamma(D/2, 1/rateIG);
    }
    
    
    for (int t=0;t<D;t++){
      for (int k=0;k<K;k++){
        float xi_nu = 0;
        float xi_de = 0;
        for (int j=0;j<N;j++){
          float tempxi = 0;
          for (int h=0;h<K;h++){
            if (h<k){
              tempxi = tempxi + a(k,t)*(adjMat(K*D*j+k*D+t)*xi(K*D*(i)+h*D+t));
            } else if (h>k){
              tempxi = tempxi + a(k,t)*(adjMat(K*D*j+k*D+t)*xi(K*D*(i-1)+h*D+t));
            } else {}
          }
          xi_nu = xi_nu + (centeredData(t,j)-tempxi)*(adjMat(K*D*j+k*D+t))*a(k,t)/sigma2(i,j);
          xi_de = xi_de + 1+adjMat(K*D*j+k*D+t)*a(k,t)*adjMat(K*D*j+k*D+t)*a(k,t)/sigma2(i,j);
        }
        xi(K*D*(i)+k*D+t) = R::rnorm(xi_nu/xi_de, sqrt(1/xi_de));
      }
    }
  }
  
  
  return List::create(
    _["a"] = a,
    _["sigma2"] = sigma2,
    _["xi"] = xi
  );
}
