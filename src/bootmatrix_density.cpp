#include<Rcpp.h>
#include <stdlib.h>
#include <time.h>
#include<iostream>
using namespace Rcpp;
using namespace std;
//[[Rcpp::export]]
NumericVector get_density_function(IntegerMatrix m1,int N){
  int a;
  int b;
  int num=m1.nrow();
  NumericVector boot_stat(N);
  IntegerMatrix x1(num,num);
  for(int k=0;k<num;k++){
    for(int q=0;q<num;q++){
      x1(k,q)=0;
    }
  }
  IntegerVector blist(num);
  for(int t=0;t<num;t++){
    blist(t)=round(R::runif(-0.49,num-0.5));}
  srand((unsigned)time(NULL));
  for(int r=0;r<N;r++){
    double density=0;
    double value=0;
    for(int i=0;i<num;i++){
      a=blist[i];
      for(int j=0;j<num;j++){
        b=blist[j];
        if(a!=b){
          x1(i,j)=m1(a,b);
          x1(j,i)=m1(a,b);
        }
        else{
          a=round(R::runif(-0.49,num-0.5));
          b=round(R::runif(-0.49,num-0.5));
          while(a==b){
            b=round(R::runif(-0.49,num-0.5));
          }
          x1(i,j)=m1(a,b);
          x1(j,i)=m1(a,b);
        }
      }
    }
    for(int w=0;w<num;w++){
      for(int z=0;z<num;z++){
        if(x1(w,z)!=0){
          value=value+1;
        }else{
          value=value+0;
        }
      }
    }
    density=value/(num*(num-1));
    boot_stat(r)=density;
  }
  return (boot_stat);
}

