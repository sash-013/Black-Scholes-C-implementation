// Black scholes model assumes no dividends (q=0)
#include <bits/stdc++.h>
using namespace std;
#define vd vector<double>
#define int double
#define nl '\n'
#define rep(i,a,b) for(int i=a;i<b;i++)
// N(0,1) density
//Returns normal probability density(PDF) of x
int f (int x){  
  int pi = 4.0*atan(1.0);
  return exp(-x*x*0.5)/sqrt(2*pi);
}

// Boole's Rule for Numerical Integration 
// calculating definite integral in a range
int Boole (int StartPt, int EndPt, int n){

  vd X(n + 1, 0.0);
  vd Y(n + 1, 0.0);

  int dx =(EndPt-StartPt)/int(n);

  rep(i,0,n+1){
      X[i]=StartPt+i*dx;
      Y[i]=f(X[i]);
  }
    
  int sum = 0;

  for(int t=0; t<=(n-1)/4;t++){
      int ind=4*t;
      sum +=(2/45.0)*(7*Y[ind]+ 32*Y[ind+1] + 12*Y[ind+2] + 32*Y[ind+3] + 7*Y[ind+4]) * dx;
  }

  return sum;
}
// N(0,1) cdf by Boole's Rule
// provides cdf of x
int Normal_dist (int x){
  return Boole (-10.0, x, 240);
}

// Black-Scholes Call Price

int BlackScholesPrice (int S, int K, int T, int r, int q, int v, char OpType) {
  
  int d = (log (S/K)+T*(r-q+0.5*v*v))/(v*sqrt(T));
  int call = S*exp(-q*T)*Normal_dist(d)-exp(-r*T)*K*Normal_dist(d-v*sqrt(T));

  if (OpType == 'C')
    return call;
  else
    // Put Parity
    return call - S*exp (-q * T) + K * exp (-r * T);
  

}

signed main (){
   //Sample Input 
   // S = 100.0;  // Stock Price
   // K = 100.0;  // Strike Price - (underlying price at which the option holder can buy or sell the underlying asset)
   // T = 1;   // Years to maturity - (time after which the option expires)
   // r = 0.05;  // Risk free interest rate - (rate of return on a risk free investment(eg - govt. bond))
   // q = 0.0;  // Dividend yield - rate of return earned on the underlying asset through dividends
   // v = 0.20;  // Yearly volatility - measure of stock price's variability. Higher volatility = Higher option price
   // OpType = 'C';  // 'C'all or 'P'ut
  int S,K,T,r,q,v;
  char OpType;  
  cout<<"Enter Stock price: ";
  cin>>S;
  cout<<"Enter Strike price: ";
  cin>>K;
  cout<<"Enter number of years to maturity: ";
  cin>>T;
  cout<<"Enter risk free interest rate: ";
  cin>>r;
  cout<<"Enter continuous dividend yield rate: ";
  cin>>q;
  cout<<"Enter yearly volatility: ";
  cin>>v;
  cout<<"Enter option type('C' for Call/'P' for Put): ";
  cin>>OpType;
  int ans = BlackScholesPrice (S, K, T, r, q, v, OpType) ;
  cout <<endl<< "Black Scholes"<<(OpType=='C'?" Call Option ":" Put Option ") <<"Price " << ans<< endl;
  return 0;
}

