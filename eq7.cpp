#define _USE_MATH_DEFINES

#include<random>
#include<cstdio>
#include<cmath>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<vector>

int main()
{
  double ti = 0.0;
  double tf = 16.0;
  double v0 = 1.7331;
  double r0 = 0.0;

  double b  = 50.0;
  double F0 = 1.0;
  double m  = 1.0;
  double k  = 1.0;
  double T  = 1.0;

  std::random_device rd;
  std::mt19937 gen(rd());
  
  
  for(int j=ti; j<=tf; j++)
    {
      double t=j/100.0;

      float sd_r=sqrt((3*k*T/m)*(1-exp(-2*b*t)));
      float sd_v=sqrt(((6*k*T)/(m*pow(b,2)))*(b*t-2*(1-exp(-b*t))/(1+exp(-b*t))));
      
      std::normal_distribution<float> dis_r{0,sd_r}; //mean and standard deviation
      std::normal_distribution<float> dis_v{0,sd_v}; //mean and standard deviation
      
      double B1 = dis_r(gen);
      double B2 = dis_v(gen);

      double v=v0*exp(-b*t)+(F0/(m*b))*(1-exp(-b*t)) + B1;
      double r=r0 + (1/b)*(v+v0-((2*F0)/(m*b)))*((1-exp(-b*t))/(1+exp(-b*t))) + (F0/(m*b))*t + B2 ;
      
      double R=r-r0-(v0/b)*(1-exp(-b*t))-(F0/(m*b))*(t-(1/b)*(1-exp(-b*t)));
      double V=v-(v0*exp(-b*t))-(F0/(m*b))*(1-exp(-b*t));
      double G=((k*T)/m)*(1-exp(-2*b*t));
      double H=((k*T)/(m*b))*pow((1-exp(-b*t)),2);
      double E=((k*T)/(m*pow(b,2)))*(2*b*t-3+4*exp(-b*t)-exp(-2*b*t));
      double w=pow(4*pow(M_PI,2)*(E*G-pow(H,2)),-3/2)*exp(-(G*pow(R,2)-2*H*R*V+E*pow(V,2))/(2*(E*G-pow(H,2))));


      //double B1 = sqrt(((3*k*T)/m)*(1-exp(-2*b*t)));
      //double B2 = sqrt(((6*k*T)/(m*pow(b,2)))*(b*t-2*(1-exp(-b*t))/(1+exp(-b*t)))) ;
      
      double tau7=b*t;
      double u7=sqrt(m/(3*k*T))*v;
      double u0=sqrt(m/(3*k*T))*v0;
      double x7=b*sqrt(m/(3*k*T))*r;
      double x0=b*sqrt(m/(3*k*T))*r0;

      double Au = exp(-tau7); //expresion analitica para <u(tau)*u(0)>
      double Ax = 2*(tau7-1+exp(-tau7)); //expresion analitica para  <x(tau)-x(0)>
      double Aux = 1 - exp(-tau7); //expresion analitica para <u(tau)*(x(tau)-x(0))>

      printf("%5.1f  %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f \n", tau7, Au, u7*u0, Aux, u7*(x7-x0), Ax, pow(x7-x0,2));

    }
  
}
