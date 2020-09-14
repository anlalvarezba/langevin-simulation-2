#define _USE_MATH_DEFINES

#include<random>
#include<cstdio>
#include<cmath>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<vector>

  double b  = 50.0;

int main()
{
  double ti = 0.0;
  double tf = 17.0;
  //int ti = 0;
  //int tf = 17;
  double v0 = 1.7331;
  double r0 = 0.0;

 

  //double b  = 50.0;
  double F0 = 0.0;
  double m  = 1.0;
  double k  = 1.0;
  double T  = 1.0;

  std::random_device rd;
  std::mt19937 gen(rd());

  int traj=100000; //numero de trayectorias para las que se realiza la simulacion


  double* table = new double [traj*(int(tf)*4)]; //eq 7
  double* average = new double [(int(tf)*4)]; //eq 7
  
  double* analytic = new double[(int(tf)*4)]; // analitica

  double* eq8 = new double [traj*(int(tf)*4)]; //eq 8
  double* promedio8 = new double [(int(tf)*4)]; //eq 8
  
  for(int c=0; c<traj; c++)
    {
  
      for(int j=ti; j<tf; j++)
    {
      double t=j/100.0;

      //////////////////////////////////////////////// Equation 7 ////////////////////////////////////////
      
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

      table[0+j*4+c*(int(tf)*4)]= tau7;
      table[1+j*4+c*(int(tf)*4)]= u7*u0;
      table[2+j*4+c*(int(tf)*4)]= u7*(x7-x0);
      table[3+j*4+c*(int(tf)*4)]= pow(x7-x0,2);

      /////////////////////////////////////// Equation 8 ///////////////////////////////////////////////////
      


      float C1 = 2*b*t - 3 + 4*exp(-b*t) - exp(-2*b*t);
      float sd_8r=sqrt((3*k*T/(m*pow(b,2)))*C1);
      float sd_8v=sqrt((((6*k*T)/(m))*(b*t*(1 - exp(-2*b*t))-2*pow(1-exp(-b*t),2)))/C1);
      
      std::normal_distribution<float> dis_8r{0,sd_8r}; //mean and standard deviation
      std::normal_distribution<float> dis_8v{0,sd_8v}; //mean and standard deviation
      
      double B3 = dis_8r(gen);
      double B4 = dis_8v(gen);


      double r8, v8;
      
      if(t>0){
	
      r8 = r0 + (v0/b)*(1-exp(-b*t)) + F0/(m*b)*(t-(1/b)*(1-exp(-b*t))) + B3 ;
      v8 = v0*(2*b*t*exp(-b*t) - 1 + exp(-2*b*t))/C1 + b*(r8 - r0)*pow(1 - exp(-b*t), 2)/C1 + (F0/(m*b))*(b*t*(1 - exp(-2*b*t)) - 2*pow(1-exp(-b*t),2)) + B4;
      }

      else{
	r8 = r0;
	v8 = v0;
      }

      
      double u8=sqrt(m/(3*k*T))*v8;
      //double u0=sqrt(m/(3*k*T))*v0;
      double x8=b*sqrt(m/(3*k*T))*r8;
      //double x0=b*sqrt(m/(3*k*T))*r0;

      eq8[0+j*4+c*(int(tf)*4)]= tau7;
      eq8[1+j*4+c*(int(tf)*4)]= u8*u0;
      eq8[2+j*4+c*(int(tf)*4)]= u8*(x8-x0);
      eq8[3+j*4+c*(int(tf)*4)]= pow(x8-x0,2);
     
      
      
      /*
      double Au = exp(-tau7); //expresion analitica para <u(tau)*u(0)>
      double Ax = 2*(tau7-1+exp(-tau7)); //expresion analitica para  <x(tau)-x(0)>
      double Aux = 1 - exp(-tau7); //expresion analitica para <u(tau)*(x(tau)-x(0))>
      */

      //printf("%5.1f  %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f \n", tau7, Au, u7*u0, Aux, u7*(x7-x0), Ax, pow(x7-x0,2));
    }
    }

  for(int i=0; i<4; i++)
    {
      for(int j=0; j<int(tf); j++)
      {
	double temporal = 0.0;
	for(int k=0; k<traj; k++)
	  {
	    temporal = temporal + table[i+j*4+k*(int(tf)*4)];
	  }
	average[i+j*4] = temporal/traj;
	//printf("%5.3f \n", average[i+j*4]);
      }
  }

    for(int i=0; i<4; i++)
    {
      for(int j=0; j<int(tf); j++)
      {
	double temporal = 0.0;
	for(int k=0; k<traj; k++)
	  {
	    temporal = temporal + eq8[i+j*4+k*(int(tf)*4)];
	  }
	promedio8[i+j*4] = temporal/traj;
	//printf("%5.3f \n", average[i+j*4]);
      }
  }


  for(int j=ti; j<tf; j++)
    {
      double t=j/100.0;
      double tau=b*t;
      analytic [j*4] = tau;
      analytic [1+j*4] = exp(-tau); //expresion analitica para <u(tau)*u(0)>
      analytic [2+j*4] = 1 - exp(-tau); //expresion analitica para <u(tau)*(x(tau)-x(0))>
      analytic [3+j*4] = 2*(tau-1+exp(-tau)); //expresion analitica para  <x(tau)-x(0)>

      
    }


  printf("%8s %25s %34s %26s \n", "\u03C4" ,"<u(\u03C4)\u2219u(0)>", "<u(\u03C4)\u2219[x(\u03C4)-x(0)]>", "<[x(\u03C4)-x(0)]\u00B2>");
  printf("\n");
  
  for(int i=0; i<int(tf); i++)
      {
	printf("%8.1f  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n", analytic[i*4], analytic[i*4+1], average[i*4+1], promedio8[i*4+1], analytic[i*4+2], average[i*4+2], promedio8[i*4+2], analytic[i*4+3], average[i*4+3], promedio8[i*4+3]);	
      }
  
  delete [] table;
  delete [] average;
  delete [] analytic;
  delete [] eq8;
  delete [] promedio8;
}
