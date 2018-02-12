#include<vector>
#include<iostream>
#include<cmath>
#include<fstream>
#include<algorithm>
#include<iomanip>

using namespace std;

//Sets the initial condition
double IC(double x2)
{
 double spot=2*exp(-1000*(x2-0.5)*(x2-0.5))-1;
 if (spot>0) {return spot;}
 else {return 0;}
}

double big_integrand(double D, double x0,
		     double x1,double t,
		     double intd_IC)
{return exp(-pow((x0-x1),2)/(4*D*t)-(1/(2*D))*intd_IC);}

double trap_integrate_IC(double x1)
{
 double xstart=0.0, numel=100;
 // Set dx, instantiate the total integral and the a and b points
 double total=0.0, a, b, dx=std::abs(xstart-x1)/numel;
 
 for (unsigned i=0;i<numel;i++)
 {
  a=xstart+i*dx;
  b=xstart+(i+1)*dx;
  total+= 0.5*(IC(a)+IC(b))*dx;
 }
 return total;
}
double trap_integrate_full(double D,
			   double x0,
			   double t)
{
 double xstart=-100.0, xend=100.0, numel=1000;
 // Set dx, instantiate the total integral and the a and b points
 double total=0.0, a, b, dx=std::abs(xstart-xend)/numel;
 double intd_IC_a=0.0, intd_IC_b=0.0;
 double full_integrand_a=0.0, full_integrand_b=0.0;
 //Loop over the trapezoids to integrate over
 for (unsigned i=0;i<numel;i++)
 {
  //Set the start and end point of the trapezoid
  a=xstart+i*dx;
  b=xstart+(i+1)*dx;
  //Get the IC integrated upto a...
  intd_IC_a=trap_integrate_IC(a);
  //... and then b
  intd_IC_b=trap_integrate_IC(b);
  //Get the integrand for x1=a...
  full_integrand_a=big_integrand(D,x0,a,t,intd_IC_a);
  //... and x1=b
  full_integrand_b=big_integrand(D,x0,b,t,intd_IC_b);
  //Add the trapezoid area to the integral total
  total+=0.5*(full_integrand_a+full_integrand_b)*dx;
 }
 //Return sum of all the trapezoid areas
 return total;
}

double phi(double D, double x0, double t)
{
  double pi=acos(-1);
  double int_big=trap_integrate_full(D,x0,t);
  double phi_val=int_big/(sqrt(4*pi*D*t));
  return phi_val;
}

double diff_phi(double D,double x0,double t)
{
 double h=0.001;
 double x0m=x0-h/2, x0p=x0+h/2;
 double diff_phi_val=(1/h)*(phi(D,x0p+h/2,t)-phi(D,x0m,t));
 return diff_phi_val;
}

int main()
{
 // Set the parameter D and the indep variables to later be looped through
 double D=1.0, dt=0.01, Ntime=100, Nx0=1000;
 double x0start=0.0, x0end=1.0;
 double dx0=(x0end-x0start)/Nx0;
 double t;
 double x0;


 //Perform the second integration to get the whole integrand
 double int_big=trap_integrate_full(D,x0,t);
 std::cout<<int_big<<std::endl;
 
 double phi_val=phi(D,x0,t);
 double diff_phi_val=diff_phi(D,x0,t);
 double u_val=-2*D*diff_phi_val/phi_val;

 for (double ntime=0;ntime<Ntime;ntime++)
 {
  t=ntime*dt;
  for (double nx0=0;nx0<Nx0;nx0++)
  {
   x0=nx0*
  }
 }
  
}
