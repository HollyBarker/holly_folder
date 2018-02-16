#include<vector>
#include<iostream>
#include<cmath>
#include<fstream>
#include<algorithm>
#include<iomanip>
#include<string>
#include<sstream>
#include<cstdio>
#include<ostream>

using namespace std;
//-------------------------------------------------------------
//-------------------------------------------------------------
//-------------------------------------------------------------
//Function prototypes
//-------------------------------------------------------------
//-------------------------------------------------------------
//-------------------------------------------------------------
double IC(double x2);
double big_integrand(double D,double x0,double x1,double t,
		     double intd_IC);
double trap_integrate_IC(double x1);
double trap_integrate_full(double D,double x0, double t);
double phi(double D, double x0, double t);
double diff_phi(double D,double x0,double t);
void output_paraview(std::ofstream &file,
		     const vector<double> & x0_vec,
		     const vector<double>u_val_vec);


//-------------------------------------------------------------
//-------------------------------------------------------------
//-------------------------------------------------------------
//Main()
//-------------------------------------------------------------
//-------------------------------------------------------------
//-------------------------------------------------------------
int main()
{
 // Set the parameter D and the max numbers of the indep variables
 //to later be looped through
 double D=1.0, dt=0.01, Ntime=10, Nx0=100;

 //Initialise the vectors to hold the coordinates and solution
 vector<double> x0_vec(Nx0), u_val_vec(Nx0);
 
 //Set the problem domain
 double x0start=0.0, x0end=1.0;
 
 //Get the dx0
 double dx0=(x0end-x0start)/Nx0;
 
 //Instantiate the indep variables
 double t, x0;

 //Instantiate the file to be written to
 ofstream file;
 char fileName[100];

 
 //Perform the second integration to get the whole integrand
 //double int_big=trap_integrate_full(D,x0,t);

 //Loop over the timesteps - in the phi functionwe divide by t
 //so we can't get the initial condition this way
 for (int ntime=1;ntime<Ntime;ntime++)
 {
  //Give the time point
  t=ntime*dt;

  //Open a file for the solution at this timestep
  sprintf(fileName,"exact_sol_step%i.dat",ntime);
  file.open(fileName);
  if(!file) return 1;
  
  //Loop over the x nodes
  for (double nx0=0;nx0<Nx0;nx0++)
  {
   //Give the x point and put into the vector
   x0=nx0*dx0;
   x0_vec[nx0]=x0;
   
   //Get the phi and dphi/dx values
   double phi_val=phi(D,x0,t);
   double diff_phi_val=diff_phi(D,x0,t);
   
   //Get u from the substitution definition and put into the vector
   double u_val=-2*D*diff_phi_val/phi_val;
   u_val_vec[nx0]=u_val;
   
   //See the values in the command line
   //std::cout<<phi_val<<"      "<<diff_phi_val<<"       "<<u_val<<std::endl;
   
   file<<x0<<"   "<<u_val<<endl;
   
  }
  //output_paraview(file,x0_vec,u_val_vec);
  file.close();
 }
 return 0;
}

//-------------------------------------------------------------
//-------------------------------------------------------------
//-------------------------------------------------------------
//Function implementations
//-------------------------------------------------------------
//-------------------------------------------------------------
//-------------------------------------------------------------


//-------------------------------------------------------------
//Function to return the initial condition
//-------------------------------------------------------------
double IC(double x2)
{
 double spot=2*exp(-1000.0*(x2-0.5)*(x2-0.5));
 return spot;
}

//-------------------------------------------------------------
//Function to return the full integrand
//-------------------------------------------------------------
double big_integrand(double D, double x0,
		     double x1,double t,
		     double intd_IC)
 
{return exp(-pow((x0-x1),2)/(4*D*t)-(1/(2*D))*intd_IC);}

//-------------------------------------------------------------
//Function to integrate the IC
//-------------------------------------------------------------
double trap_integrate_IC(double x1)
{
 double xstart=0.0, numel=100;
 // Set dx, instantiate the total integral and the a and b points
 double total=0.0, a, b, dx=std::abs(xstart-x1)/numel;
 
 for (unsigned i=0;i<numel;i++)
 {
  //Get a and b
  a=xstart+i*dx;
  b=xstart+(i+1)*dx;
  //Increment the integral by this section 
  total+= 0.5*(IC(a)+IC(b))*dx;
 }
 return total;
}

//-------------------------------------------------------------
//Function to perform the second integration
//-------------------------------------------------------------
double trap_integrate_full(double D,
			   double x0,
			   double t)
{
 //Set the integration limits (approx -inf to inf)
 double xstart=-100.0, xend=100.0;
 
 //Set the number of trapeziums used to calculate the integration
 double numel=1000;
 
 //Get dx
 double dx=std::abs(xstart-xend)/numel;
 
 //Instantiate the a and b points and the integralof the IC upto
 //x1=a and x1=b;
 double a, b, intd_IC_a=0.0, intd_IC_b=0.0;
 
 //Instantiate the full integrand for x1=a and x1=b
 double full_integrand_a=0.0, full_integrand_b=0.0;

 //Instantiate the total integral to be incremented with each
 //trapezium
 double  total=0.0; 
 
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

//-------------------------------------------------------------
//Function to give phi
//-------------------------------------------------------------
double phi(double D, double x0, double t)
{
 //Get a value for pi
 double pi=acos(-1);
 //Perform the second integral calculation
 double int_big=trap_integrate_full(D,x0,t);
 //Get the value for phi
 double phi_val=int_big/(sqrt(4*pi*D*t));
 //Return phi
 return phi_val;
}

//-------------------------------------------------------------
//Function to differentiate phi by x
//-------------------------------------------------------------
double diff_phi(double D,double x0,double t)
{
 //Define small h=dx over which to calculate the derivative
 double h=0.001;
 //Get the two x values required to find the derivative
 double  x0p=x0+h/2, x0m=x0-h/2;
 //Central difference scheme for derivative
 double diff_phi_val=(1/h)*(phi(D,x0p,t)-phi(D,x0m,t));
 //Return dphi/dx
 return diff_phi_val;
}





//-------------------------------------------------------------
//Function to write to file : unused at the moment
//-------------------------------------------------------------

void output_paraview(std::ofstream &file, 
		     const vector<double> & x0_vec,
		     const vector<double>u_val_vec)
{

 // Change the scientific format so that E is used rather than e
 file.setf(std::ios_base::uppercase);

 // Decide how many elements there are to be plotted
 unsigned long number_of_nodes=x0_vec.size();
 unsigned long number_of_elements=number_of_nodes-1;
   
   
 // File Declaration
 //------------------
   
 // Insert the necessary lines plus header of file, and 
 // number of nodes and elements
 file
  << "<?xml version=\"1.0\"?>\n"
  << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
  << "byte_order=\"LittleEndian\">\n"
  << "<UnstructuredGrid>\n" 
  << "<Piece NumberOfPoints=\""
  << number_of_nodes
  << "\" NumberOfCells=\""
  << number_of_elements
  <<"\">\n";
   
 // Point Data
 //-----------

   
 // Point data is going in here
 file << "<PointData ";
   
 // Insert just the first scalar name, since paraview reads everything
 // else after that as being of the same type. Get information from 
 // first element.
 file << "Scalars=\""
          << "V0"
          << "\">\n";
   
 // Loop over i scalar fields and j number of elements
 for(unsigned i=0;i<1;i++)
  {
   file << "<DataArray type=\"Float32\" "
            << "Name=\""
            << "V0"
            << "\" "
            << "format=\"ascii\""
            << ">\n";

   for(unsigned j=0;j<number_of_nodes;j++)
    {
     file<<u_val_vec[j]<<endl;
    }
       
   // Close of the DataArray
   file << "</DataArray>\n";
  }
   
 // Close off the PointData set 
 file  << "</PointData>\n";
   
   
 // Geometric Points
 //------------------
   
 file
  << "<Points>\n"
  << "<DataArray type=\"Float32\""
  << " NumberOfComponents=\""
  // This always has to be 3 for an unstructured grid
  << 3  << "\" "
  << "format=\"ascii\">\n";
   
 // Loop over all the elements to print their plot points
 for(unsigned i=0;i<number_of_nodes;i++)
  {
   file<<x0_vec[i]<<" 0 0"<<endl;
  }
   
 file
  << "</DataArray>\n"
  << "</Points>\n";
   
   
 // Cells
 //-------
   
 file
  << "<Cells>\n"
  << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
   
 // Make counter for keeping track of all the local elements,
 // because Paraview requires global coordinates
 unsigned counter=0;
   
 // Write connectivity with the local elements
 for(unsigned i=0;i<number_of_elements;i++)
  {
   file << i+counter << " "
              << i+1+counter 
              << std::endl;
  }
   
 file << "</DataArray>\n"
          << "<DataArray type=\"Int32\" "
          << "Name=\"offsets\" format=\"ascii\">\n";
   
 // Make variable that holds the current offset number
 unsigned offset_sum=0;
   
 // Write the offset for the specific elements
 for(unsigned i=0;i<number_of_nodes;i++)
  {
     offset_sum+=2;
     file << offset_sum << std::endl;
  }
   
 file <<"</DataArray>\n"
          <<"<DataArray type=\"UInt8\" Name=\"types\">\n";
   
 // Loop over all elements to get the type that they have
 for(unsigned i=0;i<number_of_nodes;i++)
  {
   file << "3" << std::endl;
  }
   
 file <<"</DataArray>\n"
          <<"</Cells>\n";
   
   
 // File Closure
 //-------------
 file <<"</Piece>\n"
          <<"</UnstructuredGrid>\n"
          <<"</VTKFile>";
}

