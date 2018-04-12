
#include "generic.h"
#include "meshes/one_d_mesh.h"
#include "flux_elements.h"
#include "bulk_elements.cc"
#include "bulk_elements.h"


//===========================
//DEMODRIVER FILE STARTS HERE
//===========================

//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented,
//LIC// multi-physics finite-element library, available
//LIC// at http://www.oomph-lib.org.
//LIC//
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1282 $
//LIC//
//LIC// $LastChangedDate: 2017-01-16 08:27:53 +0000 (Mon, 16 Jan 2017) $
//LIC//
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
//LIC//
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC//
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC//
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC//
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC//
//LIC//====================================================================

/// General form of the equations which are inputted here
/// \f[ 
/// \tau_{i} \frac{\partial C_{i}}{\partial t} 
/// + w_{j} \frac{\partial C_{i}}{\partial x_{j}}
/// + R_{i}(C_{j},\frac{\partial C_j}{\partial x_k})
/// + fct_{i}= 
/// D_{i}\frac{\partial^2 C_{i}}{\partial x_j^2}
/// +\frac{\partial}{\partial x_k}(F_i(C_j,\frac{\partial C_j}
/// {\partial x_l}))
/// \f]

 
//Standard C++ namespace
using namespace std;
//The oomph-lib namespace
using namespace oomph;
//Define Global variables in a namespace to keep things neat
namespace GlobalVariables
{
 //Equations for time dependent concentration, temperature and f
 //(local hydride volume fraction)
 //0 corresponds to the concentration equation
 //1 corresponds to the temperature equation
 //2 corresponds to the stress equation
 //3 corresponds to the volume fraction equation

 //Needed for the FD calculation of df/dt
 //This is set in main
 double timestep=0.0;
 
 //The timescale parameters (multiplying the time derivatives)
 //These are set in the function activator_inhibitor_tau
 Vector<double> Tau(5,1.0);

 //The diffusion parameters (multiplying the second spatial derivatives)
 //These are set in main
 Vector<double> D(5,1.0);

 //Proportionality constants for the flux BCs
 double lambda_C=-1.0;
 double lambda_T=-1.0;

 //Young's modulus (Fisher 1965) [Pa]
 double E=1.9e11;

 //Coefficient of thermal expansion
 double alpha=3.3e-5;
 
 //Temperature scaling - RT [K]
 double T0=298;
 //Concentration scaling - atmospheric concentration of hydrogen is set to
 //1.0 [mol/m^3]
 double C0=1.0;

 //Diffusivity (J.B. Condon - 1974) (assumed constant) [m^2 s^-1]
 double diffusivity=1.0;//(1.9e-6)*exp(-5820/T0);

 //Partial molal volume of hydrogen (Dutton 1977) [m^3 mol^-1]
 double Vh=7e-7;
 //Partial molal volume of hydride (Condon 1980) [m^3 mol^-1]
 double VHr=2.19e-5;
 //Partial molal volume of uranium (Condon 1980) [m^3 mol^-1]
 double VU=1.25e-5;

 //Bound concentration of hydrogen
 double Cb=1.0;
 //Terminal solid solubility limit (considered constant) [m^3 mol^-1] 
 double CTSS=0.56;
 //Density of uranium [kg m^-3]
 double rho=19.1e3;
 //Enthalpy of formation of hydride [Joules]
 double EoF=-30346;
 //Specific heat capacity at constant pressure [J kg^-1 K^-1]
 double c_p=117.2;
 //Conductivity of heat [W m^-1 K^-1]
 double kappa=22.5;
 //Heat of transport [J mol^-1]
 double Q=5400;
 //Ideal gas constant [J mol^-1 K^-1]
 double R=8.31;

 //Set the constants multiplying terms in energy conservation 
 double A1=(rho*c_p*T0*VHr)/EoF;
 double A2=(T0*VHr*kappa)/(diffusivity*EoF);
 double A3=(C0*R*T0*VHr)/EoF;
 double A4=C0*Q*VHr/EoF;
 double A5=(-2*C0*Q*VHr*Vh*E)/EoF;
 double A6=(-C0*Q*VHr*Vh*E)/EoF;
 double A7=(C0*VHr*Vh*Vh*E*E)/(R*T0*EoF);

 //Set the constants multiplying terms in stress-strain relation
 double B1=alpha*T0;
 double B2=VHr/VU;
 double B3=C0*Vh;

 //Characteristic lengthscale of the microstructure
 double char_length=3.0e-5;

 //Characteristic lengthscale of the whole domain (domain is length 0.1m)
 //This is used to nondimensionalise t=a^2 t' / D
 double a=0.1;

 //Characteristic time of reaction
 double char_time=char_length*char_length/diffusivity;
 
 //Function to set tau(C)
 void activator_inhibitor_tau(const Vector<double> &C,
			      const DenseMatrix <double> &dCdx,
			      Vector <double> &Tau)
 {
  Tau[0]=(1-C[4]);
  Tau[1]=A1;
  Tau[2]=0.0;
  Tau[3]=0.0;
  Tau[4]=1.0;
 }

 void activator_inhibitor_diff(const Vector<double> &C,
			      const DenseMatrix <double> &dCdx,
			      Vector <double> &Diff)
 {
  Diff[0]=1.0;
  Diff[1]=A2;
  Diff[2]=0.0;
  Diff[3]=0.0;
  Diff[4]=0.0;
 }
 
 //Function to set reaction term
 void activator_inhibitor_reaction(const Vector<double> &C,
				   const DenseMatrix <double> &dCdx,
				   Vector<double> &R)
 {
  //Reaction term in mass conservation and energy conservation equations
  //includes df/dt, so I've included the equation from the FD here.
  //a^2 /D comes from the non-dimensionalisation of time
  double dfdt_FD= ((1-C[4])*(C[0]-CTSS)*a*a)/((char_time+timestep)*(diffusivity)*(Cb-CTSS));
  
  //R[0]=dfdt_FD*(Cb-C[0]);
  R[0]=0.0;
  R[1]=dfdt_FD -(A3*C[1]*dCdx(0,0)*dCdx(0,0)/C[0]
		 +A4*dCdx(0,0)*dCdx(1,0)/C[1]
		 +A5*dCdx(0,0)*dCdx(2,0)
		 +A6*C[0]*dCdx(1,0)*dCdx(2,0)/(C[1]*C[1])
		 +A7*C[0]*dCdx(2,0)*dCdx(2,0)/C[1]);
  R[2]=dCdx(2,0);
  R[3]=0;//C[2]-(dCdx(3,0)-B1*C[1]-(B2*C[4]+(1-C[4])*B3*C[0]));
  R[4]=-dfdt_FD;

  //Hollyyyyy:This is for the viscous Burgers equation
  //R[0] =C[0]*dCdx(0,0);
  //R[0]=0.0;
 }

 //Function to set F term
 void activator_inhibitor_f(const Vector<double> &C,
			    const DenseMatrix <double> &dCdx,
			    Vector<double> &F)
 {
  F[0]=(Q*C[0]*dCdx(1,0))/(R*T0*C[1]*C[1])-(Vh*E*C[0]*dCdx(2,0))/(R*T0*C[1]);
  F[1]=0.0;
  F[2]=0.0;
  F[3]=0.0;
  F[4]=0.0;

  //Hollyyyyy:This is for the viscous Burgers equation
  //F[0]=0;
  //F[0]=-0.5*C[0]*C[0];
 }

 //Function to set R derivative term by hand
 /*
   void activator_inhibitor_reaction_derivative(const Vector<double> &C, 
   const DenseMatrix <double> &dCdx,
   DenseMatrix<double> &dRdC)
   {
   dRdC(0,0) = 1.0;
   }
 */

 //Function to set F derivative term by hand
 /*
   void activator_inhibitor_F_derivative(const Vector<double> &C,
   const DenseMatrix <double> &dCdx,
   DenseMatrix<double> &dFdC)
   {
   dFdC(0,0) = 1.0;
   }
 */

 void get_exact_u(const double& t,
		  const Vector<double>& x,
		  Vector<double>& u)
 {
  // Calculate the solution here...
 }
}

namespace FirstBoundaryConditions
{
 void prescribed_flux_on_outer_boundary(const Vector<double>& x,
					const Vector<double>& C,
					Vector<double>& flux)
 {
  //flux[0]=GlobalVariables::lambda_C*(1-C[0])/GlobalVariables::diffusivity;
  //flux[1]=GlobalVariables::lambda_T*(1-C[1])/GlobalVariables::kappa;
  flux[0]=0.0;
  flux[1]=0.0;
  flux[2]=0.0;
  flux[3]=0.0;
  flux[4]=0.0;
 }
 void prescribed_flux_on_inner_boundary(const Vector<double>& x,
					const Vector<double>& C,
					Vector<double>& flux)
 {
  //double N[2]={1.0,0.0};
  flux[0]=0.0;
  flux[1]=0.0;
  flux[2]=0.0;
  flux[3]=0.0;
  flux[4]=0.0;
 }
}
//======start_of_problem_class============================================
/// 1D AdvectionDiffusionReaction problem discretised with refineable
/// 1D QAdvectionDiffusionReaction elements.
/// The specific type of element is specified via the template parameter.
/// (The bit in the angle brackets)
//========================================================================
template<class ELEMENT>
class RefineableOneDAdvectionDiffusionReactionProblem : public Problem
{
public:
 /// Constructor. No arguments
 RefineableOneDAdvectionDiffusionReactionProblem();
 /// Destructor (empty)
 ~RefineableOneDAdvectionDiffusionReactionProblem() {}
 /// Set the initial condition
 void set_initial_condition();
 /// Perform nstep timesteps of size dt
 void timestep(const double &dt, const unsigned &nstep);
 /// \short Overloaded Problem's access function to the mesh.
 /// Recasts the pointer to the base Mesh object to the actual mesh type.
 /// This is required so that we can call specific RefineableMesh functions
 /*Refineable*/OneDMesh<ELEMENT>* mesh_pt()
 {
  return dynamic_cast</*Refineable*/OneDMesh<ELEMENT>*>(Problem::mesh_pt());
 }
private:
 ///Internal storage for the timestep
 double Dt;


 /// Hollyyyyy: added the following functions to the problem class
 /// \short Create Advection Diffusion flux elements on boundary b of 
 /// the Mesh pointed to by bulk_mesh_pt and add them to the Mesh 
 /// object pointed to by surface_mesh_pt
 void create_flux_elements(const unsigned &b, Mesh* const &bulk_mesh_pt,
                           Mesh* const &surface_mesh_pt);

 /// \short Delete Advection Diffusion flux elements and wipe the surface mesh
 void delete_flux_elements(Mesh* const &surface_mesh_pt);

 /// Pointer to the "bulk" mesh
 OneDMesh<ELEMENT>* Bulk_mesh_pt;

 /// Pointer to the "surface" mesh
 Mesh* Outer_surface_mesh_pt;
 Mesh* Inner_surface_mesh_pt;
 
}; // End of problem class

//=====start_of_constructor===============================================
/// Constructor for AdvectionDiffusionReaction problem:
//========================================================================
template<class ELEMENT>
RefineableOneDAdvectionDiffusionReactionProblem<ELEMENT>::
RefineableOneDAdvectionDiffusionReactionProblem()
{

 //Allocate the timestepper (second order implicit)
 add_time_stepper_pt(new BDF<2>);
 
 // Set up the mesh
 // Number of elements initially
 const unsigned n = 100;
 
 // Domain length
// const double pi=acos(-1);
 const double length = 0.1;
 
 // Build and assign the refineable mesh, need to pass in number of
 // elements, length and the timestepper
 Bulk_mesh_pt =
  new /*Refineable*/OneDMesh<ELEMENT>(n,length,Problem::time_stepper_pt());

  Bulk_mesh_pt->boundary_node_pt(0,0)
   ->make_periodic(Bulk_mesh_pt->boundary_node_pt(1,0));
 /*
 //Hollyyyyy - this was part of de-bugging
 unsigned nnode=mesh_pt()->nnode();
 
 //Loop over and set all previous time values to zero 
 for (unsigned i=0;i<nnode;i++)
 {
 time_stepper_pt()->assign_initial_positions_impulsive(
 mesh_pt()->node_pt(i));
 }*/
 
 //Output the initial mesh
 unsigned nplot=5;
 ofstream filename("initial_mesh.dat");
 this->Bulk_mesh_pt->output(filename,nplot);

  
 //----------------------------------------------
 // Set the boundary conditions for this problem.
 //----------------------------------------------
 
 // The ID of the leftmost boundary (OUTER)
 unsigned left_boundary_id=0;
 
 // The ID of the rightmost boundary (INNER)
 unsigned right_boundary_id=1;
 
 //Hollyyyyy: below will need to be changed when going into
 //2D. Creates flux elements on the inner and outer boundaries.

 // Create "surface mesh" that will contain only the prescribed-flux 
 // elements. The constructor just creates the mesh without
 // giving it any elements, nodes, etc.
 Outer_surface_mesh_pt = new Mesh;

 // Create prescribed-flux elements from all elements that are 
 // adjacent to boundary 0, but add them to a separate mesh.
 // Note that this is exactly the same function as used in the 
 // single mesh version of the problem, we merely pass different Mesh pointers.
 create_flux_elements(left_boundary_id,Bulk_mesh_pt,Outer_surface_mesh_pt);
 
 // Repeat for the 'inner' boundary
 Inner_surface_mesh_pt= new Mesh;
 create_flux_elements(right_boundary_id,Bulk_mesh_pt,Inner_surface_mesh_pt);
 
 // Add the sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Outer_surface_mesh_pt);
 add_sub_mesh(Inner_surface_mesh_pt);

 // Combine all submeshes into a single Mesh
 build_global_mesh();

// Create/set error estimator (default)
 //mesh_pt()->spatial_error_estimator_pt() = new Z2ErrorEstimator;

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here.
 //Get the number of boundaries
 unsigned n_bound = Bulk_mesh_pt->nboundary();
 unsigned n_node = 0;
 //Loop over number of boundaries
 for(unsigned b=0;b<n_bound;b++)
 {
  //Set pinned outer boundary Dirichlet conditions
  if (b==0)
  {
   //Get the number of nodes on the boundary
   n_node = Bulk_mesh_pt->nboundary_node(b);
   //Loop over the boundary nodes
   for (unsigned n=0;n<n_node;n++)
   {
    
    //Set scaled temperature on outer boundary to 1 (RT) for all time
    Bulk_mesh_pt->boundary_node_pt(b,n)->set_value(1,1.0);
    
    //Set stress on outer boundary to 0 for all time
    Bulk_mesh_pt->boundary_node_pt(b,n)->set_value(2,0.0);

    //Set displacement on outer boundary to 0 for all time
    Bulk_mesh_pt->boundary_node_pt(b,n)->set_value(3,0.0);

    //Set volume fraction on outer boundary to 0 for all time
    Bulk_mesh_pt->boundary_node_pt(b,n)->set_value(4,0.0);
    
    //Pin temperature on outer boundary
    Bulk_mesh_pt->boundary_node_pt(b,n)->pin(1);

    //Pin stress on outer boundary
    Bulk_mesh_pt->boundary_node_pt(b,n)->pin(2);

    //Pin displacement on outer boundary
    Bulk_mesh_pt->boundary_node_pt(b,n)->pin(3);
    
    //Pin volume fraction on outer boundary
    Bulk_mesh_pt->boundary_node_pt(b,n)->pin(4);
   }
  }
  //Set pinned inner boundary Dirichlet conditions
  if (b==1)
  {
   //Get the number of nodes on the boundary
   n_node = Bulk_mesh_pt->nboundary_node(b);
   //Loop over the boundary nodes
   for (unsigned n=0;n<n_node;n++)
   {
    
    //Set scaled temperature on inner boundary to 1 (RT) for all time
    Bulk_mesh_pt->boundary_node_pt(b,n)->set_value(1,1.0);
    
    //Set stress on inner boundary to 0 for all time
    Bulk_mesh_pt->boundary_node_pt(b,n)->set_value(2,0.0);

    //Set displacement on inner boundary to 0 for all time
    Bulk_mesh_pt->boundary_node_pt(b,n)->set_value(3,0.0);
    
    //Set volume fraction on outer boundary to 0 for all time
    Bulk_mesh_pt->boundary_node_pt(b,n)->set_value(4,0.0);
    
    //Pin temperature on inner boundary
    Bulk_mesh_pt->boundary_node_pt(b,n)->pin(1);

    //Pin stress on inner boundary
    Bulk_mesh_pt->boundary_node_pt(b,n)->pin(2);

    //Pin displacement on inner boundary
    Bulk_mesh_pt->boundary_node_pt(b,n)->pin(3);
    
    //Pin volume fraction on outer boundary
    Bulk_mesh_pt->boundary_node_pt(b,n)->pin(4);

   }
  }
 }

 //----------------------------------------------
 //----------------------------------------------
 //----------------------------------------------
 //----------------------------------------------

//Hollyyyyy: viscous Burgers equation
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// CHANGE THE PREVIOUS REFERENCES
// WITHIN THE PROBLEM CONSTRUCTOR FROM
// BULK_MESH_PT BACK TO MESH_PT().
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 //
 // Set zero Dirichlet conditions for the boundaries...
/*
// The ID of the leftmost boundary
unsigned left_boundary_id=0;
 
// The ID of the rightmost boundary
unsigned right_boundary_id=1;

// Which node do we want on each boundary (1D mesh means there's
// only one node on each boundary)
unsigned i_node=0;
mesh_pt()->boundary_node_pt(left_boundary_id,i_node)->pin(0);
mesh_pt()->boundary_node_pt(right_boundary_id,i_node)->pin(0);
*/
 // ...or make the domain periodic by setting the values at the left-hand boundary
 // equal to those on the right
 //mesh_pt()->boundary_node_pt(left_boundary_id,i_node)
 //->make_periodic(mesh_pt()->boundary_node_pt(right_boundary_id,i_node));
 
 //----------------------------------------------
 //----------------------------------------------



 //Hollyyyyy
 //Trying to pin different things to see which cause the non-converging residuals
 n_node=Bulk_mesh_pt->nnode();
 for (unsigned n=0;n<n_node;n++)
 {
  //pinning the concentration to a 0.1 value at all nodes
  //Bulk_mesh_pt->node_pt(n)->set_value(0,0.1);
  //Bulk_mesh_pt->node_pt(n)->pin(0);
  //pinning the temperature to a RT at all nodes
  //Bulk_mesh_pt->node_pt(n)->set_value(1,1.0);
  //Bulk_mesh_pt->node_pt(n)->pin(1);
  //pinning the stress to a zero value at all nodes
  Bulk_mesh_pt->node_pt(n)->set_value(2,0.0);
  Bulk_mesh_pt->node_pt(n)->pin(2);
  //pinning the displacement to a zero value at all nodes
  Bulk_mesh_pt->node_pt(n)->set_value(3,0.0);
  Bulk_mesh_pt->node_pt(n)->pin(3);
  //pinning the volume fraction to a zero value at all nodes
  Bulk_mesh_pt->node_pt(n)->set_value(4,0.0);
  Bulk_mesh_pt->node_pt(n)->pin(4);
 }
 
 
 // Loop over the elements to set up element-specific things that cannot
 // be handled by the (argument-free!) ELEMENT constructor: Pass pointer
 // to source function
 //Get the number of elements in the bulk
 const unsigned n_element = Bulk_mesh_pt->nelement();
 //Loop over the bulk elements
 for(unsigned i=0;i<n_element;i++)
 {
  // Upcast from GeneralisedElement to the present element
  ELEMENT *elem_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(i));

  //Set the diffusion coefficients
  elem_pt->diff_fct_pt() = &GlobalVariables::activator_inhibitor_diff;
  //Set the reaction terms
  elem_pt->reaction_fct_pt() = &GlobalVariables::activator_inhibitor_reaction;
  //And their derivatives
  elem_pt->reaction_deriv_fct_pt() = 0;
  //Set the F terms
  elem_pt->f_fct_pt() = &GlobalVariables::activator_inhibitor_f;
  //And their derivatives
  elem_pt->f_deriv_fct_pt() = 0;
  //Set the timescales
  elem_pt->tau_fct_pt() = &GlobalVariables::activator_inhibitor_tau;

 }
//---------------------------------------------------------------------
 // Loop over the flux elements to pass pointer to prescribed flux function
 unsigned n_surface_element=Outer_surface_mesh_pt->nelement();
 for(unsigned e=0;e<n_surface_element;e++)
 {
  // Upcast from GeneralisedElement to AdvectionDiffusionReaction flux element
  AdvectionDiffusionReactionFluxElement<ELEMENT> *el_pt = 
   dynamic_cast< AdvectionDiffusionReactionFluxElement<ELEMENT>*>(
    Outer_surface_mesh_pt->element_pt(e));

  // Set the pointer to the prescribed flux function
  el_pt->flux_fct_pt() = 
   &FirstBoundaryConditions::prescribed_flux_on_outer_boundary;
  
 }


 // Loop over the flux elements to pass pointer to prescribed flux function
 n_surface_element=Inner_surface_mesh_pt->nelement();
 for(unsigned e=0;e<n_surface_element;e++)
 {
  // Upcast from GeneralisedElement to AdvectionDiffusion flux element
  AdvectionDiffusionReactionFluxElement<ELEMENT> *el_pt = 
   dynamic_cast< AdvectionDiffusionReactionFluxElement<ELEMENT>*>(
    Inner_surface_mesh_pt->element_pt(e));

  // Set the pointer to the prescribed flux function
  el_pt->flux_fct_pt() = 
   &FirstBoundaryConditions::prescribed_flux_on_inner_boundary;
 }
 
//--------------------------------------------------------------------- 
 // Set up equation numbering scheme
 cout << "Number of equations: " << assign_eqn_numbers() << std::endl;
 
} // End of constructor

//============start_of_create_flux_elements==============================
/// Create AdvectionDiffusion Flux Elements on the b-th boundary of 
/// the Mesh object pointed to by bulk_mesh_pt and add the elements 
/// to the Mesh object pointeed to by surface_mesh_pt.
//=======================================================================
template<class ELEMENT>
void RefineableOneDAdvectionDiffusionReactionProblem<ELEMENT>::
create_flux_elements(const unsigned &b, Mesh* const &bulk_mesh_pt,
                     Mesh* const &surface_mesh_pt)
{
 // Get the number ofbulk elements adjacent to boundary b
 unsigned n_element = bulk_mesh_pt->nboundary_element(b);

 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0;e<n_element;e++)
 {
  // Get pointer to the bulk element that is adjacent to boundary b
  ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
   bulk_mesh_pt->boundary_element_pt(b,e));
   
  // Find the index of the face of element e along boundary b
  int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);

  // Build the corresponding prescribed-flux element
  AdvectionDiffusionReactionFluxElement<ELEMENT>* flux_element_pt = new 
   AdvectionDiffusionReactionFluxElement<ELEMENT>(bulk_elem_pt,face_index);

  //Add the flux element to the surface mesh
  surface_mesh_pt->add_element_pt(flux_element_pt);

 } //end of loop over bulk elements adjacent to boundary b
} // end of create_flux_elements

//============start_of_delete_flux_elements==============================
/// Delete Advection Diffusion Flux Elements and wipe the surface mesh
//=======================================================================
template<class ELEMENT>
void RefineableOneDAdvectionDiffusionReactionProblem<ELEMENT>::
delete_flux_elements(Mesh* const &surface_mesh_pt)
{
 // How many surface elements are in the surface mesh
 unsigned n_element = surface_mesh_pt->nelement();

 // Loop over the surface elements
 for(unsigned e=0;e<n_element;e++)
 {
  // Kill surface element
  delete surface_mesh_pt->element_pt(e);
 }

 // Wipe the mesh
 surface_mesh_pt->flush_element_and_node_storage();

} // end of delete_flux_elements


//=====================================================================
/// Set the initial conditions 
//=====================================================================
template<class ELEMENT>
void RefineableOneDAdvectionDiffusionReactionProblem<ELEMENT>::
set_initial_condition()
{
 // Backup time in global Time object
 double backed_up_time=time_pt()->time();
 //Set the initial concentrations of the reagent
 //Get the number of nodes in the bulk
 unsigned n_node = Bulk_mesh_pt->nnode();
 
 // Set continuous times at previous timesteps:
 // How many previous timesteps does the timestepper use?
 int nprev_steps=time_stepper_pt()->nprev_values();
 Vector<double> prev_time(nprev_steps+1);
 for (int t=nprev_steps;t>=0;t--)
 {
  prev_time[t]=time_pt()->time(unsigned(t));
 }

 // Loop over current & previous timesteps
 for (int t=nprev_steps;t>=0;t--)
 {
  // Continuous time
  double time=prev_time[t];
  cout << "setting IC at time =" << time << std::endl;
  //Loop over the nodes in the bulk
  for(unsigned n=0;n<n_node;n++)
  {
   //Local pointer to the node in the bulk
   Node* nod_pt = Bulk_mesh_pt->node_pt(n);   

   const double pi=3.141592654;
   double x=nod_pt->x(0);
   double ICconc=sin(2*pi*x/0.1);

   //Intial small concentration value everywhere
   double initial_concentration=ICconc;
   //Initial temperature:scaled such that RT=1.0
   double initial_temperature=1.0;
  
   //Set the initial concentration of hydrogen
   nod_pt->set_value(t,0,initial_concentration);
   //Set the initial temperature everywhere
   nod_pt->set_value(t,1,initial_temperature);
   //Set the initial displacement to 0
   nod_pt->set_value(t,3,0.0);
   //Set the initial hydride volume fraction to zero
   nod_pt->set_value(t,4,0.0);
  }//Finish loop over bulk nodes

 }//Finish loop over previous timesteps
 
 // Reset backed up time for global timestepper
 time_pt()->time()=backed_up_time;
 //Document the initial solution
 ofstream filename("RESLT/initial.dat");
 //Plot the solution with 5 points per element
 Bulk_mesh_pt->output(filename,5);
 filename.close();
 //Set the initial values impulsive
 //i.e. assume that the solution has been at the initial condition for all
 //previous times
 //assign_initial_values_impulsive(Dt);
 initialise_dt(Dt);
}
//====================================================================
/// Timestep the problem for nstep timesteps of length dt
//===================================================================
template<class ELEMENT>
void RefineableOneDAdvectionDiffusionReactionProblem<ELEMENT>::timestep(
 const double &dt, const unsigned &nstep)
{
 //Set the problem's Dt to use in the inital condition
 Dt = dt;
 //Maximum adaptation for the first timestep
 unsigned max_adapt = 0;
 //Take the first timestep
 bool first = true;
 //Set the initial condition
 set_initial_condition();

 //Solve the first step (one solve is done before adapting the mesh)
 //Spatial adaptivity does not work here
 unsteady_newton_solve(Dt);
 //newton_solve();

 //Output the result
 {
  unsigned i=0;
  char file1[100];
  sprintf(file1,"RESLT/step%i.dat",i+1);
  ofstream out1(file1);
  Bulk_mesh_pt->output(out1,5);
  out1.close();
 }
 //Zero rounds of adaptation per timestep
 max_adapt = 0;
 //This is not the first timestep, so we shouldn't use the initial conditions
 first = false;
 //Loop over timesteps
 for(unsigned i=1;i<nstep;i++)
 {
  std::cout<<"timestep number  "<<i<<std::endl;
  //Take a timestep
  unsteady_newton_solve(dt,max_adapt,first);
  //Output the result
  char file1[100];
  sprintf(file1,"RESLT/step%i.dat",i+1);
  ofstream out1(file1);
  Bulk_mesh_pt->output(out1,5);
  out1.close();
 }
}
//======start_of_main=====================================================
/// Driver code for 1D AdvectionDiffusionReaction problem
//======================================================================== 
/// General form of the equations which are inputted here
/// Hollyyyyy : I know the indexing is wrong.
/// \f[ 
/// \tau_{i} \frac{\partial C_{i}}{\partial t} 
/// + w_{j} \frac{\partial C_{i}}{\partial x_{j}} = 
/// D_{i}\frac{\partial^2 C_{i}}{\partial x_j^2}
/// +\frac{\partial}{\partial x_k}(F_i(C_j,\frac{\partial C_j}
/// {\partial x_l}))
/// - R_{i}(C_{j},\frac{\partial C_j}{\partial x_k}) - fct_{i}
/// \f]

int main()
{
 //Set the timestep
 double dt = 0.0001;
 GlobalVariables::timestep=dt;

 //Set the number of timesteps to be taken
 unsigned nstep=10;
 
 //Set up the problem
 //------------------
 // DREIGIAU: There's an inherent problem with the 1D elements. Doesn't
 // seem to be a problem with 2D elements but we can't even initialise
 // a new object of the type used to template the problem below.
 // Create the problem with 1D three-node refineable elements from the
 // RefineableLineAdvectionDiffusionReactionElement family.
 
 RefineableOneDAdvectionDiffusionReactionProblem<
  /*Refineable*/QAdvectionDiffusionReactionElement<5,1,3> > problem;
 //Take four levels of uniform refinement to start things off
 //for(unsigned i=0;i<4;i++) { problem.refine_uniformly(); }
 //Now timestep the problem
 problem.timestep(dt,nstep);

} // End of main
