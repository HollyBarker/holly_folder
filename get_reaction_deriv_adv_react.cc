inline virtual void get_reaction_adv_diff_react(const unsigned& ipt,
                                                 const Vector<double> &C,
						 const DenseMatrix <double> &dCdx,
                                                 Vector<double>& R) const
  {
   //If no wind function has been set, return zero
   if(Reaction_fct_pt==0)
    {
     for(unsigned r=0;r<NREAGENT;r++) {R[r]= 0.0;}
    }
   else
    {
     // Get reaction terms
     (*Reaction_fct_pt)(C, dCdx ,R);
    }
  }

 /// \short Get the derivatives of the reaction terms with respect to the 
 /// concentration variables. If no explicit function pointer is set,
 /// these will be calculated by finite differences

 virtual void get_reaction_deriv_adv_diff_react(const unsigned& ipt,
						const Vector<double> &s,
                                                const Vector<double> &C,
						const DenseMatrix <double> &dCdx,
                                                DenseMatrix<double> &dRdC)
  const
  {
   //If no reaction pointer set, return zero
   if(Reaction_fct_pt==0)
    {
     for(unsigned r=0;r<NREAGENT;r++)
      {
       for(unsigned p=0;p<NREAGENT;p++)
        {
         dRdC(r,p) = 0.0;
        }
      }
    }
   else
    {
     //If no function pointer get finite differences
     if(Reaction_deriv_fct_pt==0)
      {
       //Local copy of the unknowns and their derivatives
       Vector<double> C_local = C;
       DenseMatrix<double> dCdx_local= dCdx;
       //Finite differences
       Vector<double> R(NREAGENT), R_plus_C(NREAGENT), R_minus_C(NREAGENT);
       Vector<double> R_plus_dCdx(NREAGENT), R_minus_dCdx(NREAGENT);
       //Get the initial reaction terms
       //(*Reaction_fct_pt)(C,R);
       const double fd_step = GeneralisedElement::Default_fd_jacobian_step;
       //Now loop over all the reagents
       //Holly - want to check the order of these loops is ok, and the +=
       
       for(unsigned p=0;p<NREAGENT;p++)
        {
         //Store the old value
         double old_var_C = C_local[p];
         //Increment the value
         C_local[p] += fd_step;
         //Get the new values
         (*Reaction_fct_pt)(C_local,dCdx_local,R_plus_C);
         //Reset the values
         C_local[p] = old_var_C;
         //Decrement the values
         C_local[p] -= fd_step;
	 
         //Get the new values
         (*Reaction_fct_pt)(C_local,dCdx_local, R_minus_C);
	 
         //Assemble the column of the jacobian
         for(unsigned r=0;r<NREAGENT;r++)
          {
           dRdC(r,p) = (R_plus_C[r] - R_minus_C[r])/(2.0*fd_step);
          }

         //Reset the value
         C_local[p] = old_var_C; 
	 }
       
         for (unsigned p=0;p<NREAGENT;p++)
	 {
	  for(unsigned i=0;i<DIM;i++)
	  {
	   double old_var_dCdx= dCdx_local(p,i);
	   //Increment the value
	   dCdx_local(p,i) += fd_step;
	   //Get the new value
	   (*Reaction_fct_pt)(C_local,dCdx_local,R_plus_dCdx);
	   //Reset the value
	   dCdx_local(p,i)=old_var_dCdx;

	   //Decrement the value
	   dCdx_local(p,i)-=fd_step;
	   //Get the new value
	   (*Reaction_fct_pt)(C_local,dCdx_local,R_minus_dCdx);
	   //Reset the value
	   dCdx_local(p,i)=old_var_dCdx;
	 
	   //Find out how many nodes there are in the element
	   const unsigned n_node = nnode();
	   //Set up memory for the shape and test functions
	   Shape psi(n_node);
	   DShape dpsidx(n_node,DIM);
	 
	   //Call the derivatives of the shape and test functions
	   dshape_eulerian(s,psi,dpsidx);
	   
	   for (unsigned k=0;k<n_node;k++)
	   {
	    for(unsigned r=0;r<NREAGENT;r++)
	    {
	     dRdC(r,p)+=(R_plus_dCdx[r] - R_minus_dCdx[r])*dpsidx(k,i)/(2.0*fd_step);
	    }
	   }
	  }
	 }
      }
     //Otherwise get the terms from the function
     else
      {
       //Holly- tells you when it's using the analytical dRdC
       //std::cout<<"not doing it"<<std::endl;
       (*Reaction_deriv_fct_pt)(C,dCdx,dRdC);
      }
    }
  }
