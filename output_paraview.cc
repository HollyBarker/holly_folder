void output_paraview(std::ofstream &file, 
		     const vector<double> & x0_vec,
		     const vector<double>u_val_vec) const
{

 // Change the scientific format so that E is used rather than e
 file.setf(std::ios_base::uppercase);

 // Decide how many elements there are to be plotted
 unsigned long number_of_nodes=x0_vec.size();
   
   
 // File Declaration
 //------------------
   
 // Insert the necessary lines plus header of file, and 
 // number of nodes and elements
 file_out 
  << "<?xml version=\"1.0\"?>\n"
  << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
  << "byte_order=\"LittleEndian\">\n"
  << "<UnstructuredGrid>\n" 
  << "<Piece NumberOfPoints=\""
  << number_of_nodes
  << "\" NumberOfCells=\""
  << number_of_nodes-1
  <<"\">\n";
   
 // Point Data
 //-----------

 // Check the number of degrees of freedom 
 unsigned ndof = fe_pt->nscalar_paraview();
   
 // Point data is going in here
 file << "<PointData ";
   
 // Insert just the first scalar name, since paraview reads everything
 // else after that as being of the same type. Get information from 
 // first element.
 file << "Scalars=\""
          << "V0"
          << "\">\n";
   
 // Loop over i scalar fields and j number of elements
 for(unsigned i=0;i<ndof;i++)
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
   file<<x0_vec[i]<<endl;
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
   file_out << i+counter << " "
              << i+1+counter 
              << std::endl;
  }
   
 file_out << "</DataArray>\n"
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
   
 file_out <<"</DataArray>\n"
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
