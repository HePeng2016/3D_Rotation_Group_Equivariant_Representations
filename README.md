# Introduction
# Installation
     include("./equivalentFeatures.jl")
     using  .equivalentFeatures 
# Initial
     equivalentFeatures.setN(N_new);
*N_new* is the orbital quantum number *l* of the input spherical harmonic representation *Y(lm)*.   
    If this command is missing, then the default orbital quantum number *l* is 3.  
     
     equivalentFeatures.Initial();
   Once this function had been carried out, the tables for the Clebsch Gordan coefficients and the Winger 3j coefficients were generated. And the data structures for storing the relations between coefficients and variables were generated.
# Usage
    equivalentFeatures.CtoS_Encode(V1);

  V1 is the vector that records the multipole moment terms of the multipole expansion in the Cartesian system of coordinates.
  The order of the terms in this vector is determined by the length of the term, if the terms are equally long, they are ordered according to the lexicographical order. 
  
  e.g.
  
   *x<y<z*
   
   *xx<xy<xz<yy<yz<zz*
   
V1 is formatted as: 

N=1 [monopole]

[1]
    
N=2 [monopole,dipole]

[1,x,y,z]
    
N=3 [monopole,dipole,quadruple]

[1,x,y,z,xx,xy,xz,yy,yz,zz]

   

