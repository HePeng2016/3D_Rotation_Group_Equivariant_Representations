# Introduction
# Installation
     include("./equivalentFeatures.jl")
     using  .equivalentFeatures 
# Initial
     equivalentFeatures.setN(N_new);
*N_new* is the maximum orbital quantum number *l* of the input spherical harmonic representation *Y(lm)*.   
    If this command is missing, then the default maximum orbital quantum number *l* is 3.  
     
     equivalentFeatures.Initial();
   Once this function had been carried out, the tables for the Clebsch Gordan coefficients and the Winger 3j coefficients were generated. And the data structures for storing the relations between coefficients and variables were generated.
# Usage
    equivalentFeatures.CtoS_Encode(V1);

  This function changes the cartesian coordination terms to the spherical harmonic tensor terms. 
  
  *V1* is the vector that records the multipole moment terms of the multipole expansion in the Cartesian system of coordinates.
  The order of the terms in this vector is determined by the length of the term, if the terms are equally long, they are ordered according to the lexicographical order. 
  
  e.g.
  
   *x<y<z*,
   *xx<x* ,
   *xx<xy<xz<yy<yz<zz*
   
*V1* is formatted as: 

N=1: [monopole]

*[1]*
    
N=2: [monopole,dipole]

*[1,x,y,z]*
    
N=3: [monopole,dipole,quadruple]

*[1,x,y,z,xx,xy,xz,yy,yz,zz]*

The output is a spherical harmonic tensor stored as a complex vector.

e.g.

N=2:

*[Y(l=0,m=0),Y(l=1,m=-1),Y(l=1,m=0),Y(l=1,m=1),Y(l=2,m=-2),Y(l=2,m=-1),Y(l=2,m=0),Y(l=2,m=1),Y(l=2,m=2)]*

*Y(l,m)* is the term for spherical multipole moment in the spherical coordinate system. *l* is the orbital quantum number and *m* is the azimuthal quantum number. 

     equivalentFeatures.SelfProduct(V1)

This function will change a spherical harmonic tensor to a rotation equivariant embedding vector. 
e.g. 

    V1 =[1.0,0.0043477849927746155,0.0,0.9999905483381614,0.11]; 
    V2 =[1.0,0.772027518982468,0.33454525822573616,0.5404192632877276];
    include("./equivalentFeatures.jl")
    using  .equivalentFeatures
    equivalentFeatures.setN(2); 
    equivalentFeatures.Initial(); 
    S1 = equivalentFeatures.CtoS_Encode(V1);
    S2 = equivalentFeatures.CtoS_Encode(V2);
    E1 = equivalentFeatures.SelfProduct(S1); 
    E2 = equivalentFeatures.SelfProduct(S2);
    Loss = sum(abs.(E1 - E2))   
    
V1 and V2 are two identical spherical harmonic tensors with different rotations. 

     equivalentFeatures.W3jProduct(V1,V2,V3) 

This function transforms three spherical harmonic tensors into a vector that is invariant to rotation. These three spherical harmonic tensors could be the tensors that represent the density of two nodes with different centers and the tensor that is derived from the subtraction between the coordinates of two nodes.

      equivalentFeatures.DecodeMatrix(V1,V2)
      
This function will return a matrix that can be used to convert the invariant coding calculated by the W3jProduct function into one of the original spherical harmonic tensors (V3).  And the V1, V2 are two other original spherical harmonic tensors.

e.g.
      
       W3 = equivalentFeatures.W3jProduct(V1,V2,V3)
       M  = equivalentFeatures.DecodeMatrix(V1,V2) 
       norm(M*W3-V3,2);

 V1, V2, V3 are three original spherical harmonic tensors.
 W3 is the invariant encoding derived from the W3j product of these three tensors. 
 M is the matrix that converts the invariant coding (W3) into the original spherical harmonic tensor (V3).


      equivalentFeatures.ProductEncode(V1,V2)


This function will return the invariant coding of the tensor product of two spherical harmonic tensors (V1, V2).

      V1 =[1.0,0.0043477849927746155,0.0,0.9999905483381614,0.11]; 
      V2 =[1.0,0.772027518982468,0.33454525822573616,0.5404192632877276];
      include("./equivalentFeatures.jl")
      using  .equivalentFeatures
      equivalentFeatures.setN(2); 
      equivalentFeatures.Initial();
      S1 = equivalentFeatures.CtoS_Encode(V1);
      S2 = equivalentFeatures.CtoS_Encode(V2);
      equivalentFeatures.ProductEncode(S1,S2);

     
