# Introduction

This software package is used to represent the spherical harmonic tensors in a rotation equivariance or invariance way. This package is implemented in two different versions. One version uses Julia. The other version uses C++. They have the same functions.
# Installation (julia)
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

e.g.

      V1 =[1.0,0.0043477849927746155,0.0,0.9999905483381614,0.11]; 
      V2 =[1.0,0.772027518982468,0.33454525822573616,0.5404192632877276];
      include("./equivalentFeatures.jl")
      using  .equivalentFeatures
      equivalentFeatures.setN(2); 
      equivalentFeatures.Initial();
      S1 = equivalentFeatures.CtoS_Encode(V1);
      S2 = equivalentFeatures.CtoS_Encode(V2);
      equivalentFeatures.ProductEncode(S1,S2);

# Installation (c++)


   Eigen package is needed for install. 

       cd  include
       make
   
   rotation.o will be generated. 
   
   In your project file: 
   
       #include "include/rotation.h"
       equivalentFeatures Test; 
       Test.Initial();

   In your Makefile file: 

        objects = include/rotation.o 
        CC = g++  -std=c++11  -g  -I include 
        $(CC)  -o rotation yourproject.cpp $(objects)
  
# Usage   
    std::vector<std::complex<double>> RStoCS_Encode(const std::vector<double>& V1, int n); 
    std::vector<double> CStoRS_Encode(const std::vector<std::complex<double>>& V1, int n);

RStoCS_Encode will convert a real spherical harmonic tensor into a complex spherical harmonic tensor, CStoRS_Encode will convert a complex spherical harmonic tensor back to a real spherical harmonic tensor. 
    std::vector<std::complex<double>> IncreaseDegree(const std::vector<std::complex<double>>& V, int n);
IncreaseDegree function can increase the degree of spherical harmonic tensor by n. 
e.g.

    std::vector<double> V1 = { 1.0,0.2820947917738782,0.1558348923076441,-0.41829957442835725};
    std::vector< std::complex<double> > V2 = Test.RStoCS_Encode(V1,2);
    std::vector< std::complex<double> > V3 = Test.IncreaseDegree(V2,1);
    std::vector<double> V4 = Test.CStoRS_Encode(V3,3);
V1 is a real spherical harmonic tensor to be converted as a complex spherical harmonic V2 via RStoCS_Encode function, and the complex spherical harmonic V2 is increased by one degree using IncreaseDegree function, increased spherical harmonic tensor V3 is converted back as a  complex spherical harmonic tensor V4 using CStoRS_Encode. 

    std::vector<double> SelfProduct(const std::vector<std::complex<double>>& V, int n, int n2); 
 V is a spherical harmonic tensor to be encoded. 
 Two shells of V will be Clebsch–Gordan (CG) producted to generate a set of spherical harmonic tensors,the maximum degree of spherical harmonic tensors in this set is n, and the maximum size of the set is n2. The norm of each spherical harmonic tensor from this set derived from the shells (CG) production of V is calculated as return.  
 
    std::vector<double> SelfProductPairwise(const std::vector<std::complex<double>>& V, int n, int n2);
   
 V is a spherical harmonic tensor to be encoded. 
 Two shells of V will be  Clebsch–Gordan (CG)  producted to generate a set of spherical harmonic tensors, the maximum degree of spherical harmonic tensors in this set is n, and the maximum size of the set is n2. For this set derived from the shells (CG) production of V, the norm of each spherical harmonic tensor and the product of two spherical harmonic tensors are caculated as return. 

 
    std::vector<std::complex<double>> W3jProduct(const std::vector<std::complex<double>>& V1, const std::vector<std::complex<double>>& V2,const std::vector<std::complex<double>>& V3,int n1, int n2); 
   
V1, V2, V3 are spherical harmonic tensors with different orientation to be encoded with wigner 3J. Only shell n1 for V1 is selected, and for three shells wigner production from V1,V2,V3 with degrees s1,s2,s3 satisfied that s1 = n1, and s1 < abs(s2-s3)+n2. 

    std::vector<double> W3jProductCToR(const std::vector<std::complex<double>>& InvariantV, int n1,int d2, int d3, int n2);

The definition of n1 and n2 is the same of W3jProduct, d2,d3 are the dimensions of V2 and V3 that are as inputs of W3jProduct.This funcation will convert the spherical harmonic tensor that is the result of W3jProduct into a real vector. 

     std::vector<std::complex<double>> W3jProductRToC(const std::vector<double>& InvariantV, int n1,int d2, int d3, int n2);

The definition of n1,n2,d2,d3 is the same of W3jProductCToR, the input is a real vector converted by W3jProductCToR, the output is spherical harmonic tensor that is the same as input of W3jProductCToR. This function is the inverse process of the W3jProductCToR function.

      std::vector<std::complex<double>> W3jProductCompact(const std::vector<std::complex<double>>& V1,const std::vector<std::complex<double>>& V2,const std::vector<std::complex<double>>& V3,int n);

 V1, V2, V3 are spherical harmonic tensors with different orientation to be encoded with wigner 3J.  The length of result of wigner 3j product for these three spherical harmonic should be the same as V1, and this result is a rotation invariance vector. 
     
     Eigen::MatrixXcd DecodeMatrixCompact(const std::vector<std::complex<double>>& V2, const std::vector<std::complex<double>>& V3, int n);

V2, V3 are spherical harmonic tensors that are rotational equivariantly same as input of W3jProductCompact. This matrix will return a matrix M.  This matrix will recover an invariance vector yielded by W3jProductCompact into an spherical harmonic tensor that is rotational equivariantly same as the V1 that is the input of W3jProductCompact, the orientation of this spherical harmonic tensor is provied by V2, V3.  
e.g.

       equivalentFeatures Test;
       Test.Initial();

       std::vector<std::complex<double>> V1 = {
            {0.2820947917738782, 0.0}, {0.1558348923076441, -0.08714600779676653},
            {-0.41829957442835725, 0.0}, {-0.1558348923076441, -0.08714600779676653},
            {0.05400984229015042, -0.08789344893420385}, {-0.29831935733068915, 0.16682618799224708},
            {0.3780897816045623, 0.0}, {0.29831935733068915, 0.16682618799224708},
            {0.05400984229015042, 0.08789344893420385}
       };
      std::vector<std::complex<double>> V2 = {
        {0.2820947917738782, 0.0}, {-0.13150538452459107, 0.20714627161985197},
        {0.34398535815919695, 0.0}, {0.13150538452459107, 0.20714627161985197},
        {-0.08289419650884133, -0.17630513479494964}, {-0.20702023603219083, 0.32609668569054673},
        {0.1535728097664055, 0.0}, {0.20702023603219083, 0.32609668569054673},
        {-0.08289419650884133, 0.17630513479494964}
    };
     std::vector<std::complex<double>> V3 = {
        {0.2820947917738782, 0.0}, {-0.08323545825287791, -0.286716367230279},
        {0.2458934345224644, 0.0}, {0.08323545825287791, -0.286716367230279},
        {-0.24360329182248358, 0.15445613562862262}, {-0.09366656951197919, -0.32264781266424436},
        {-0.07575460360437422, 0.0}, {0.09366656951197919, -0.32264781266424436},
        {-0.24360329182248358, -0.15445613562862262}
    };
     std::vector<std::complex<double>> V1_ = {
        {0.2820947917738782, 0.0}, {0.1294158205116877, -0.22846248270670816},
        {0.3175614797663446, 0.0}, {-0.1294158205116877, -0.22846248270670816},
        {-0.11470689480471548, -0.19135797511692804}, {0.18808081324411152, -0.3320259406721129},
        {0.08429133903919787, 0.0}, {-0.18808081324411152, -0.3320259406721129},
        {-0.11470689480471548, 0.19135797511692804}
    };
    std::vector<std::complex<double>> V2_ = {
        {0.2820947917738782, 0.0}, {-0.11794378575971845, 0.12037014272293214},
        {-0.42653604644964793, 0.0}, {0.11794378575971845, 0.12037014272293214},
        {-0.0018711923316910401, -0.09188362988533172}, {0.23022904373214756, -0.23496534959002102},
        {0.4056684557059928, 0.0}, {-0.23022904373214756, -0.23496534959002102},
        {-0.0018711923316910401, 0.09188362988533172}
    };
     std::vector<std::complex<double>> V3_ = {
        {0.2820947917738782, 0.0}, {-0.04734136912624448, 0.3026803946881253},
        {0.22587421614524694, 0.0}, {0.04734136912624448, 0.3026803946881253},
        {-0.28921882314977415, -0.09274049507548826}, {-0.04893693653854026, 0.3128817678849446},
        {-0.11318592970238676, 0.0}, {0.04893693653854026, 0.3128817678849446},
        {-0.28921882314977415, 0.09274049507548826}
    };

    std::vector<std::complex<double>> V__ = Test.W3jProductCompact(V1,V2,V3,3); 
    Eigen::MatrixXcd M = Test.DecodeMatrixCompact(V2_,V3_,3); 
    std::vector<std::complex<double>> result(M.rows(), 0);

    for (size_t i = 0; i < M.rows(); ++i) {
        for (size_t j = 0; j < M.cols(); ++j) {
            result[i] += M(i,j) * V__[j];
          }
    }


    double norm = 0.0;
    for (size_t i = 0; i < result.size(); ++i) {
        norm += (std::pow(result[i] - V1_[i], 2)).real();
    }
    norm = std::sqrt(norm);
  V1, V2, V3  are same as  V1_, V2_, V3_  in rotation equivalention.  V1 is embedded as rotation invariant vector V__ using V2, V3.  V2_, V3_ generated a auxiliary matrix that can convert invariant vector V__ into a spherical harmonic tensor that is the same as V1_ and rotational equivalent to V1. 
  
       std::vector<double> W3jProductCompactCToR(const std::vector<std::complex<double>>& InvariantV, int n);
       std::vector<std::complex<double>> W3jProductCompactRToC(const std::vector<double>& InvariantV, int n);
    
W3jProductCompactRToC will change the result of W3jProductCompact into real vector, and W3jProductCompactRToC will recover the this real vector into original spherical harmonic tensor. 

       std::vector< std::vector<double> > ReferencesExtract(const  std::vector<double>& VR);
VR is a real spherical harmonic tensor, two reference vectors are extracted from dipole and quadrupole shells of this real spherical harmonic tensor using ReferencesExtract.  

       std::vector<std::vector<std::complex<double>>> SelfProductMatrix(const std::vector<std::complex<double>>& V, int n, int n2);
V is a complex spherical harmonic tensor, the pair of shells from this spherical harmonic tensor are producted as a set of spherical harmonic tensors, the maximum size of each set is n2, in each set, only spherical harmonic tensors with degree n are returned.   

The reference vectors for the W3jProductCompact and DecodeMatrixCompact functions can be generated by either the ReferencesExtract or the SelfProductMatrix functions.
e.g.

       std::vector<std::complex<double>> V_input = {
        {0.5641895835477564, 0.0},
        {0.02432950778305304, 0.12000026382308544},
        {-0.0743142162691603, 0.0},
        {-0.02432950778305304, 0.12000026382308544},
        {-0.028884354218690914, -0.2641985837291535},
        {-0.5053395933628799, 0.4929228736827938},
        {0.5316625913709678, 0.0},
        {0.5053395933628799, 0.4929228736827938},
        {-0.028884354218690914, 0.2641985837291535}
    };
    std::vector<std::complex<double>> V_output = {
        {0.5641895835477564, 0.0},
        {0.046180362258809796, -0.5151788499369871},
        {0.563454914288809, 0.0},
        {-0.046180362258809796, -0.5151788499369871},
        {-0.35831018662719905, -0.036901839488305416},
        {0.09441424373213234, -0.6546737533363572},
        {0.008536735434823653, 0.0},
        {-0.09441424373213234, -0.6546737533363572},
        {-0.35831018662719905, 0.036901839488305416}
    };
    /* reference vectors extraction via ReferencesExtract function*/
    std::vector<std::complex<double>>  V_input_ =  V_input;
    std::vector< std::vector<double> > V__ = Test.ReferencesExtract(Test.CStoRS_Encode(V_input_,3));
    std::vector<double> v1_r = {1};
    v1_r.insert(v1_r.end(), V__[0].begin(), V__[0].end()); 
    std::vector<double> v2_r = {1};
    v2_r.insert(v2_r.end(), V__[1].begin(), V__[1].end()); 
    std::vector<std::complex<double>>  v1 = Test.RStoCS_Encode(v1_r,2); 
    std::vector<std::complex<double>>  v2 = Test.RStoCS_Encode(v2_r,2); 
    v1 = Test.IncreaseDegree(v1,1);
    v2 = Test.IncreaseDegree(v2,1);
    std::vector<std::complex<double>> V_Encode = Test.W3jProductCompact(V_output,v1,v2,3);//Embedding as rotation-invariant coding 
    Eigen::MatrixXcd M = Test.DecodeMatrixCompact(v1,v2,3);//Generate an auxiliary matrix for recovering embedded. 
    std::vector<std::complex<double>> V_Re(M.rows(), 0);
    for (size_t i = 0; i < M.rows(); ++i) {
        for (size_t j = 0; j < M.cols(); ++j) {
            V_Re[i] += M(i,j) * V_Encode[j];
          }
    }// V_Re=M*V_Encode;
    double norm = 0.0;
    for (size_t i = 0; i < V_Re.size(); ++i) {
        norm += (std::pow(V_Re[i] - V_output[i], 2)).real();
    }
    norm = std::abs(norm);
    /* reference vectors extraction via SelfProductMatrix function*/
    std::vector<std::vector<std::complex<double>>>  V_ = Test.SelfProductMatrix(V_input,2,2);  
    v1.resize(1); 
    v2.resize(1);
    v1[0] = {1.0, 0.0}; 
    v2[0] = {1.0, 0.0};
    v1.insert(v1.end(), V_[0].begin(), V_[0].end());
    v2.insert(v2.end(), V_[1].begin(), V_[1].end());
    v1 = Test.IncreaseDegree(v1,1);
    v2 = Test.IncreaseDegree(v2,1);
    V_Encode = Test.W3jProductCompact(V_output,v1,v2,3);//Embedding as rotation-invariant coding 
    M = Test.DecodeMatrixCompact(v1,v2,3);//Generate an auxiliary matrix for recovering embedded. 
    equivalentFeatures.DecodeMatrixCompact(v1,v2,3);
    std::fill(V_Re.begin(), V_Re.end(),0); 
    for (size_t i = 0; i < M.rows(); ++i) {
        for (size_t j = 0; j < M.cols(); ++j) {
            V_Re[i] += M(i,j) * V_Encode[j];
          }
    }// V_Re=M*V_Encode; 
    norm = 0.0;
    for (size_t i = 0; i < V_Re.size(); ++i) {
        norm += (std::pow(V_Re[i] - V_output[i], 2)).real();
    }
    norm = std::abs(norm);

     std::vector<double> ProductEncode(const std::vector<std::complex<double>>& V1, const std::vector<std::complex<double>>& V2, int n);
 V1 and V2 are spherical harmonic tensors to be encoded, one shell of V1 and one shell of V2 will be Clebsch–Gordan (CG) producted to generate a set of spherical harmonic tensors, the maximum size of each set is n, the return is the list for the norm of each spherical harmonic tensors from these sets.    

    std::vector<double> ProductEncodePairwise(const std::vector<std::complex<double>>& V1, const std::vector<std::complex<double>>& V2, int n);

 V1 and V2 are spherical harmonic tensors to be encoded,  one shell of V1 and one shell of V2 will be Clebsch–Gordan (CG) producted to generate a set of spherical harmonic tensor, the maximum size of each set is n, the return is the list for the norm of each spherical harmonic tensors and pairwise products among spherical harmonic tensors from these sets.        
 
