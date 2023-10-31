# Introduction
# Installation
     include("./equivalentFeatures.jl")
     using  .equivalentFeatures 
# Initial
     equivalentFeatures.setN( N_new )
*N_new* is the orbital quantum number *l* of the input spherical harmonic representation *Y(lm)*.   
    If this command is missing, then the default orbital quantum number *l* is 3.  
     
     equivalentFeatures.Initial();
   Once this function had been carried out, the tables for the Clebsch Gordan coefficients and the Winger 3j coefficients were generated. And the data structures for storing the relations between coefficients and variables were generated.
# Usage
    equivalentFeatures.CtoS_Encode( V1 )

