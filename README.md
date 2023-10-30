# Introduction
# Installation
     include("./equivalentFeatures.jl")
     using  .equivalentFeatures 
# Initial
     equivalentFeatures.setN( N_new )
*N_new* is the orbital quantum number *l* of the input spherical harmonic representation *Y(lm)*.   
    If this command is missed, the default orbital quantum number *l* is 3.  
     
     equivalentFeatures.Initial();
