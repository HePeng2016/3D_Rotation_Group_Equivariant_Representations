module  equivalentFeatures

          using SphericalHarmonicExpansions
          using WignerSymbols
          using DataStructures
          using MultivariatePolynomials
          using LinearAlgebra
          

          struct Index
            m1::Int16
            m2::Int16
          end

          N=3;

          #CGTableV = fill(Vector{Float64}(undef,0),(d1,d2,d3));
          export CGTableC
          export CGTableI
          export W3JTableC
          export W3JTableI
          export SelfPI
          export WignerPI
          export CtoS_C

          function setN(N_new)
            global N = N_new;
          end 


          function Initial()

             d1 = N 
             d2 = N 
             d3 = abs(d1+d2-2)-abs(d1-d2)+1

             global CGTableC = fill(Vector{ Vector{Float64} }(undef,0),(d1,d2,d3));
             global CGTableI = fill(Vector{ Vector{Index} }(undef,0),(d1,d2,d3));
             


             for J1 = 1:d1 
               for J2 = J1:d2
                  for j3 = (J2-J1):(J1+J2-2)
                    j1 = J1-1;
                    j2 = J2-1;
                    J  = j3-(J2-J1)+1;
              
                    # CGTableV[J1,J2,J+1]=Vector{Float64}(undef,2*J+1);
                    CGTableC[J1,J2,J]=Vector{ Vector{Float64} }(undef,2*j3+1);
                    CGTableI[J1,J2,J]=Vector{ Vector{Index} }(undef,2*j3+1);

                    for i = 1:(2*j3+1)  
                      CGTableC[J1,J2,J][i]=Vector{Float64}(undef,0);
                      CGTableI[J1,J2,J][i]=Vector{Index}(undef,0);
                    end 
                   
                    for m1 = -j1:j1
                       for m2 = -j2:j2
                          m3=m2+m1;
                          if abs(m3) <= j3
                             i  = m3+j3+1; 
                             push!(CGTableI[J1,J2,J][i],Index(m1,m2));
                             push!(CGTableC[J1,J2,J][i],clebschgordan(Float64,j1,m1,j2,m2,j3,m1+m2));
                          end
                       end
                    end 
                  end
               end
             end



             d3 = N
             global W3JTableC = fill(Vector{Float64}(undef,0),(d1,d2,d3));
             global W3JTableI = fill(Vector{Index}(undef,0),(d1,d2,d3));

             for J1 = 1:d1
                for J2 = 1:d2 #for J2 = J1:d2
                   for J3 = abs(J2-J1):abs(J1+J2-2)
                      if J3 < d3
                         j1 = J1-1;
                         j2 = J2-1;
                         j3 = J3;
                         J3 = J3+1;

                         W3JTableC[J1,J2,J3]=Vector{Float64}(undef,0);
                         W3JTableI[J1,J2,J3]=Vector{Index}(undef,0);


                         for m1 = -j1:j1
                            for m2 = -j2:j2
                               m3 =-(m1+m2);
                               if abs(m3) <= j3
                                  push!(W3JTableI[J1,J2,J3],Index(m1,m2));
                                  push!(W3JTableC[J1,J2,J3],wigner3j(Float64,j1,j2,j3,m1,m2,m3));
                               end
                            end
                         end

                      end
                   end
                end 
             end


             dictT = Dict{Index,Float64}();
             empty!(dictT);
             List = Any[];
             AlwaysZero=true;

             for J1 = 1:d1 
                for J2 = J1:d2
                   for j3 = (J2-J1):(J1+J2-2)

                      J3 = j3+1;
                      if J1 == J2
                         AlwaysZero=true;
                         for I = 1:length(CGTableC[J1,J2,J3])

                                  empty!(dictT);

                                  for I2 = 1:length(CGTableC[J1,J2,J3][I])
                                     
                                        m1 = max(CGTableI[J1,J2,J3][I][I2].m1,CGTableI[J1,J2,J3][I][I2].m2);
                                        m2 = min(CGTableI[J1,J2,J3][I][I2].m1,CGTableI[J1,J2,J3][I][I2].m2);
                                        M_Index = Index(m1,m2);
                                        if haskey(dictT,M_Index) == false
                                          dictT[M_Index] = 0.0;
                                        end 
                                          dictT[M_Index] = dictT[M_Index]+CGTableC[J1,J2,J3][I][I2];
                                  end

                                  for (k, v) in dictT
                                     if v != 0 
                                        AlwaysZero=false;
                                     end
                                  end
                         end
                      end 

                      if (AlwaysZero != true) && (length(CGTableC[J1,J2,J3])!=0)
                         push!(List,[J1,J2,J3]);
                      end

                   end
                end 
             end

             global SelfPI  = Array{Int}(undef,length(List),3);

             for I = 1:length(List)
               SelfPI[I,:]=List[I];
             end

             List = Any[];
                for J1 = 1:d1
                  for J2 = 1:d2 # for J2 = J1:d2
                      for J3 = abs(J2-J1):abs(J1+J2-2)
                         if J3 < d3
                        
                            J3 = J3+1;
                            if (length(W3JTableI[J1,J2,J3]) != 0)
                               push!(List,[J1,J2,J3]);
                            end
                         end
                      end
                   end 
                end

             global WignerPI = Array{Int}(undef,length(List),3);

             for I = 1:length(List)
               WignerPI[I,:]=List[I];
             end

             @polyvar x y z
             global CtoS_C = Vector{ Vector{Float64} }(undef,N*N);# Cartesian to spherical
             I_= 0;
             for J = 1:N
               j=J-1;
               dict = DictInitial(j);
               for m =-j:j
                 I_ = I_+1; 
                 p=rlylm(j,m,x,y,z);
                 CtoS_C[I_] =  zeros(length(dict));
                 List = collect(terms(p));
                   for I = 1:length(List)
                     i1 = degree(monomial(List[I]),x);
                     i2 = degree(monomial(List[I]),y);
                     i3 = degree(monomial(List[I]),z);
                     CtoS_C[I_][dict[[i1,i2,i3]]]=coefficient(List[I]);
                   end 
               end
             end
          end

         function DictInitial(N)

             s = Stack{Int}();
             push!(s,2);
             Index = 0;
             dict = Dict{Vector{UInt32},Int64}();

             if N<=0
               dict[[0,0,0]]=1.0;
               return dict;  
             end

             while length(s)!= 0
                n = first(s);
                if length(s) < N
                  if n >= 0
                     push!(s,n);
                  else
                     pop!(s);
                     if(length(s)== 0)
                       break;
                     end
                     n=pop!(s);
                     push!(s,n-1);
                  end 
                elseif length(s) == N
                  Index = Index+1;
                  x=0;y=0;z=0;
                  for i in Iterators.reverse(s)
                     if i == 2
                       x=x+1; 
                     elseif i == 1
                       y=y+1;
                     elseif i == 0
                       z=z+1;
                     end
                  end
                     dict[[x,y,z]]=Index;
                  n=pop!(s);
                  n=n-1;
                    if n>=0
                      push!(s,n);
                    else 
                      if(length(s)== 0)
                        break;
                      end
                      n=pop!(s);
                      push!(s,n-1);
                    end 
                end
             end
             return dict; 
         end

         function SelfProduct(V)
             RInvariantV = Complex.(zeros(size(SelfPI)[1]));

             for I =1:size(SelfPI)[1]
                 
                 TempI = CGTableI[SelfPI[I,1],SelfPI[I,2],SelfPI[I,3]];
                 TempC = CGTableC[SelfPI[I,1],SelfPI[I,2],SelfPI[I,3]];
                 J1 = SelfPI[I,1]-1;
                 J2 = SelfPI[I,2]-1;

                 for I_1 = 1:length(TempI)
                   Value = Complex.(0);
                   for I_2 = 1:length(TempI[I_1])
                      m1 = TempI[I_1][I_2].m1;
                      m2 = TempI[I_1][I_2].m2; 
                      I1 = (J1^2 + m1+J1+1);
                      I2 = (J2^2 + m2+J2+1);
                      Value=Value+V[I1]*V[I2]*TempC[I_1][I_2];
                   end
                   RInvariantV[I] = RInvariantV[I] + Value*conj(Value);
                 end 
             end
             return RInvariantV;
         end

         function  W3jProduct(V1,V2,V3)
             RInvariantV = Complex.(zeros(size(WignerPI)[1]));

             for I =1:size(WignerPI)[1]
                TempI = W3JTableI[WignerPI[I,1],WignerPI[I,2],WignerPI[I,3]];
                TempC = W3JTableC[WignerPI[I,1],WignerPI[I,2],WignerPI[I,3]];
                J1 = WignerPI[I,1]-1;
                J2 = WignerPI[I,2]-1;
                J3 = WignerPI[I,3]-1;
                Value = Complex.(0);
                for I_1 = 1:length(TempI)
                   m1 = TempI[I_1].m1;
                   m2 = TempI[I_1].m2;
                   m3 = -(m1+m2);
                   I1 = (J1^2 + m1+J1+1);
                   I2 = (J2^2 + m2+J2+1);
                   I3 = (J3^2 + m3+J3+1);
                   Value = Value + V1[I1]*V2[I2]*V3[I3]*TempC[I_1];
                end
                RInvariantV[I]=Value;
             end 
             return RInvariantV;
         end  



         function  W3jProduct(V1,V2,V3,n)
             
            RInvariantVSize = 0;

            for I in 1:n
               RInvariantVSize = RInvariantVSize +(min((n-1),(n+I-2))-abs(n-I))+1;
            end

            RInvariantV = Complex.(zeros(RInvariantVSize));
            Index =1;
             for I =1:size(WignerPI)[1]
               if WignerPI[I,1] != n 
                  continue; 
               end 
               if WignerPI[I,2] > n || WignerPI[I,3] > n 
                  continue; 
               end 
                TempI = W3JTableI[WignerPI[I,1],WignerPI[I,2],WignerPI[I,3]];
                TempC = W3JTableC[WignerPI[I,1],WignerPI[I,2],WignerPI[I,3]];
                J1 = WignerPI[I,1]-1;
                J2 = WignerPI[I,2]-1;
                J3 = WignerPI[I,3]-1;
                Value = Complex.(0);
                for I_1 = 1:length(TempI)
                   m1 = TempI[I_1].m1;
                   m2 = TempI[I_1].m2;
                   m3 = -(m1+m2);
                   I1 = (J1^2 + m1+J1+1);
                   I2 = (J2^2 + m2+J2+1);
                   I3 = (J3^2 + m3+J3+1);
                   Value = Value + V1[I1]*V2[I2]*V3[I3]*TempC[I_1];
                end
                RInvariantV[Index]=Value;
                Index = Index+1;
             end 
             return RInvariantV;
         end  

         function  CtoS_Encode( V1 ) 
            V2 = zeros(length(CtoS_C));
            Base_ = 0;
            PreviousLength = 0;     
            for I = 1:length(CtoS_C)

               if length(CtoS_C[I]) != PreviousLength
                  Base_ = Base_+length(CtoS_C[I]);
                  PreviousLength = length(CtoS_C[I]);
               end
               Start = Base_ - length(CtoS_C[I]);
               for I2 = 1:length(CtoS_C[I])
                  V2[I] =V2[I]+V1[Start+I2]*CtoS_C[I][I2];
               end  
            end 
            V3 = Complex.(zeros(length(CtoS_C)));

            for I = 1:N
               Start = (I-1)*(I-1);
               for I2 = 1:I
                 if I2 == I
                   V3[Start+I2] = V2[Start+I2];
                 else
                   V3[Start+I2]   = (-V2[Start+I2]*im + V2[Start+2*I-I2])/(2^0.5);
                   V3[Start+2*I-I2] = ((-1)^(I-I2))*(V2[Start+I2]*im + V2[Start+2*I-I2])/(2^0.5);
                 end 
               end 
            end 

            return V3;
         end


         function  CStoRS_Encode( V1 )

            V3 = zeros(length(V1));

            for I = 1:N
               Start = (I-1)*(I-1);
               for I2 = 1:I
                 if I2 == I
                   V3[Start+I2] = real(V1[Start+I2]);
                 else
                   V3[Start+I2]   = ((-1)^(I-I2))*(imag(V1[Start+2*I-I2]))*(2^0.5);
                   V3[Start+2*I-I2] = ((-1)^(I-I2))*(real(V1[Start+2*I-I2]))*(2^0.5);
                 end 
               end 
           end 

           return V3;
         end

         function  CStoRS_Encode(V1,n)
            V3 = zeros(length(V1));
            for I2 = 1:n
                 if I2 == n
                   V3[I2] = real(V1[I2]);
                 else
                   V3[I2]   = ((-1)^(n-I2))*(imag(V1[2*n-I2]))*(2^0.5);
                   V3[2*n-I2] = ((-1)^(n-I2))*(real(V1[2*n-I2]))*(2^0.5);
                 end 
            end 
            return V3;
         end

         function  RStoCS_Encode(V1,n)
            V3 = Complex.(zeros(length(V1))); 
            for I2 = 1:n
                 if I2 == n
                   V3[I2] = V1[I2];
                 else
                   V3[I2]     = (-V1[I2]*im + V1[2*n-I2])/(2^0.5);
                   V3[2*n-I2] = ((-1)^(n-I2))*(V1[I2]*im + V1[2*n-I2])/(2^0.5);
                 end 
            end 
            return V3; 
         end 

         function  RStoCS_Encode( V1 )

           V3 = Complex.(zeros(length(V1))); 

           for I = 1:N
               Start = (I-1)*(I-1);
               for I2 = 1:I
                 if I2 == I
                   V3[Start+I2] = V1[Start+I2];
                 else
                   V3[Start+I2]   = (-V1[Start+I2]*im + V1[Start+2*I-I2])/(2^0.5);
                   V3[Start+2*I-I2] = ((-1)^(I-I2))*(V1[Start+I2]*im + V1[Start+2*I-I2])/(2^0.5);
                 end 
               end 
           end 

           return V3;
         end






         function DecodeMatrix(V2,V3)
   
               LTm =  Complex.(zeros(size(WignerPI)[1],N*N)); 
               pseudoInput = zeros(size(WignerPI)[1]);

               for I =1:size(WignerPI)[1]
                  TempI = W3JTableI[WignerPI[I,1],WignerPI[I,2],WignerPI[I,3]];
                  TempC = W3JTableC[WignerPI[I,1],WignerPI[I,2],WignerPI[I,3]];
                  J1 = WignerPI[I,1]-1;
                  J2 = WignerPI[I,2]-1;
                  J3 = WignerPI[I,3]-1;
                     for I_1 = 1:length(TempI)
                        m1 = TempI[I_1].m1;
                        m2 = TempI[I_1].m2;
                        m3 = -(m1+m2);
                        I1 = (J1^2 + m1+J1+1);
                        I2 = (J2^2 + m2+J2+1);
                        I3 = (J3^2 + m3+J3+1);
                        LTm[I,I1] = LTm[I,I1]+V2[I2]*V3[I3]*TempC[I_1];
                     end
               end

               for I =1:size(WignerPI)[1]
                  pseudoInput[I]=(LTm[I,:]'LTm[I,:]).re;
               end

               G_inv=LinearAlgebra.pinv(LTm); 
               return pseudoInput,G_inv;
         end







         function DecodeMatrix(V2,V3,n)

            VSize = 0;
            for I in 1:n
               VSize = VSize +(min((n-1),(n+I-2))-abs(n-I))+1;
            end
   
            LTm =  Complex.(zeros(VSize,2*(n-1)+1)); 
            # pseudoInput = zeros(VSize);
               Index = 1;  
               for I =1:size(WignerPI)[1]
                  if WignerPI[I,1] != n
                     continue; 
                  end 
                  if WignerPI[I,2] > n || WignerPI[I,3] > n 
                     continue; 
                  end 
                  TempI = W3JTableI[WignerPI[I,1],WignerPI[I,2],WignerPI[I,3]];
                  TempC = W3JTableC[WignerPI[I,1],WignerPI[I,2],WignerPI[I,3]];
                  J1 = WignerPI[I,1]-1;
                  J2 = WignerPI[I,2]-1;
                  J3 = WignerPI[I,3]-1;
                     for I_1 = 1:length(TempI)
                        m1 = TempI[I_1].m1;
                        m2 = TempI[I_1].m2;
                        m3 = -(m1+m2);
                        I1 = (m1+J1+1);
                        I2 = (J2^2 + m2+J2+1);
                        I3 = (J3^2 + m3+J3+1);
                        LTm[Index,I1] = LTm[Index,I1]+V2[I2]*V3[I3]*TempC[I_1];
                     end
                  Index = Index+1;
               end

              #for I =1:VSize
              #    pseudoInput[I]=(LTm[I,:]'LTm[I,:]).re;
              #end

              G_inv=LinearAlgebra.pinv(LTm); 
              return G_inv;
         end
   end


