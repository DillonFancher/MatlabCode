function LU_Factorization(A)

%Written by: Dillon Fancher
%2/20/14

%This code performs an LU factorization with partial Pivoting on a square
%matrix A, that is an element of the Real(NxN) plane. It computes the
%pivoted version of L in real time, whilst computing U in place so no
%pivoting is actually performed in the LU-Factorization for loop, for 
%efficiency purposes. The goal of using a pivot vector and looping over its
%entries inside the for loops (i.e. U(Piv(i), j) instead of U(i, j)) was to
%help memory in not forming the Permutation Matrix, however Due to the need
%for specific analysis of the Infinity Norms of U, L and the Error,
%(Error = LU - A) we form the Permutation matrix at the end for easy
%analysis, and it is not detrimental because we are not working with
%extremely large matrices here. Also I computed U and L seperately instead
%of ontop of A in order for the analysis of the norms and error to be easy
%to handle after the Factorization instead of pulling out the entries of A
%and throwing them into L and U afterwards.

%start--------------------------------------------------------------------
    %Initialize Parameters to Start with
    A = rand([160, 160]);
    %Consider the worst case scenario of having an extremely small initial
    %pivot, so we assign A(1,1) to be tiny
    A(1,1) = 10e-7;
    N = size(A, 2);
    L = eye(size(A));
    U = A;
    %This is the pivot vector we will be using in order to calculate U in
    %place
    Piv = [1:N];

%LU Factorization Without Pivoting
for k=1:N-1  
    
    %Check current pivot against potentially greater values in the column
    %and note the value of the new pivot
    Curr_Piv = abs(U(k,k));
    New_Piv = k;
    for p = k+1:N
       if(Curr_Piv < abs(U(Piv(p),k)))
           Curr_Piv = abs(U(Piv(p),k));
           New_Piv = p;
       end
    end
 
    %Adjust the pivot vector and entries of L
    if(New_Piv ~= k)
        tempPiv = Piv(k);
        Piv(k) = Piv(New_Piv);
        Piv(New_Piv) = tempPiv;
        %Update L because we go through L in place during the
        %elimination of column entries below
        for L1_ent = 1: k-1
             tempL = L(k,L1_ent);
             L(k, L1_ent) = L(New_Piv, L1_ent);
             L(New_Piv, L1_ent) = tempL;
        end
    end
        
    %Handle the elimination of all column values below current pivot and
    %update in the L matrix and U matrix  
    for i = k+1 : N
        L(i,k)=U(Piv(i),k)/U(Piv(k),k);
        for j=1:N
             U(Piv(i),j)=U(Piv(i),j)-L(i,k)*U(Piv(k),j);
        end
    end
    
end

%Calculate the Permutation Matrix in order to calculate L*U properly and in
%order to re-order A. This was in essence an extra step because the
%computational efficiency of the code relied on not having to actually
%having to form the pivot matrix, but for the sake of analysis it was
%necessary.
Pivot = zeros(size(A));
for npiv = 1:N
    Pivot(npiv, Piv(npiv)) = 1;
end

%Note the infinity norms for L, U and the error
L_norm = norm(L, Inf)
U_Norm = norm(U, Inf)
Error = L*Pivot*U - Pivot*A;
Error_Norm = norm(Error, Inf)

%end----------------------------------------------------------------------
 