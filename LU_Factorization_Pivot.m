function LU_Factorization(A)

%A = [2 1 1 4; 3 2 0 1; 4 3 3 3; 5 1 2 0];
% = rand([160, 160]);
%A(1,1) = 10e-7;
N = size(A, 2);
L = zeros(size(A));
U = A;
Piv = [1:N];


%LU Factorization Without Pivoting
for k=1:N-1  
    
    %Check current pivot against potentially greater values in the column
    %and note the value of the new pivot
    Curr_Piv = abs(U(k,k));
    New_Piv = k;
    for p = k+1:N
       if(Curr_Piv < abs(U(p,k)))
           Curr_Piv = abs(U(p,k));
           New_Piv = p;
       end
    end
   
    %Initialize PIVOTER, so only 1 Row swap happens per iteration
    
    PIVOTER = eye(N);
    tempPiv = 0;
        if(New_Piv ~= k)
            tempPiv = Piv(k);
            Piv(k) = Piv(New_Piv);
            Piv(New_Piv) = tempPiv;
            %Create Pivot Matrix to Swap Rows of L
            PIVOTER(k,k) = 0;
            PIVOTER(k,New_Piv) = 1;
            PIVOTER(New_Piv,New_Piv) = 0;
            PIVOTER(New_Piv,k) = 1;
            
            %Keep track of pivot changes
            %Piv(k,k) = 0;
            %Piv(k,New_Piv) = 1;
            %Piv(New_Piv,New_Piv) = 0;
            %Piv(New_Piv,k) = 1;
            
            %Use PIVOTER to swap rows of L
            L = PIVOTER*L;
            
            %Use PIVOTER to swap rows of U
            U = PIVOTER*U;
        end
    
    %Handle the elimination of all column values below current pivot and
    %update in the L matrix and U matrix
    for i = k+1 : N
        L(i,k)=U(i,k)/U(k,k);
        for j=1:N
             U(i,j)=U(i,j)-L(i,k)*U(k,j);
        end
    end
    
end

%Add the 1 povit along the main diagonal of L
for lpiv = 1:N
   L(lpiv,lpiv) = 1; 
end

Pivot = zeros(size(A));
for npiv = 1:N
    Pivot(npiv, Piv(npiv)) = 1;
end

%Note the infinity norms for L, U and the error
L_norm = norm(L, Inf)
U_Norm = norm(U, Inf)

Error = L*U - Pivot*A;
Error_Norm = norm(Error, Inf)
 