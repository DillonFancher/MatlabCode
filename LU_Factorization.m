function [LU] = LU_Factorization

%A = [2 1 1; 4 3 3; 8 7 9];
A = rand([40,40]);
N = size(A, 2);
L = eye(size(A));
U = A;

%LU Factorization Without Pivoting
for k=1:N
    for i = k+1 : N
        L(i,k)=U(i,k)/U(k,k);
        for j=1:N
             U(i,j)=U(i,j)-L(i,k)*U(k,j);
        end
    end
end

L_norm = norm(L, Inf)
U_Norm = norm(U, Inf)

Error = L*U - A;

Error_Norm = norm(Error, Inf)
 