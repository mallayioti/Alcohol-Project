function Delta=der_matrix_exp(D,E,q1,s,k)
matrix=[(q1*D+E)*s, D*s; zeros(k), (q1*D+E)*s ];
matrix_2=expm(matrix);
Delta=matrix_2(1:k,(k+1):(2*k));
end
