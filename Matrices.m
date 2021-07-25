function [A_0,A_1,B_0,C_1] = matrices(n)

%  A = -q_1*A_0 - A_1;  B =  q_2*B_0


M_1 = Build_M_1(n);
K_1 = Build_K_1(n);
L_1 = Build_L_1(n);
B_1 = Build_B_1(n);
C_1 = Build_C_1(n);

A_0 = M_1\K_1;
A_1 = M_1\L_1;
B_0 = M_1\B_1;

end

