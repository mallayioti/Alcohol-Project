% K matrix - HPC

q_0=[1,1];
sigma_1=0.001;
sigma_2=0.001;
C =[0,0,0,1];
D =[ 39.6000,-50.4000,14.4000,-3.6000;-25.2000,46.8000,-28.8000,7.2000;7.2000,-28.8000,46.8000,-25.2000; -3.6000,14.4000,-50.4000,39.6000];
E =[0,0,0,-0.4000;0,0,0,0.8000;0,0,0,-2.8000;0,0,0,10.4000];
F =[10.4000;-2.8000;0.8000;-0.4000];
dim=3;
S=1;
m=60;
p=5;

A=@(q1) -q1*D - E;
B=@(q2) q2*F;
pol_=@(u) 0;
psi_=@(q1,q2,s) 0;

for i=1:(p+1)
    pol{i}=@(u) chebyshevT(i-1,u)/sqrt(integral(@(u)chebyshevT(i-1,u).^2,0,1));
    psi{i}=@(q1,q2,s)integral(@(u)arrayfun(@(U)C*expm(A(q1)*(s-U))*B(q2)*pol{i}(U),u),0,s);
    pol_=@(u)[pol(u);pol{i}(u)];
    psi_=@(q1,q2,s)[psi_(q1,q2,s);psi{i}(q1,q2,s)];
    prl_psi_1{i}=@(q1,q2,s)integral(@(u)arrayfun(@(U)C*der_matrix_exp(-D,-E,q1,s-U,dim+1)*B(q2)*pol{i}(U),u),0,s);
    prl_psi_2{i}=@(s,q1,q2)integral(@(u)arrayfun(@(U)(1/q2)*C*expm(A(q1)*(s-U))*B(q2)*pol{i}(U),u),0,s);   
end

mu=@(u)-0.5*u*(u-1)*(u^2-u+5);
f_mu=@(s)integral(@(u)arrayfun(@(U)C*expm(A(q_0(1))*(s-U))*B(q_0(2))*mu(U),u),0,s);

for i=1:(p+1)
    for j=1:(p+1)
        if (i<=j)
            G{i,j}=@(q1,q2)integral(@(s)(psi{i}(q1,q2,s)*psi{j}(q1,q2,s)),0,S,'ArrayValued', true);
            partial_G_1{i,j}=@(q1,q2)integral(@(s)(prl_psi_1{i}(q1,q2,s)*psi{j}(q1,q2,s)+psi{i}(q1,q2,s)*prl_psi_1{j}(q1,q2,s)),0,S,'ArrayValued', true);
            partial_G_2{i,j}=@(q1,q2)2/q2*G{i,j}(q1,q2);
        else
            G{i,j}=@(q1,q2)G{j,i}(q1,q2);
            partial_G_1{i,j}=@(q1,q2)partial_G_1{j,i}(q1,q2);
            partial_G_2{i,j}=@(q1,q2)partial_G_2{j,i}(q1,q2);
        end
        G_0(i,j)=G{i,j}(q_0(1),q_0(2));
        p_1_G_0(i,j)=partial_G_1{i,j}(q_0(1),q_0(2));
        p_2_G_0(i,j)=partial_G_2{i,j}(q_0(1),q_0(2));
    end 
    i
end

for i=1:(p+1)
    Z{i}=@(q1,q2)integral(@(s)(psi{i}(q1,q2,s)*f_mu(s)),0,S,'ArrayValued', true);
    partial_Z_1{i}=@(q1,q2)integral(@(s)(prl_psi_1{i}(q1,q2,s)*f_mu(s)),0,S,'ArrayValued', true);
    partial_Z_2{i}=@(q1,q2) (1/q2)*Z{i}(q1,q2);
    Z_0(i)=Z{i}(q_0(1),q_0(2));
    p_1_Z_0(i)=partial_Z_1{i}(q_0(1),q_0(2));
    p_2_Z_0(i)=partial_Z_2{i}(q_0(1),q_0(2));
end

beta_0=inv(G_0)*transpose(Z_0);
partial_beta_1_0=inv(G_0)*transpose(p_1_Z_0)-inv(G_0)*p_1_G_0*inv(G_0)*transpose(Z_0);
partial_beta_2_0=inv(G_0)*transpose(p_2_Z_0)-inv(G_0)*p_2_G_0*inv(G_0)*transpose(Z_0);
K=([transpose(partial_beta_1_0);transpose(partial_beta_2_0)])
lim_cov=sigma_1^2*transpose(K)*K;
disp(K)
disp(lim_cov)
