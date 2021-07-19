function [K,G] = K_matrix(C,D,E,F,dim,p,S,q,bas,mu)

q1=q(1);
q2=q(2);
A=-q1*D - E;
B=q2*F;
pol_=@(u) 0;
psi_=@(s) 0;

for i=1:p
    pol{i}=@(u)bas{i+1}(u);
    psi{i}=@(s)integral(@(u)arrayfun(@(U)C*expm(A*(s-U))*B*pol{i}(U),u),0,s);
    prl_psi_1{i}=@(s)integral(@(u)arrayfun(@(U)C*der_matrix_exp(-D,-E,q1,s-U,dim)*B*pol{i}(U),u),0,s);
    prl_psi_2{i}=@(s)integral(@(u)arrayfun(@(U)(1/q2)*C*expm(A*(s-U))*B*pol{i}(U),u),0,s);   
end

f_mu=@(s)integral(@(u)arrayfun(@(U)C*expm(A*(s-U))*B*mu(U),u),0,s);

for i=1:p
    for j=1:p
        if (i<=j)
            G(i,j)=integral(@(s)(psi{i}(s)*psi{j}(s)),0,S,'ArrayValued', true);
            partial_G_1(i,j)=integral(@(s)(prl_psi_1{i}(s)*psi{j}(s)+psi{i}(s)*prl_psi_1{j}(s)),0,S,'ArrayValued', true);
            partial_G_2(i,j)=2/q2*G(i,j);
        else
            G(i,j)=G(j,i);
            partial_G_1(i,j)=partial_G_1(j,i);
            partial_G_2(i,j)=partial_G_2(j,i);
        end
        j
    end 
    i
end

Z=zeros(1,p);
partial_Z_1=zeros(1,p);
partial_Z_2=zeros(1,p);
for i=1:p
    Z(i)=integral(@(s)(psi{i}(s)*f_mu(s)),0,S,'ArrayValued', true);
    partial_Z_1(i)=integral(@(s)(prl_psi_1{i}(s)*f_mu(s)),0,S,'ArrayValued', true);
    partial_Z_2(i)=(1/q2)*Z(i);
    i
end

G=G./S
Z=Z./S
partial_G_1=partial_G_1./S
partial_G_2=partial_G_2./S
partial_Z_1=partial_Z_1./S
partial_Z_2=partial_Z_2./S

partial_beta_1=inv(G)*transpose(partial_Z_1)-inv(G)*partial_G_1*inv(G)*transpose(Z);
partial_beta_2=inv(G)*transpose(partial_Z_2)-inv(G)*partial_G_2*inv(G)*transpose(Z);
K=([transpose(partial_beta_1);transpose(partial_beta_2)]);