function [K,lim_cov,beta_0,beta_n,t_BrAC,e_BrAC,t]=delta_method_general_case_simulations(q_0,sigma_1,sigma_2,C,D,E,F,dim,S,M,m,p,Gamma,q_est)

%Input:
%q_0 = 1x2 vector of true q parameter
%sigma_1 = standard deviation of errors in calibration experiment
%sigma_2 = standard deviation of errors in field experiment
%C,D,E,F = diffusion  odel matrices
%dim = dimension of diffusion model matrices
%S = duration of field experiment
%M = number of experiment repetitions
%m = number of observations in calibration experiment
%p = higher degree of polynomial basis
%Gamma = 2x2 matrix 
%q_est = q estimators

%Output:
%K
%lim_cov = covariance matrix of scaled betas
%beta_0=true parameter beta 
%beta_n = beta estimators
%t_BrAC= true BrAC field experiment
%e_BrAC = nxM matrix with columns the reconstructed BrAC values for each estimated beta
%t = times of field experiment observations


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

mu=@(u) -0.05*u*(u-1);
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
        j
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
    i
end

beta_0=inv(G_0)*transpose(Z_0);
partial_beta_1_0=inv(G_0)*transpose(p_1_Z_0)-inv(G_0)*p_1_G_0*inv(G_0)*transpose(Z_0);
partial_beta_2_0=inv(G_0)*transpose(p_2_Z_0)-inv(G_0)*p_2_G_0*inv(G_0)*transpose(Z_0);
K=([transpose(partial_beta_1_0);transpose(partial_beta_2_0)])
lim_cov=sigma_1^2*transpose(K)*K;


n=floor(m^(3/2))+1;
k=n-1;
t=[0:1/k:S];
epsilon=normrnd(0,sigma_2,[n,1]);
for i=1:M
    sum_G=zeros(p);
    sum_Z=zeros(p,1);
    eta=zeros(p,1);
    Y=zeros(n,1);
    X=zeros(n,2);
    for z=1:n
        Y(z,1)=f_mu(t(z))+epsilon(z);
        X(z,:)=transpose(psi_(q_m(i,1),q_m(i,2),t(z)));
        sum_G=sum_G+psi_(q_m(i,1),q_m(i,2),t(z))*transpose(psi_(q_m(i,1),q_m(i,2),t(z)));
        sum_Z=sum_Z+(psi_(q_m(i,1),q_m(i,2),t(z))*f_mu(t(z)));
        eta=eta+psi_(q_m(i,1),q_m(i,2),t(z))*epsilon(z);
    end
    sum_G=sum_G./n;
    sum_Z=sum_Z./n;
    eta=eta./n;
    beta_n(i,:)=transpose(inv(sum_G)*(sum_Z+eta)); 
end

%reconstructed BrAC curves
true_BrAC=@(u) -0.2*u*(u-1);
for i=1:M
    est_BrAC{i}=@(u) beta_n(i,:)*pol_(u)
    e_BrAC(1,i)=est_BrAC{i}(0);
end

t_BrAC(1)=true_BrAC(0);
for j=1:M
    for i=1:n-1
        t_BrAC(i+1)=true_BrAC(i/k);
        e_BrAC(i+1,j)=est_BrAC{j}(i/k);
    end
end
