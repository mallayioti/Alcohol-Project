function [q_est,y,Gamma]=q_estimator_simulations(n,m,sigma,q_1,q_2,C,D,E,F,sim_BrAC)

%q parameter least squares estimator for simulations using equally spaced observations

%Input:
%n = number of q estimators to be obtained
%m = number of BrAC/TAC observations -1
%sigma = st.deviation of errors added to calculated TAC 
%q_1,q_2 = true parameters
%C, D, E, F = diffusion model matrices
%sim_BrAC= simulated BrAC - (m+1)x2 vector with time in first column and BrAC value in second column
 

%Output:
%q_est = q estimators
%y = calculated TAC values

 A=-q_1*D - E;
 B=q_2;
 A_inv=inv(A);
 A_q=@(q) -q(1)*D-E;
 A_inv_q=@(q) inv(A_f(q(1)));
 B_q=@(q) q(2);
 t=[1/m:1/m:1];
 J=@(q) 0;
 sum2=@(q)0;
 sum1=0;
 y=zeros(m+1,1);
 
 for z=1:n
    epsilon=normrnd(0,sigma,m+1,1); 
    y(1)=epsilon(1);
    for j=1:m
        s=[0:1/m:t(j)];
        for i=1:j
            mu_1=sim_BrAC(i,2);
            mu_2=sim_BrAC(i+1,2);
            mu=interp1([s(i),s(i+1)],[mu_1,mu_2],(s(i)+s(i+1))/2);
            sum1=sum1+integral(@(u)arrayfun(@(U)C*expm(A*(t(j)-U))*B*F*mu,u),s(i),s(i+1));
            sum2=@(q) sum2(q)-integral(@(u)arrayfun(@(U)C*expm(A_q(q)*(t(j)-U))*B_q(q)*F*mu,u),s(i),s(i+1));
        end
        y(j+1)=sum1+epsilon(j+1);
        J=@(q) J(q)+(y(j+1)+sum2(q))^2;
        sum2=@(q)0;
        sum1=0;
    end
    q_est(z,:)=fmincon(J,[1,1],-eye(2),zeros(2,1),[],[],[0.9,0.9],[1.1,1.1]);
    J=@(q) 0;
 end 
 
 Gamma=inv(cov(q_est)*(m+1)/sigma^2);
    
 