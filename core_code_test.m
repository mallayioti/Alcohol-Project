function [est_BrAC,beta_n,TAC,q_est]=core_code_test(q,sigma_c,sigma_f,Gamma,m,n,dim,C,D,E,F,mu,S,T,bas,p,exp_num)

% q=[0.63,0.78];
% sigma_c=0.003;
% sigma_f=0;
% Gamma=0.01*eye(2);
% m=60;
% n=floor(m^(3/2))+1;
% dim=4;
% [D,E,F,C]=matrices(dim-1);
% mu=@(u)Bspl{3}(u);
% S=1;
% T=1;
% bas=Bspl;
% p=5;
%exp_num=100;

inv_G=inv(Gamma);
A= -q(1)*D - E;
B= q(2)*F;
TAC(:,1)=[0:S/(n-1):S];
t=TAC(:,1);

    for w=1:exp_num
        w
        %generate q_m and new field session 
        q_est(w,:)=mvnrnd(q,sigma_c^2*inv_G/m);
        A_est= -q_est(w,1)*D - E;
        B_est= q_est(w,2)*F;
        epsilon=normrnd(0,sigma_f,[n,1]);
        for j=1:n
            TAC(j,2)=integral(@(u)arrayfun(@(U)C*expm(A*(TAC(j,1)-U))*B*mu(U),u),0,TAC(j,1))+epsilon(j);
        end
    
        pol_=@(u)0;
        psi_=@(u)0;
        for i=1:p
            pol{i}=@(u)bas{i+1}(u);
            psi{i}=@(s)integral(@(u)arrayfun(@(U)C*expm(A_est*(s-U))*B_est*pol{i}(U),u),0,s);
            pol_=@(u)[pol_(u);pol{i}(u)];
            psi_=@(s)[psi_(s);psi{i}(s)];   
        end
        
        %calculate beta 
        clear v X
        for i=1:length(TAC)
            v=psi_(t(i));
            v=v(2:(p+1));
            X(i,:)=transpose(v);
        end
        beta_n(w,:)=inv(transpose(X)*X)*transpose(X)*TAC(:,2);
        
        %calculate reconstructed BrAC curve
        est_BrAC_=@(u) 0;

        for i=1:p
            est_BrAC_=@(u) est_BrAC_(u)+beta_n(w,i)*pol{i}(u);
        end
        
        for j=1:200
            est_BrAC(j,w)=est_BrAC_(T*(j-1)/199);
        end
    end
