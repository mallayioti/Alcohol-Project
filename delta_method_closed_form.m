function [K,lim_cov,beta_n,t]=delta_method_closed_form(q0,sigma_1,sigma_2,S,M,m,Gamma)

%Input:
%q_0 = 1x2 vector of true q parameter
%sigma_1 = standard deviation of errors in calibration experiment
%sigma_2 = standard deviation of errors in field experiment
%S = duration of field experiment
%M = number of experiment repetitions
%m = number of observations in calibration experiment
%Gamma = 2x2 matrix 

%Output:
%K
%lim_cov = covariance matrix of scaled betas
%beta_n = beta estimators
%t= times of field experiment observations


f1=@(q1,q2) sqrt(3)*((11*q2)/(30*q1^2) - (3*q2)/(20*q1) + (3*q2*(1/q1^2 - (exp(-q1)*(q1 + 1))/q1^2))/(5*q1^2) + (3*q2*(exp(-q1) - 1))/(5*q1^3) + (3*q2*(exp(-1) - 1))/(5*q1^2) - (q2*(2/q1^3 - (exp(-q1)*(q1^2 + 2*q1 + 2))/q1^3))/(5*q1^2) - (3*q2*(2*exp(-1) - 1))/(5*q1) - (3*q2*(exp(- q1 - 1) - 1))/(5*q1^2*(q1 + 1)));
f2=@(q1,q2) (309*5^(1/2)*q2)/(100*q1) - (29*5^(1/2)*q2)/(10*q1^2) + (11*5^(1/2)*q2)/(3*q1^3) + (3*5^(1/2)*q2)/q1^4 - (18*5^(1/2)*q2)/(5*q1^5) + (16*5^(1/2)*q2)/(5*q1^6) - (9*5^(1/2)*q2)/(5*(q1^3 + q1^2)) - (24*5^(1/2)*q2)/(5*(q1^4 + q1^3)) - (3*5^(1/2)*q2*exp(-q1))/(5*q1^3) - (5^(1/2)*q2*exp(-q1))/q1^4 + (2*5^(1/2)*q2*exp(-q1))/(5*q1^5) - (16*5^(1/2)*q2*exp(-q1))/(5*q1^6) - (42*5^(1/2)*q2*exp(-1))/(5*q1) + (39*5^(1/2)*q2*exp(-1))/(5*q1^2) - (24*5^(1/2)*q2*exp(-1))/(5*q1^3) + (9*5^(1/2)*q2*exp(-q1)*exp(-1))/(5*(q1^3 + q1^2)) + (24*5^(1/2)*q2*exp(-q1)*exp(-1))/(5*(q1^4 + q1^3));
f3=@(q1,q2) 3*(q2^2*(q1 - exp(-2*q1)/2 - 2*q1*exp(-q1) - q1^2 + q1^3/3 + 1/2))/q1^5;
f4=@(q1,q2) sqrt(3)*(q1*(8*5^(1/2)*q2^2*exp(-q1) - (19*5^(1/2)*q2^2)/2 + (3*5^(1/2)*q2^2*exp(-2*q1))/2) + 4*5^(1/2)*q2^2 + q1^2*(5*5^(1/2)*q2^2 + 2*5^(1/2)*q2^2*exp(-q1)) - 5^(1/2)*q1^3*q2^2 - 8*5^(1/2)*q2^2*exp(-q1) + 4*5^(1/2)*q2^2*exp(-2*q1))/q1^6;
f5=@(q1,q2) (q2^2*exp(-2*q1)*(3840*exp(q1) - 2880*exp(2*q1) - 720*q1 + 1200*q1*exp(2*q1) + 480*q1^2*exp(q1) + 180*q1^3*exp(q1) - 345*q1^2*exp(2*q1) + 110*q1^3*exp(2*q1) - 30*q1^4*exp(2*q1) + 6*q1^5*exp(2*q1) + 1440*q1*exp(q1) - 135*q1^2 - 960))/(6*q1^7);
f6=@(q1,q2) sqrt(3)*((7*q2)/(15*q1^3) - (3*q2)/(5*(q1^4 + 2*q1^3 + q1^2)) - (9*q2)/(20*q1^2) - (6*q2)/(5*(q1^4 + q1^3)) + (9*q2)/(5*q1^4) - (12*q2)/(5*q1^5) + (2*q2)/q1^6 + (6*q2*exp(-1))/(5*q1^2) - (6*q2*exp(-1))/(5*q1^3) - (q2*exp(-q1))/(5*q1^3) - (2*q2*exp(-q1))/(5*q1^4) + (2*q2*exp(-q1))/(5*q1^5) - (2*q2*exp(-q1))/q1^6 + (3*q2*exp(-q1)*exp(-1))/(5*(q1^3 + 2*q1^2 + q1)) + (6*q2*exp(-q1)*exp(-1))/(5*(q1^4 + q1^3)) + (6*q2*exp(-q1)*exp(-1))/(5*(q1^4 + 2*q1^3 + q1^2)));
f7=@(q1,q2) (29*5^(1/2)*q2)/(5*q1^3) - (309*5^(1/2)*q2)/(100*q1^2) - (11*5^(1/2)*q2)/q1^4 - (12*5^(1/2)*q2)/q1^5 + (18*5^(1/2)*q2)/q1^6 - (96*5^(1/2)*q2)/(5*q1^7) + (18*5^(1/2)*q2)/(5*(q1^4 + q1^3)) + (72*5^(1/2)*q2)/(5*(q1^5 + q1^4)) + (9*5^(1/2)*q2)/(5*(q1^4 + 2*q1^3 + q1^2)) + (24*5^(1/2)*q2)/(5*(q1^5 + 2*q1^4 + q1^3)) + (3*5^(1/2)*q2*exp(-q1))/(5*q1^3) + (14*5^(1/2)*q2*exp(-q1))/(5*q1^4) + (18*5^(1/2)*q2*exp(-q1))/(5*q1^5) + (6*5^(1/2)*q2*exp(-q1))/(5*q1^6) + (96*5^(1/2)*q2*exp(-q1))/(5*q1^7) + (42*5^(1/2)*q2*exp(-1))/(5*q1^2) - (78*5^(1/2)*q2*exp(-1))/(5*q1^3) + (72*5^(1/2)*q2*exp(-1))/(5*q1^4) - (18*5^(1/2)*q2*exp(-q1)*exp(-1))/(5*(q1^4 + q1^3)) - (72*5^(1/2)*q2*exp(-q1)*exp(-1))/(5*(q1^5 + q1^4)) - (42*5^(1/2)*q2*exp(-q1)*exp(-1))/(5*(q1^4 + 2*q1^3 + q1^2)) - (48*5^(1/2)*q2*exp(-q1)*exp(-1))/(5*(q1^5 + 2*q1^4 + q1^3)) - (9*5^(1/2)*q2*exp(-q1)*exp(-1))/(5*(q1^3 + 2*q1^2 + q1));
f8=@(q1,q2) 2*3*(q2^2*exp(-2*q1)*(6*q1 - 15*exp(2*q1) - 24*q1*exp(2*q1) + 12*q1^2*exp(q1) + 18*q1^2*exp(2*q1) - 4*q1^3*exp(2*q1) + 48*q1*exp(q1) + 15))/(12*q1^6);
f9=@(q1,q2) sqrt(3)*-(q1*(32*5^(1/2)*q2^2*exp(-q1) - (95*5^(1/2)*q2^2)/2 + (31*5^(1/2)*q2^2*exp(-2*q1))/2) + 24*5^(1/2)*q2^2 + q1^2*(20*5^(1/2)*q2^2 + 16*5^(1/2)*q2^2*exp(-q1) + 3*5^(1/2)*q2^2*exp(-2*q1)) - q1^3*(3*5^(1/2)*q2^2 - 2*5^(1/2)*q2^2*exp(-q1)) - 48*5^(1/2)*q2^2*exp(-q1) + 24*5^(1/2)*q2^2*exp(-2*q1))/q1^7;
f10=@(q1,q2) -(q2^2*exp(-2*q1)*(26880*exp(q1) - 20160*exp(2*q1) - 6240*q1 + 7200*q1*exp(2*q1) + 3840*q1^2*exp(q1) + 1200*q1^3*exp(q1) + 180*q1^4*exp(q1) - 1725*q1^2*exp(2*q1) + 440*q1^3*exp(2*q1) - 90*q1^4*exp(2*q1) + 12*q1^5*exp(2*q1) + 12480*q1*exp(q1) - 2115*q1^2 - 270*q1^3 - 6720))/(6*q1^8);
f_mu=@(s) (3*s)/5 + (3*exp(-s))/5 - s^2/5 - 3/5;
psi_1=@(q1,q2,s) sqrt(3)*(q2*(exp(-q1*s) + q1*s - 1))/q1^2;
psi_2=@(q1,q2,s) -(5^(1/2)*q2*(8*exp(-q1*s) - 8) + 5^(1/2)*q1^2*q2*(- 4*s^2 + 3*s) + 5^(1/2)*q1*q2*(8*s + 3*exp(-q1*s) - 3))/q1^3;
psi=@(q1,q2,s)[psi_1(q1,q2,s);psi_2(q1,q2,s)];


Gn=@(q1,q2)[f3(q1,q2),f4(q1,q2);f4(q1,q2),f5(q1,q2)];
Zn=@(q1,q2)[f1(q1,q2);f2(q1,q2)];
G0=Gn(q0(1),q0(2));
Z0=Zn(q0(1),q0(2));
par_1_Z0=transpose([f6(q0(1),q0(2)),f7(q0(1),q0(2))]);
par_2_Z0=Z0/q0(2);
par_1_G0=[f8(q0(1),q0(2)),f9(q0(1),q0(2));f9(q0(1),q0(2)),f10(q0(1),q0(2))];
par_2_G0=2*G0/q0(2);
beta_0=inv(G0)*Z0;
par_1_beta=(inv(G0)*par_1_Z0)-(inv(G0)*par_1_G0*inv(G0)*Z0);
par_2_beta=(inv(G0)*par_2_Z0)-(inv(G0)*par_2_G0*inv(G0)*Z0);
K=[transpose(par_1_beta);transpose(par_2_beta)];
lim_cov=sigma_1^2*transpose(K)*K;

n=floor(m^(3/2))+1;
k=n-1;
t=[0:1/k:S];
epsilon=normrnd(0,sigma_2,[n,1]);

for i=1:M
    q_m(i,:)=mvnrnd([1,1],(1/m)*sigma_1^2*inv(Gamma),1);
    sum_G=zeros(2);
    sum_Z=zeros(2,1);
    eta=zeros(2,1);
    Y=zeros(n,1);
    X=zeros(n,2);
    for z=1:n
        Y(z,1)=f_mu(t(z))+epsilon(z);
        X(z,:)=transpose(psi(q_m(i,1),q_m(i,2),t(z)));
        sum_G=sum_G+psi(q_m(i,1),q_m(i,2),t(z))*transpose(psi(q_m(i,1),q_m(i,2),t(z)));
        sum_Z=sum_Z+(psi(q_m(i,1),q_m(i,2),t(z))*f_mu(t(z)));
        eta=eta+psi(q_m(i,1),q_m(i,2),t(z))*epsilon(z);
    end
    sum_G=sum_G./n;
    sum_Z=sum_Z./n;
    eta=eta./n;
    beta_n(i,:)=transpose(inv(sum_G)*(sum_Z+eta)); 
end