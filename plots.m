% 1. bivariate normal plot for scaled betas

mean_ = [0 0];
Sigma = sigma^2*transpose(K)*inv(Gamma)*K;
x1 = -0.1:.0002:0.1;
x2 = -0.1:.0002:0.1;
[X1,X2] = meshgrid(x1,x2);
X = [X1(:) X2(:)];
%Evaluate the pdf of the normal distribution at the grid points.
y = mvnpdf(X,mean_,Sigma);
y = reshape(y,length(x2),length(x1));
%level curves
rng('default')  
r= mvnrnd(mean_,Sigma,1000);
pdf=mvnpdf(r,mean_,Sigma);
N1=0;
N2=0;
N3=0;
N4=0;
for i=1:1000
if pdf(i)>5500
    N1=N1+1;
    N2=N2+1;
    N3=N3+1;
    N4=N4+1;
elseif (pdf(i)>2790) & (pdf(i)<=5500)
    N2=N2+1;
    N3=N3+1;
    N4=N4+1;
elseif (pdf(i)>1100) & (pdf(i)<=2790)
    N3=N3+1;
    N4=N4+1;
elseif (pdf(i)>300) & (pdf(i)<=1100)
    N4=N4+1;
end
end
%create a contour plot of the multivariate normal distribution 
hold on
contour(x1,x2,y,[300 1100 2790 5500])
colorbar('Ticks',[300 1100 2790 5500],'TickLabels',{'95%','84%','60%','20%'})
title('beta_n distribution')
xlabel('beta_1')
ylabel('beta_2')
scatter(sqrt(m)*(beta_n(:,1)-beta_0(1)),sqrt(m)*(beta_n(:,2)-beta_0(2)))
hold off

% 2. marginal distribution of betas (expression 64)
x=[-0.1:0.0001:0.1];
mean_=0;
Sigma=sigma^2*transpose(K)*inv(Gamma)*K;
y1 = normpdf(x,mean_,sqrt(Sigma(1,1)));
y2 = normpdf(x,mean_,sqrt(Sigma(2,2)));
subplot(1,2,1)
hold on
title('Marginal distributions of beta_n (expression 64)')
histogram(sqrt(m)*(beta_n(:,1)-beta_0(1)),5);
plot(x,y1)
xlabel('beta_n(1)')
ylabel('density')
hold off
subplot(1,2,2)
hold on
title('Marginal distributions of beta_n (expression 64)')
histogram(sqrt(m)*(beta_n(:,2)-beta_0(2)),5);
plot(x,y2)
xlabel('beta_n(2)')
ylabel('density')
hold off


% 3. Plot for expression 65

x=[-0.1:0.0001:0.1];
mean_=0;
obs=101;
s=(obs-1)/k; %k= number of field experiment observations -1 
Sigma=sigma^2*fi_(s)*transpose(K)*inv(Gamma)*K*transpose(fi_(s));
y = normpdf(x,mean_,sqrt(Sigma(1,1)));
hold on
title('Expression 65 (s=0.216)')
histogram(sqrt(m)*(e_BrAC(obs,:)-t_BrAC(obs)));
plot(x,y)
hold off

% 4. reconstructed BrAC curve with confidence bounds

fi_=@(s)[sqrt(3)*s,sqrt(80)*(s^2-0.75*s)];

for i=1:1000
    R1=@(s)sigma_1*fi_(s)*transpose(K)*Gamma^(-1/2)*r(:,i);
    R2=@(s)-sigma_1*fi_(s)*transpose(K)*Gamma^(-1/2)*r(:,i);
    min1R(i)=fmincon(R1,0.5,[],[],[],[],0,1);
    min2R(i)=fmincon(R2,0.5,[],[],[],[],0,1);
    R1min=R1(min1R(i));
    R2min=R2(min2R(i));
 
    if R1min<R2min
        RHS_fmin(i)=-R1min;
        RHS_fmin_time(i)=min1R(i);
    else
        RHS_fmin(i)=-R2min;
        RHS_fmin_time(i)=min2R(i);
    end
end

counter=0;
for i=1:M
if RHS_fmin(i)<=0.055
    counter=counter+1;
end
end

cb=0.055/sqrt(m);
x = [0:(1/k):1]; %k= total number of field experiment obs -1 (n-1)                 
y1 = [transpose(e_BrAC(:,1))];
y0=t_BrAC;

xconf = [x x(end:-1:1)];          
y1conf = [y1+cb y1(end:-1:1)-cb];
figure
p = fill(xconf,y1conf,'red');
p.FaceColor = [1 0.8 0.8];      
p.EdgeColor = 'none';           
hold on
for i=1:M
    plot(x,transpose(e_BrAC(:,i)))
    hold on
end
plot(x,transpose(t_BrAC),'LineWidth',4)
title('Confidence bounds (without using fmincon)')
xlabel('time(in hours)')
ylabel('BrAC')
hold off