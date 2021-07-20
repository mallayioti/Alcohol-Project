
% construct plot of reconstructed BrAC curves and confidence bounds

fi_=@(u)[Bspl{2}(u),Bspl{3}(u),Bspl{4}(u),Bspl{5}(u),Bspl{6}(u)];
r1=mvnrnd([0,0],eye(2),1000);
r2=mvnrnd([0,0,0,0,0],eye(5),1000);
for i=1:1000
    R1=@(s)sigma*fi_(s)*transpose(K)*Gamma^(-1/2)*transpose(r1(i,:))+sqrt(rho)*sigma*fi_(s)*G^(-1/2)*transpose(r2(i,:));
    R2=@(s)-sigma*fi_(s)*transpose(K)*Gamma^(-1/2)*transpose(r1(i,:))+sqrt(rho)*sigma*fi_(s)*G^(-1/2)*transpose(r2(i,:));
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
for i=1:1000
if RHS_fmin(i)<=0.0665
    counter=counter+1;
end
end

cb=0.0665/sqrt(m);
x=y;
y1 =transpose(spl(:,1));
xconf = [x x(end:-1:1)];          
y1conf = [y1+cb y1(end:-1:1)-cb];
figure
p = fill(xconf,y1conf,'red');
p.FaceColor = [1 0.8 0.8];      
p.EdgeColor = 'none';           
hold on
for i=1:100
    plot(transpose([0:1/199:1]),est_BrAC(:,i),'LineWidth',2)
    hold on
end
scatter(transpose(y),spl(:,1),'LineWidth',2)
title('Reconstructed BrAC curve p=5 with 95% Confidence bounds')
xlabel('Time (in hours)')
ylabel('BrAC')
ylim([0,0.8])
hold off