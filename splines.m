function [Bspl]=splines(degree,N)

n=200;
y=linspace(0,1,n);
intknots=(1:N)/(N+1);
Bspl=bs_function(degree,intknots,[0,1]);
for i=1:N+degree+1
    Bspl{i}=@(s)Bspl{i}(s);
end