function B=basis_function(degree,i,knots,K)

syms x
if degree==0
    B=1.*((x>=knots(i)) & (x<knots(i+1)));
else
    if (knots(degree+i)-knots(i)==0)
        alpha=0*x;
    else
        alpha=(x-knots(i))/(knots(degree+i)-knots(i));
    end
    if knots(i+degree+1)-knots(i+1)==0
        alpha2=x-x+1;
    else
        alpha2=(knots(i+degree+1)-x)/(knots(i+degree+1)-knots(i+1));
    end
    B=(alpha)*basis_function(degree-1,i,knots)+(alpha2)*basis_function(degree-1,i+1,knots);
end

B = matlabFunction(B);

