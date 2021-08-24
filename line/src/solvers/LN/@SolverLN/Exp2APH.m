% This function is responsible for fitting APH to exponential distribution
function [alpha,T] = Exp2APH(lambda)
m1 = 1/lambda;
m2 = 2/(lambda^2);
m3 = 6/(lambda^3);
[alpha,T] = getAPH2(m1,m2,m3);
end