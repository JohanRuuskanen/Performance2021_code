% This function is used to return CDF points of given APH distribution 
function points = getCdfEntry(alpha,T)
    e = ones(2,1);
    mean = alpha*((-T)^(-1))*e;
    F = [];
    i = 1;
    for x = 0:1:15*mean
      F(i) = 1-alpha*expm(T.*x)*e;
      i = i+1;
    end
    x = 0:1:15*mean;
    points = [F',x'];
end