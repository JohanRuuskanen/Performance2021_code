% This function belongs to the simplification part of AutoWfAPH tool, 
% it is capable of convolving any number of activities in parallel structure.
function [alpha_2, T_2] = convAP(s)
[x,y] = Simplify(s{1},s{2},s{3},s{4},1,1,2);
[xx,yy] = reduceOrder(x,y);
for i = 1:1:(length(s)/2)-2
    [x,y] = Simplify(xx,yy,s{3+2*i},s{4+2*i},1,1,2);
    [xx,yy] = reduceOrder(x,y);
end
alpha_2 = xx;
T_2 = yy;
end