% This function is responsible for getting first three moments from CDF.
function [m1,m2,m3] = getMoments(data)
[row,col] = size(data);
m1 = 0; % the first moment
m2 = 0; % the second moment
m3 = 0; % the third moment
for i = 1:1:row-1
    bin1 = ((data(i+1,2)-data(i,2))/2+data(i,2))*((data(i+1,1)-data(i,1)));
    bin2 = ((data(i+1,2)-data(i,2))/2+data(i,2))^2*((data(i+1,1)-data(i,1)));
    bin3 = ((data(i+1,2)-data(i,2))/2+data(i,2))^3*((data(i+1,1)-data(i,1)));
    m1 = m1+bin1;
    m2 = m2+bin2;
    m3 = m3+bin3;
end
end