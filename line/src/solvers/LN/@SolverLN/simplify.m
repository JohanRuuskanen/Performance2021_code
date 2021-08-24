% This function belongs to the simplification part of AutoWfAPH tool, 
% it is responsible for convolving basic patterns.
function [alpha_2,T_2] = Simplify(a1,T1,a2,T2,p1,p2,number)
e = ones(2,1);
zerom = zeros(2,2);
zerom1 = zeros(2,4);

if number == 1 % sequence structure
    alpha_2 = [a1,(1-a1*e)*a2];
    T_2 = [T1,(-T1*e)*a2;zerom,T2];

elseif number == 2 % parallel strucutre
    alpha_2 = [kron(a1,a2),(1-a2*e)*a1,(1-a1*e)*a2];
    Tr1 = [kron(T1,eye(2))+kron(eye(2),T2),kron(eye(2),-T2*e),kron(-T1*e,eye(2))];
    Tr2 = [zerom1,T1,zerom];
    Tr3 = [zerom1,zerom,T2];
    T_2 = [Tr1;Tr2;Tr3];
elseif number == 3 % branch structure
    alpha_2 = [p1*a1,p2*a2];
    T_2 = [T1,zerom;zerom,T2];
else % loop structure
    alpha_2 = a1;
    T_2 = [T1(1,1),T1(1,2);-a1(1)*p1*T1(2,2),T1(2,2)-T1(2,2)*a1(2)*p1];
end
end