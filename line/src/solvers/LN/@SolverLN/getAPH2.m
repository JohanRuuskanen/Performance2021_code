% This function is responsible for fitting APH distributions to given moments.
function [alpha_2, T_2] = getAPH2(m1,m2,m3)

if m2<1.5*m1^2
    m2 = 1.5*m1^2;
else
    m2 = m2;
end

scv = (m2/(m1^2))-1;
if scv == 0.5
    satisfy = 2;
elseif 0.5<scv && scv<1 
    if 6*(m1^3)*scv<m3 && m3<3*(m1^3)*(3*scv-1+sqrt(2)*(1-scv)^(3/2))
        satisfy = 1;
        m2 = 3*(m1^3)*(3*scv-1+sqrt(2)*(1-scv)^(3/2));
        m3 = 6*(m1^3)*scv;
    elseif m3<min(3*(m1^3)*(3*scv-1+sqrt(2)*(1-scv)^(3/2)),6*(m1^3)*scv)
        satisfy = 0;
    elseif m3>max(6*(m1^3)*scv,3*(m1^3)*(3*scv-1+sqrt(2)*(1-scv)^(3/2)))
        satisfy = 0;
    else
        satisfy = 1;
        m3 = m3;
    end
elseif scv == 1
    satisfy = 3;
    % one dimensional exponential distribution with lambda = 1/m1
elseif scv>1 
    satisfy = 1;
    if m3<=(3/2)*m1^3*(1+scv)^2
        m3 = (3/2)*m1^3*(1+scv)^2;
    else 
        m3 = m3;
    end
else
    satisfy = 1;
    m3 = m3;
end

if satisfy == 2
    c = 0.75*m1^4;
    d = 0.5*m1^2;
    b = 1.5*m1^3;
    a = 0;
    mu = (-b+6*m1*d+sqrt(a))/(b+sqrt(a));
    lambda1 = (b-sqrt(a))/c;
    lambda2 = (b+sqrt(a))/c;
elseif satisfy == 1
    c = 3*m2^2-2*m1*m3;
    d = 2*m1^2-m2;
    b = 3*m1*m2-m3;
    a = b^2-6*c*d;

    if c>0
        mu = (-b+6*m1*d+sqrt(a))/(b+sqrt(a));
        lambda1 = (b-sqrt(a))/c;
        lambda2 = (b+sqrt(a))/c;
    elseif c<0
        mu = (b-6*m1*d+sqrt(a))/(-b+sqrt(a));
        lambda1 = (b+sqrt(a))/c;
        lambda2 = (b-sqrt(a))/c;
    else
        mu = 1/(2*scv);
        lambda1 = 1/(scv*m1);
        lambda2 = 2/m1;
        % one dimensional exponential distribution with lambda = 1/m1
    end
else
    mu = 1/(2*scv);
    lambda1 = 1/(scv*m1);
    lambda2 = 2/m1;
end
alpha_2 = [mu,1-mu];
T_2 = [-lambda1,lambda1;0,-lambda2];
end