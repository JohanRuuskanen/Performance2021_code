function [r,transpose,order] = compute_feasible_interleave_reorder(uv, dv)
        
    tol = 1e-6; % tolerance for strict inequalities

    k = length(uv);

    s = sdpvar(k,1);
    c = sdpvar(k,1);
    %alpha = sdpvar(k,1);
    %beta = sdpvar(k,1);
    x = binvar(k,1); % transposition
    p = intvar(k,1); % position
    z = cell(k,k);
    for i = 1:k
        for j = 1:k
            z{i,j} = binvar(1,1);
        end
    end
    %z = binvar(k,k); % z(i,j) = 1 if and only if p(i) <= p(j)

    feas_s = -uv .* x + uv + tol <= s <= dv - uv .* x - tol;
    feas_c = uv.* x + tol <= c <= dv - uv + uv .* x - tol;
    fit_d = s + c == dv;
    pos_s = s >= tol;
    pos_c = c >= tol;
    %pos_alpha = alpha >= tol;
    %pos_beta = beta >= tol;
    min_p = p >= 1;
    max_p = p <= k;

    constraints = [feas_s, feas_c, fit_d, ...
                   pos_s, pos_c, ...
                   %pos_alpha, pos_beta, ...
                   min_p, max_p ...
                  ];
    for i = 1:k
        for j = (i+1):k
            constraints = [constraints, p(i) ~= p(j)];
        end
    end

    %for i = 1:k
        %constraints = [constraints, sum(alpha(i:k)) == s(i)];
    %end
    %for i = 1:k
        %constraints = [constraints, sum(beta(1:i)) == c(i)];
    %end
    
    L = 1e6;
    for i = 1:k
        for j = 1:k
            constraints = [constraints, s(i) >= s(j) - (1-z{i,j}) * L];
            constraints = [constraints, c(i) <= c(j) + (1-z{i,j}) * L];
            constraints = [constraints, L*z{i,j} >= p(j) - p(i)];
        end
    end
    
    for i = 1:k
        constraints = [constraints, z{i,i} == 1];
        for j = (i+1):k
            constraints = [constraints, z{i,j}+z{j,i} == 1];
        end
    end

    %options = sdpsettings('solver','cbc');
    options = sdpsettings();
    diag = solvesdp(constraints, 0, options);

    if diag.problem ~= 0
        error('Unable to compute feasible interleaving: %s', yalmiperror(diag.problem));
    end

    r = zeros(2,k);
    r(1,:) = double(s)
    r(2,:) = double(c)
    order = zeros(k,1);
    for j = 1:k
        order(round(double(p(j)))) = j;
    end
    order
    precede = zeros(k,k);
    for i = 1:k
        for j = 1:k
            precede(i,j) = round(double(z{i,j}));
        end
    end
    precede

    transpose = double(x);

end