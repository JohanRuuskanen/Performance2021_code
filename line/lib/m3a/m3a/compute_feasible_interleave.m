function [r,transpose] = compute_feasible_interleave(uv, dv)
% feasible interleave MILP

    tol = 1e-6; % tolerance for strict inequalities

    k = length(uv);

    s = sdpvar(k,1);
    c = sdpvar(k,1);
    alpha = sdpvar(k,1);
    beta = sdpvar(k,1);
    x = binvar(k,1);

    feas_s = -uv .* x + uv + tol <= s <= dv - uv .* x - tol;
    feas_c = uv.* x + tol <= c <= dv - uv + uv .* x - tol;
    fit_d = s + c == dv;
    pos_s = s >= tol;
    pos_c = c >= tol;
    pos_alpha = alpha >= tol;
    pos_beta = beta >= tol;

    constraints = [feas_s, feas_c, fit_d, ...
                   pos_s, pos_c, pos_alpha, pos_beta];

    for i = 1:k
        constraints = [constraints, sum(alpha(i:k)) == s(i)];
    end
    for i = 1:k
        constraints = [constraints, sum(beta(1:i)) == c(i)];
    end

    sdpsettings('solver','gurobi');
    diag = solvesdp(constraints, 0);

    if diag.problem ~= 0
        error('Unable to compute feasible interleaving: %s', yalmiperror(diag.problem));
    end

    r = zeros(2,k);
    r(1,:) = double(s);
    r(2,:) = double(c);

    transpose = double(x);

end