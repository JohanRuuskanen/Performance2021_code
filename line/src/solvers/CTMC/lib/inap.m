function [x, pi, Q, it, xprev, rcatC3] = inap(ag, options)

actions = ag.actions;
roles = ag.roles;

TOL = options.iter_tol;
MAXITER = options.iter_max;

if nargin<3
    TOL = 1e-4;
elseif nargin<4
    MAXITER = 1000;
end
A=size(roles,1);
ACT=roles(:,1);
PSV=roles(:,2);
M=max(roles(:));
Aa={actions{1:end-1,1}};
Pb={actions{1:end-1,2}};
L={actions{end,1:end}};
N = zeros(1,M); for k=1:M N(k) = length(L{k}); end

x = rand(A,1);
[pi,Q] = compute_equilibrium(x);

it = 1;
xprev = [];
while it <= MAXITER
    % save old iteration result
    piprev = pi;
    xprev(it,1:A) = x;
    
    % for each action
    for a = 1:A
        LAMBDA = pi{ACT(a)}*Aa{a}./pi{ACT(a)};
        LAMBDA = LAMBDA(find(LAMBDA));
        
        %         LAMBDA = zeros(length(P{k}));
        %         for alpha = 1:length(P{k})
        %             for beta = 1:length(P{k})
        %                 LAMBDA(alpha, beta) = P{k,a}(alpha, beta) * pi{k}(alpha) / pi{k}(beta);
        %             end
        %         end
        %         P{j,a}(find(P{j,a})) = mean(LAMBDA(find(LAMBDA)));
        
        x(a) = mean(LAMBDA);
    end
    %x'
    
    [pi,Q] = compute_equilibrium(x);
    
    
    % convergence check
    MAXERR = 0;
    for k=1:M
        MAXERR = max([ MAXERR , norm(pi{k} - piprev{k}, 1) ]);
        if MAXERR < TOL
            %line_printf('requested tolerance reached\n');
            if decideC3()
                return
            end
        end
    end
    it = it + 1;
end
%line_printf('maximum number of iterations reached\n');
if nargout == 6
    rcatC3 = decideC3;
end
return

    function C3 = decideC3()
        % verify if RCAT's 3rd condition is satisfied within tolerance
        C3 = false;
        if max(LAMBDA(find(LAMBDA))) - min(LAMBDA(find(LAMBDA))) < TOL
            C3 = true;
        end
    end

    function [pi,Q] = compute_equilibrium(xsol)
        Q = cell(1,M);
        pi = cell(1,M);
        for k=1:M
            Q{k} = L{k} - diag(L{k}*e(N(k))); %eps action
            for c = 1:A
                if PSV(c) == k
                    Q{k} = Q{k}  + double(xsol(c))*Pb{c} - diag(Pb{c}*e(N(k)));
                elseif ACT(c) == k
                    Q{k} = Q{k}  + Aa{c} - diag(Aa{c}*e(N(k)));
                end
            end
            pi{k} = ctmc_solve(makeinfgen(Q{k}));
        end
    end

end
