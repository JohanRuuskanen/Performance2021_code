function [G,lG]= pfqn_mmint2(L,N,Z)
% [G,LOGG] = PFQN_PNC2(L,N,Z)

nonzeroClasses = find(N);
% repairmen integration
order = 15;
% below we use a variable substitution u->u^2 as it seems to be
% numerically better
f= @(u) (2*u'.*exp(-(u.^2)').*prod((Z(nonzeroClasses)+L(nonzeroClasses).*repmat(u(:).^2,1,length(nonzeroClasses))).^N(nonzeroClasses),2))';
%f= @(u) (exp(-u').*prod((Z(nonzeroClasses)+L(nonzeroClasses).*repmat(u(:),1,length(nonzeroClasses))).^N(nonzeroClasses),2))';

p = 1-10^-order;
exp1prctile = -log(1-p)/1; % cutoff for exponential term
w = warning ;
warning off;
lG = log(integral(f,0,exp1prctile,'AbsTol',10^-order)) - sum(factln(N));
G = exp(lG);
warning(w);
end
