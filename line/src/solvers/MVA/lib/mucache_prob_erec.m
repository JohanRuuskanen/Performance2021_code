function prob = mucache_prob_erec(gamma,m)
[n,h]=size(gamma);
E = mucache_erec(gamma, m);
for i=1:n
    for j=1:h
        Ei = mucache_erec(gamma(setdiff(1:n,i),:),oner(m,j));
        prob(i,1+j) = m(j) * gamma(i,j) * Ei / E;
    end
    prob(i,1) = abs(1 - sum(prob(i,2:end)));
end
end