function [ph, phases] = refreshMarkovianService(self)
% [PH, PHASES] = REFRESHPHSERVICE()

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

M = self.getNumberOfStations();
K = self.getNumberOfClasses();
ph = cell(M,K);
phases = zeros(M,K);
for i=1:M
    if i == self.getIndexSourceStation
        ph_i = self.stations{i}.getMarkovianSourceRates();
    else
        switch class(self.stations{i})
            case 'Fork'
                mu_i = cell(1,K);
                phi_i = cell(1,K);
                for r=1:K
                    mu_i{r} = NaN;
                    phi_i{r} = NaN;
                end
                ph_i = Coxian(mu_i,phi_i).getRepresentation;
            case 'Join'
                mu_i = cell(1,K);
                phi_i = cell(1,K);
                for r=1:K
                    mu_i{r} = NaN;
                    phi_i{r} = NaN;
                    ph_i{r} = Coxian(mu_i{r},phi_i{r}).getRepresentation;
                end
            otherwise
                ph_i = self.stations{i}.getMarkovianServiceRates();
        end
    end
    for r=1:K
        ph{i,r} = ph_i{r};
        if isempty(ph{i,r}) % fluid fails otherwise
            phases(i,r) = 1;
        elseif any(isnan(ph{i,r}{1}(:))) || any(isnan(ph{i,r}{2}(:))) % disabled
            phases(i,r) = 0;
        else
            phases(i,r) = length(ph{i,r}{1});
        end
    end
end
if ~isempty(self.qn) %&& isprop(self.qn,'mu')
    self.qn.setMAPService(ph, phases);
end
end
