function [ph, mu, phi, phases] = refreshServicePhases(self)
% [PH, MU, PHI, PHASES] = REFRESHSERVICEPHASES()
% Obtain information about phases of service and arrival processes.

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

M = self.getNumberOfStations();
K = self.getNumberOfClasses();

mu = cell(M,K);
phi = cell(M,K);
phases = zeros(M,K);
for i=1:M
    if i == self.getIndexSourceStation
        [~,mu_i,phi_i] = self.stations{i}.getMarkovianSourceRates();
    else
        switch class(self.stations{i})
            case 'Fork'
                mu_i = cell(1,K);
                phi_i = cell(1,K);
                for r=1:K
                    mu_i{r} = NaN;
                    phi_i{r} = NaN;
                end
            case 'Join'
                mu_i = cell(1,K);
                phi_i = cell(1,K);
                for r=1:K
                    mu_i{r} = NaN;
                    phi_i{r} = NaN;
                end
            otherwise
                [~,mu_i,phi_i] = self.stations{i}.getMarkovianServiceRates();
        end
    end
    for r=1:K
        mu{i,r} = mu_i{r};
        phi{i,r} = phi_i{r};
        if isnan(mu_i{r}) % disabled
            phases(i,r) = 0;
        else
            phases(i,r) = length(mu_i{r});
        end
    end
end

if ~isempty(self.qn) %&& isprop(self.qn,'mu')
    self.qn.setCoxService(mu, phi, phases);
end
[ph, phases] = refreshMarkovianService(self);
end