function [lt] = refreshLST(self)
% [LT] = REFRESHLAPLST()
% Refresh the Laplace-Stieltjes transforms in the NetworkStruct object

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

M = self.getNumberOfStations();
K = self.getNumberOfClasses();
lt = cell(M,K);
for i=1:M
    for r=1:K
        if i == self.getIndexSourceStation
            if  isa(self.stations{i}.input.sourceClasses{r}{end},'Disabled')
                lt{i,r} = [];
            else
                lt{i,r} = @(s) self.stations{i}.arrivalProcess{r}.evalLST(s);
            end
        else
            switch class(self.stations{i})
                case {'Fork'}
                    lt{i,r} = [];
                case {'Join'}
                    lt{i,r} = [];
                otherwise
                    lt{i,r} = @(s) self.stations{i}.serviceProcess{r}.evalLST(s);
            end
        end
    end
end
if ~isempty(self.qn) %&& isprop(self.qn,'mu')
    self.qn.setLSTs(lt);
end
end
