function addJobClass(self, customerClass)
% ADDJOBCLASS(CUSTOMERCLASS)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
if self.doChecks
    if sum(cellfun(@(x) strcmp(x.name,customerClass.name), {self.classes{1:end}}))>0
        line_error(mfilename,'A class with name %s already exists.\n', customerClass.name);
    end
end
nClasses = length(self.classes);
customerClass.index = nClasses+1;
self.classes{end+1,1} = customerClass;
self.setUsedFeatures(class(customerClass)); % open or closed
end
