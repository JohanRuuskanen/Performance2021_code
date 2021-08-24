function classprio = refreshPriorities(self)
% CLASSPRIO = REFRESHPRIORITIES()

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

K = self.getNumberOfClasses();
classprio = zeros(1,K);
for r=1:K
    classprio(r) = self.getClassByIndex(r).priority;
end
if ~isempty(self.qn)
    self.qn.setPrio(classprio);
end
end
