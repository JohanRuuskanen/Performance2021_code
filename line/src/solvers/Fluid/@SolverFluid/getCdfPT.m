function CDc = getCdfPT(self)
% CDC = GETCDFPT()

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

T0 = tic;
if ~exist('R','var')
    R = self.model.getAvgRespTHandles;
end
qn = self.getStruct;
if ~exist('withRefStat','var')
    withRefStat = false(1,qn.nclasses);
elseif numel(withRefStat) == 1
    val = withRefStat;
    withRefStat = false(1,qn.nclasses);
    withRefStat(:) = val;
end
% ptSpec = struct(); % passage time specification
% ptSpec.starts = false(qn.nnodes,qn.nclasses,qn.nnodes,qn.nclasses);
% ptSpec.completes = false(qn.nnodes,qn.nclasses,qn.nnodes,qn.nclasses);
% for r=1:qn.nclasses
%     if withRefStat(r)
%         % starts when arriving to ref
%         ptSpec.starts(:,:,qn.refstat(r),r) = true;
%     else % ref station excluded
%         % starts when leaving ref
%         ptSpec.starts(qn.refstat(r),r,:,:) = true;
%     end
%     % completes when arriving to ref
%     if R{qn.refstat(r),r}.class.completes
%         % class switch to r is right after departure from station i
%         ptSpec.completes(:,:,qn.refstat(r),r) = true;
%     end
% end
% options = self.getOptions;
% options.psgtime = ptSpec;
completes = false(qn.nnodes,qn.nclasses);
for i=1:qn.nstations
    for r=1:qn.nclasses
        if R{i,r}.class.completes
            completes(i,r) = true;
        end
    end
end
CDc = solver_fluid_RT(qn, self.result.solverSpecific.odeStateVec, options, completes);
runtime = toc(T0);
self.setDistribResults(CDc, runtime);
end
