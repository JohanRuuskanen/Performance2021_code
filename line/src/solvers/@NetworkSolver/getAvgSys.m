function [CNchain,XNchain] = getAvgSys(self,R,T)
% [CNCHAIN,XNCHAIN] = GETAVGSYS(SELF,R,T)

% Return average system metrics at steady state
%
% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

qn = self.model.getStruct();
if nargin < 3
    R = self.model.getAvgRespTHandles;
    T = self.model.getAvgTputHandles;
end
[~,~,RN,TN] = self.getAvg([],[],R,T);

refstats = qn.refstat;
completes = true(1,qn.nclasses);
for r=1:qn.nclasses
    completes(r) = T{refstats(r),r}.class.completes;
end

%if any(isinf(qn.njobs')) % if the model has any open class
% TODO: this could be optimised by computing the statistics
% only for open chains

% compute chain visits
alpha = zeros(qn.nstations,qn.nclasses);
CNclass = zeros(1,qn.nclasses);
for c=1:qn.nchains
    inchain = find(qn.chains(c,:));
    for r=inchain
        CNclass(r)=0;
        for i=1:qn.nstations
            if ~isempty(RN) && ~(isinf(qn.njobs(r)) && i==qn.refstat(r)) % not empty and not source
                CNclass(r) = CNclass(r) + qn.visits{c}(i,r)*RN(i,r)/qn.visits{c}(qn.refstat(r),r);
            end
        end
    end
end

for c=1:qn.nchains
    inchain = find(qn.chains(c,:));
    completingclasses = qn.chains(c,:) & completes;
    for i=1:qn.nstations
        for k=inchain % for all classes within the chain (a class belongs to a single chain, the reference station must be identical for all classes within a chain )
            alpha(i,k) = alpha(i,k) + qn.visits{c}(i,k)/sum(qn.visits{c}(qn.refstat(k),completingclasses));
        end
    end
end
alpha(~isfinite(alpha))=0;
%end

% compute average chain metrics
CNchain = zeros(1,qn.nchains);
XNchain = zeros(1,qn.nchains);
for c=1:qn.nchains
    inchain = find(qn.chains(c,:));
    completingclasses = find(qn.chains(c,:) & completes);
    if ~isempty(TN)
        XNchain(c) = 0;
        % all classes in same chain must share the same refstation, so we use the first one
        ref = refstats(inchain(1));
        % we now compute the incoming system throughput to the
        % reference station from completing classes
        for i=1:qn.nstations
            for r=completingclasses(:)'
                for s=inchain(:)'
                    if ~isnan(TN(i,r))
                        XNchain(c) = XNchain(c) + qn.rt((i-1)*qn.nclasses + r, (ref-1)*qn.nclasses + s )*TN(i,r);
                    end
                end
            end
        end
    end
    
    % If this is a closed chain we simply apply Little's law
    nJobsChain = sum(qn.njobs(find(qn.chains(c,:)))); %#ok<FNDSB>
    %if ~isinf(nJobsChain)
    %    CNchain(c) = nJobsChain / XNchain(c);
    %else % if this is an open chain
    if isinf(nJobsChain)
        if length(inchain) ~= length(completingclasses)
            line_error(mfilename,'Edge-based chain definition not yet supported for open queueing networks.');
            %else
            % we use nan sum to disregard response at stations where
            % the class is not defined
            %    CNchain(c) = sumfinite(alpha(refstats(inchain(1)),inchain).*CNclass(inchain));
        end
    end
    CNchain(c) = sumfinite(alpha(refstats(inchain(1)),inchain).*CNclass(inchain));
end
end
