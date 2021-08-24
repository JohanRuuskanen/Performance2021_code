function [QNn,UNn,RNn,TNn,ANn] = getAvgNode(self, Q, U, R, T, A)
% [QNN,UNN,RNN,TNN] = GETNODEAVG(Q, U, R, T, A)

% Compute average utilizations at steady-state for all nodes
if nargin == 1 % no parameter
    if isempty(self.model.handles) || ~isfield(self.model.handles,'Q') || ~isfield(self.model.handles,'U') || ~isfield(self.model.handles,'R') || ~isfield(self.model.handles,'T') || ~isfield(self.model.handles,'A')
        self.resetResults(); % reset in case there are partial results saved
    end
    [Q,U,R,T] = self.model.getAvgHandles;
elseif nargin == 2
    handlers = Q;
    Q=handlers{1};
    U=handlers{2};
    R=handlers{3};
    T=handlers{4};
end
[QN,UN,RN,TN] = self.getAvg(Q,U,R,T);

qn = self.model.getStruct; % must be called after getAvg

I = self.model.getNumberOfNodes;
M = self.model.getNumberOfStations;
R = self.model.getNumberOfClasses;
C = self.model.getNumberOfChains;
QNn = zeros(I,R);
UNn = zeros(I,R);
RNn = zeros(I,R);
TNn = zeros(I,R);
ANn = zeros(I,R);

for ist=1:M
    ind = qn.stationToNode(ist);
    QNn(ind,:) = QN(ist,:);
    UNn(ind,:) = UN(ist,:);
    RNn(ind,:) = RN(ist,:);
    %TNn(ind,:) = TN(ist,:); % this is built later from ANn
    %ANn(ind,:) = TN(ist,:);
end

%if any(qn.isstatedep(:,3)) || any(qn.nodetype == NodeType.Cache)
%    line_warning(mfilename,'Node-level metrics not available in models with state-dependent routing. Returning station-level metrics only.');
%    return
%end

% update tputs for all nodes but the sink and the joins
for ind=1:I
    for c = 1:C
        inchain = find(qn.chains(c,:));
        for r = inchain
            %anystat = find(qn.visits{c}(:,r));
            refstat = qn.refstat(c);
            %if ~isempty(anystat)
            if qn.nodetype(ind) ~= NodeType.Source
                %if qn.isstation(ind)
                %    ist = qn.nodeToStation(ind);
                %    ANn(ind, r) =  TN(ist,r);
                %else
                ANn(ind, r) =  (qn.nodevisits{c}(ind,r) / sum(qn.visits{c}(refstat,inchain))) * TN(refstat,r);
                %end
            end
            %end
        end
    end
end

for ind=1:I
    for c = 1:C
        inchain = find(qn.chains(c,:));
        for r = inchain
            anystat = find(qn.visits{c}(:,r));
            if ~isempty(anystat)
                if qn.nodetype(ind) ~= NodeType.Sink && qn.nodetype(ind) ~= NodeType.Join
                    for s = inchain
                        for jnd=1:I
                            if qn.nodetype(ind) ~= NodeType.Source
                                TNn(ind, s) = TNn(ind, s) + ANn(ind, r) * qn.rtnodes((ind-1)*R+r, (jnd-1)*R+s);
                            else
                                ist = qn.nodeToStation(ind);
                                TNn(ind, s) = TN(ist,s);
                            end
                        end
                    end
                end
            end
        end
    end
end

ANn(isnan(ANn)) = 0;
TNn(isnan(TNn)) = 0;

end
