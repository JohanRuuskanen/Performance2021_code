function space = spaceLocalVars(qn, ind)
% SPACE = SPACELOCALVARS(QN, IND)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

% Generate state space for local state variables

%ind: node index
%ist = qn.nodeToStation(ind);
%isf = qn.nodeToStateful(ind);

space = [];

switch qn.nodetype(ind)
    case NodeType.Cache
        space = State.spaceCache(qn.varsparam{ind}.nitems,qn.varsparam{ind}.cap);
end

switch qn.routing(ind)
    case RoutingStrategy.ID_RRB
        space = State.decorate(space, qn.varsparam{ind}.outlinks(:));
end
end
