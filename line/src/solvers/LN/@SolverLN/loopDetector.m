% This is a loop structure detector, it can automatically detect basic
% loop structures in a workflow network. 

% Matrix A contains the link information of the workflow, the first column
% represents the start node, the second column represents the end node, the
% third column represents the transition probability.

% Service is a list containing the service nodes.

% Router is a list containing the Router nodes.

% Chain is a cell, each element of chain is a loop structure in the workflow. 
function [chain] = loopDetector(A,Service,Router)
%% Step1: the goal is to detect loop structures 'roughly'. This goal is achieved by detecting the rows where Router nodes are connected with each other.

WF_1 = ismember(A,Router); % WF_ is a matrix or a list which describes how the nodes are connected in the workflow. 
                           % The elements of matrix WF_1 are either 1 or 0, 1 is Router node, 0 is other node.
indx_1 = find(WF_1(:,1)==1 & WF_1(:,2)==1); % indx_1 is the index of rows where Router nodes are connected with each other.

%% Step2: the goal is to further detect loop structures accurately.

for i = 1:1:length(indx_1)
    WF_2 = ismember(A,A(indx_1(i),1));
    indx_2 = find(WF_2(:,2)==1); % indx_2 is the index of rows where a node (denote NODE1) connects to a Router node.
    WF_3 = ismember(A,A(indx_1(i),2));
    indx_3 = find(WF_3(:,1)==1); % indx_3 is the index of rows where a Router node connects a node (denote NODE2).
    % Judge whether NODE1 and NODE2 are the same service node. If they are,
    % then chain{i} is loop structure.
    if A(indx_2,1) == A(indx_3,2) & find(Service == A(indx_2,1))>0
        chain{i} = [A(indx_2,1)];
    else
        chain{i} = [];
    end
end
if length(indx_1)>0
    id = cellfun('length',chain);
    chain(id==0) = [];
else
    chain = {}; % This means there is no loop structure in the workflow network.
end

end