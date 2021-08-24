% This is a branch structure detector, it can automatically detect basic
% loop structures in a workflow network. 

% Matrix A contains the link information of the workflow, the first column
% represents the start node, the second column represents the end node, the
% third column represents the transition probability.

% Service is a list containing the service nodes.

% Router is a list containing the Router nodes.

% Chain is a cell, each element of chain is a branch structure in the workflow. 
function [chain] = branchDetector(A,Service,Join)
%% Step1: the goal is to detect branch structures 'roughly'. This goal is achieved by detecting the rows where the transition probabilities are less than one.

indx_1= find(A(:,3)<1); % indx_1 is the index of rows where the transition probabilities are less than one.
WF_1 = [];
for i = 1:1:length(indx_1)
    WF_1 = [WF_1,A(indx_1(i),1)];
end
%% Step2: the goal is to further detect branch structures accurately.

WF_2 = unique(WF_1); % The elements of WF_2 list are the Router nodes having transition probabilities less than one.
for j = 1:1:length(WF_2)
    WF_3 = ismember(A,WF_2(j));
    indx_2 = find(WF_3(:,1) == 1);
    LIST = [];
    % LIST is the list whose elements are connected by Router nodes of WF_2 list.
    for k = 1:1:length(indx_2)
        LIST = [LIST,A(indx_2(k),2)];
    end
    % Judge whether all of the elements in LIST are service nodes. If they are, construct LIST_1.
    if isequal(intersect(LIST,Service),LIST) == 1
        LIST_1 = [];
        % LIST_1 is the list whose elements are the nodes after those of list.
        for p = 1:1:length(LIST)
            WF_4 = ismember(A,LIST(p));
            indx_3 = find(WF_4(:,1) == 1);
            LIST_1 = [LIST_1,A(indx_3,2)];
        end
        % Judge whether all the elements of LIST_1 are the same Router node.
        % If they are, then chain{j} is branch structure.
        if length(unique(LIST_1)) == 1 & find(Join == unique(LIST_1))>0
            chain{j} = LIST;
        else
            chain{j} = [];
        end
    else
        chain{j} = [];
    end
end
if length(indx_1)>0
    id = cellfun('length',chain);
    chain(id==0) = [];
else
    chain = {}; % This means there is no branch structure in the workflow network.
end

end