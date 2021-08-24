% This is a parallel structure detector, it can automatically detect basic
% parallel structures in a workflow network. 

% Matrix A contains the link information of the workflow, the first column
% represents the start node, the second column represents the end node, the
% third column represents the transition probability.

% Service is a list containing the service nodes.

% Fork is a list containing the Fork nodes.

% Join is a list containing the Join nodes.

% Chain is a cell, each element of chain is a parallel structure in the workflow. 

function [chain] = paraDetector(A,Service,Fork,Join)
%% Step1: the goal is to detect parallel structures 'roughly'. This goal is achieved by detecting the Fork nodes.

WF_1 = ismember(A,Fork); % WF_ is a matrix or a list which describes how the nodes are connected in the workflow.
                         % The elements of matrix WF_1 are either 1 or 0, 1 is Fork node, 0 is other node.
indx_1 = find(WF_1(:,1)==1); % indx_1 is the index of rows which start with Fork node.

%% Step2: extract indx_1 rows from WF_1 and use these rows to build a new matrix WF_2.

WF_2 = zeros(length(indx_1),2);
for i = 1:1:length(indx_1)
    WF_2(i,1) = A(indx_1(i),1);
    WF_2(i,2) = A(indx_1(i),2);
end

%% Step3: Recognize the fork nodes which connect service nodes.

WF_3 = ismember(WF_2,Service); % The elements of matrix WF_3 are either 1 or 0, 1 is service node, 0 is other node.
indx_2 = find(WF_3(:,2)==0); % indx_2 is the index of rows where Fork nodes connect non-service nodes.
WF_4 = zeros(1,length(indx_2));
for i = 1:1:length(indx_2)
    WF_4(i) = WF_2(indx_2(i),1);
end
% WF_4 is a list of Fork nodes which connect non-service nodes. WF_5 is gained
% by deleting WF_4 from Fork list, so all the fork nodes of WF_5 connect
% service nodes.
WF_5 = setxor(Fork,WF_4);

%% Step4: based on step3, extract rows (which start with nodes of WF_5) from original matrix A and use these rows to build a new matrix WF_7.

WF_6 = ismember(A,WF_5);
indx_3 = find(WF_6(:,1)==1);
WF_7 = zeros(length(indx_3),2);
for i = 1:1:length(indx_3)
    WF_7(i,1) = A(indx_3(i),1);
    WF_7(i,2) = A(indx_3(i),2);
end

%% Step5: the goal is to further detect parallel structures accurately.

for k = 1:1:length(WF_5)
    % Initialize chain{k} and LIST.chain{k} will be the list containing parallel service nodes.
    % LIST will be the list whose elements are connected by parallel service nodes.
    chain{k} = [];
    LIST = [];
    WF_8 = ismember(WF_7,WF_5(k)); 
    indx_4 = find(WF_8(:,1)==1); % indx_4 is the index of rows which start with a certain Fork node of WF_5.
    for i = 1:1:length(indx_4)
        chain{k} = [chain{k},WF_7(indx_4(i),2)]; % chain{k} is the list containing parallel service nodes.
        WF_9 = ismember(A,WF_7(indx_4(i),2));
        indx_5 = find(WF_9(:,1)==1);
        LIST = [LIST,A(indx_5,2)]; % LIST is the list whose elements are connected by parallel service nodes.
    end
    % Judge if the chain{k} is the parallel structure. If not, chain{k} will be deleted from the cell.
    if length(unique(LIST)) == 1 && find(Join == unique(LIST))>0
        chain{k} = chain{k};
    else
        chain{k} = [];
    end
end
if length(WF_5)>0
    id = cellfun('length',chain);
    chain(id==0) = [];
else
    chain = {}; % This means there is no parallel structure in the workflow network.
end

end