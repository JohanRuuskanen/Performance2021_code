% This is a sequence structure detector, it can automatically detect basic
% sequence structures in a workflow network. 

% Matrix A contains the link information of the workflow, the first column
% represents the start node, the second column represents the end node, the
% third column represents the transition probability.

% Service is a list containing the service nodes.

% Chain is a cell, each element of chain is a sequence structure in the workflow. 

function [chain] = seqDetector(A,Service)
%% Step1: the goal is to detect sequence structures 'roughly'.

WF_1 = ismember(A,Service); % WF_ is a matrix or a list which describes how the nodes are connected in the workflow.
                            % The elements of matrix WF_1 are either 1 or 0, 1 is service node, 0 is other node.
indx_1=find(WF_1(:,1)==1 & WF_1(:,2)==1); % indx_1 is the index of rows where service nodes are connected each other.

%% Step2: extract indx_1 rows from WF_1 and use these rows to build a new matrix WF_2.

WF_2 = zeros(length(indx_1),2); 
for i = 1:1:length(indx_1)
    WF_2(i,1) = A(indx_1(i),1);
    WF_2(i,2) = A(indx_1(i),2);
end

%% Step3: the goal is to find the number of occurences of each service node, either 1 or 2. 

WF_3 = WF_2(:); % Flatten the matrix WF_2 and get list WF_3, WF_4 is unique(WF_3),then count number can be gained by histc(WF_3,WF_4). 
WF_4 = unique(WF_3);
Count = histc(WF_3,WF_4); 

%% Step4: the goal is to further detect sequence structures accurately.

% In previous step3,the number of occurences of each service node are gained. If
% the node is start or end node of the sequence, then the number of that
% node will be 1. So, the number of sequences is equal to the half number of 1. 
for j = 1:1:length(find(Count==1))/2 
    [m,n] = size(WF_2);
    first = WF_2(1,1); % Initialize first and last, they are the start and end node of the sequence respectively. 
     last = WF_2(1,2);
    chain{j} = [first,last]; % chain{j} is one sequence structure, this is expanded in the search process below.
    indx_2 = [1]; % indx_2 is a list, the element is the order number of each element in chain{j}.
    difference = 100; % difference serves as a indicator. When difference is greater than zero, it means the sequence chain{j} should be continued to expand.
    while difference > 0
        l = length(chain{j});
        for i = 2:1:m
            if WF_2(i,1) == last
                last = WF_2(i,2);
                chain{j} = [chain{j},last];
                indx_2 = [indx_2,i];
            elseif WF_2(i,2) == first
                first = WF_2(i,1);
                chain{j} = [first,chain{j}];
                indx_2 = [indx_2,i];
            else
                chain{j} = chain{j};
            end
        end
        difference = length(chain{j})-l;
    end
    WF_2(indx_2,:) = []; % Delete the rows what we have searched to simplfy the computation.
end
if length(indx_1) == 0 % This means there is no sequence structure in the workflow network.
    chain = {};
end

end