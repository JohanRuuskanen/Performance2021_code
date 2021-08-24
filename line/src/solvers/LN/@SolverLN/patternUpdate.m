% This function is responsible for simplifying and updating workflow
% network after recognizing the basic workflow patterns by the functions
% seqDetector, paraDetector, loopDetector and branchDetector.

% Matrix A contains the link information of the workflow, the first column
% represents the start node, the second column represents the end node, the
% third column represents the transition probability.

% Service is a list containing the service nodes.

% Fork is a list containing the Fork nodes.

% Join is a list containing the Join nodes.

% Router is a list containing the Router nodes.

% param is a cell, each row of the cell represents the PH parameters (alpha
% and T) of a service node. 

% The output arguments are A and param, matrix A represents the updated
% workflow network, param represents the PH parameters of service nodes of
% the updated workflow.
function [A,param] = patternUpdate(A,Service,Fork,Join,Router,param)
%% Step1: simplify the sequence patterns and update the workflow network

s = seqDetector(A,Service);
p = paraDetector(A,Service,Fork,Join);
l = loopDetector(A,Service,Router);
b = branchDetector(A,Service,Join);

if length(s)>0
    WF_S1 = ismember(A,Service); % The elements of matrix WF_S1 are either 1 or 0, 1 is service node, 0 is other node.
    indx_S1=find(WF_S1(:,1)==1 & WF_S1(:,2)==1); % indx_S1 is the index of rows where service nodes are connected each other.
    % Update the matrix A by substituting new rows for the rows of sequence
    % patterns. This means the workflow network is updated by simplifying
    % the sequence patterns. 
    A(indx_S1,:) = [];   
    for i = 1:1:length(s)
        last = length(s{i});
        A(A == s{i}(last)) = s{i}(1);
    end
    
    % So far, A has been updated. Next, param will be updated.

    for i = 1:1:length(s)
        s_list = cell(1,2*length(s{i}));
        % The goal of this for loop is to construct a list containing the
        % PH parameters of each service of sequence patterns. This list  
        % will be fed to the function AS to calc the PH parameters of
        % the sequence patterns.
        for j = 1:1:length(s{i})
            indx_S2 = find(Service == s{i}(j)); % indx_S2 is the index of service node s{i}(j) in the Service list.
            s_list{1,2*j-1} = param{indx_S2,1};
            s_list{1,2*j} = param{indx_S2,2};
        end
        indx_S3 = find(Service == s{i}(1)); % indx_S3 is the index of service node s{i}(1) in the Service list.
        [param{indx_S3,1},param{indx_S3,2}] = convAS(s_list);
    end
    % param has been updated.
end

%% Step2: simplify the parallel patterns and update the workflow network

if length(p)>0
    for i = 1:1:length(p)
        WF_P1 = ismember(A,p{i}(1));
        indx_P1 = find(WF_P1(:,1) == 1);
        jnode = A(indx_P1,2);
        indx_P2 = find(WF_P1(:,2) == 1);
        fnode = A(indx_P2,1);
        % Update the matrix A by substituting new rows for the rows of parallel
        % patterns. This means the workflow network is updated by simplifying
        % the parallel patterns.
        for j = 1:1:length(p{i})
            WF_P2 = ismember(A,p{i}(j));
            indx_P3=find(WF_P2(:,1) == 1 | WF_P2(:,2) == 1);
            A(indx_P3,:) = [];
        end
        A(A == jnode | A == fnode) = p{i}(1);
    end
    
    % So far, A has been updated. Next, param will be updated.
    
    for i = 1:1:length(p)
        p_list = cell(1,2*length(p{i}));
        % The goal of this for loop is to construct a list containing the
        % PH parameters of each service of parallel patterns. This list  
        % will be fed to the function AS to calc the PH parameters of
        % the parallel patterns.
        for j = 1:1:length(p{i})
            indx_P4 = find(Service == p{i}(j)); % indx_P4 is the index of service node p{i}(j) in the Service list.
            p_list{1,2*j-1} = param{indx_P4,1};
            p_list{1,2*j} = param{indx_P4,2};
        end
        indx_P5 = find(Service == p{i}(1)); % indx_P5 is the index of service node p{i}(1) in the Service list.
        [param{indx_P5,1},param{indx_P5,2}] = convAP(p_list);
    end
    % param has been updated.
end

%% Step3: simplify the loop patterns and update the workflow network

if length(l)>0
    for i = 1:1:length(l)
        WF_L1 = ismember(A,l{i});
        indx_L1 = find(WF_L1(:,1) == 1);
        nextRouter = A(indx_L1,2);
        indx_L2 = find(WF_L1(:,2) == 1);
        previousRouter = A(indx_L2,1);
        % Update the matrix A by substituting new rows for the rows of loop
        % patterns. This means the workflow network is updated by simplifying
        % the loop patterns.
        indx_L3=find(WF_L1(:,1) == 1 | WF_L1(:,2) == 1);
        A(indx_L3,:) = [];
        indx_L4 = find(A(:,1) == nextRouter & A(:,2) == previousRouter);
        A(indx_L4,:) = [];
        A(A == nextRouter | A == previousRouter) = l{i};
    end
    
    % So far, A has been updated. Next, param will be updated.
    
    for i = 1:1:length(l)
        indx_L5 = find(Service == l{i});
        alpha_L = param{indx_L5,1};
        T_L = param{indx_L5,2};
        
        WF_L2 = ismember(A,l{i});
        indx_L6 = find(WF_L2(:,1) == 1);
        p = A(indx_L6,3);
        A(indx_L6,3) = 1; % change the transition probability after simplifying the loop patterns.
        
        [param{indx_L5,1},param{indx_L5,2}] = Simplify(alpha_L,T_L,0,0,1-p,0,4);
    end
    % param has been updated.
end

%% Step4: simplify the branch patterns and update the workflow network

if length(b)>0
    for i = 1:1:length(b)
        WF_B1 = ismember(A,b{i}(1));
        indx_B1 = find(WF_B1(:,1) == 1);
        jRouter = A(indx_B1,2);
        indx_B2 = find(WF_B1(:,2) == 1);
        fRouter = A(indx_B2,1);
        prob = [];
        % Update the matrix A by substituting new rows for the rows of
        % branch patterns. This means the workflow network is updated 
        % by simplifying the branch patterns.
        for j = 1:1:length(b{i})
            WF_B2 = ismember(A,b{i}(j));
            indx_B3 = find(WF_B2(:,2) == 1);
            prob = [prob,A(indx_B3,3)];
            indx_B4 = find(WF_B2(:,1) == 1 | WF_B2(:,2) == 1);
            A(indx_B4,:) = [];
 
        end
        A(A == jRouter | A == fRouter) = b{i}(1);
        probcell{i} = prob;
    end
   
    % So far, A has been updated. Next, param will be updated.
    
    for i = 1:1:length(b)
        indx_B5 = find(Service == b{i}(1));
        alpha_B1 = param{indx_B5,1};
        T_B1 = param{indx_B5,2};
        p_B1 = probcell{i}(1);
        
        indx_B6 = find(Service == b{i}(2));
        alpha_B2 = param{indx_B6,1};
        T_B2 = param{indx_B6,2};
        probcell{i};
        p_B2 = probcell{i}(2);
        
        [alpha_B3,T_B3] = Simplify(alpha_B1,T_B1,alpha_B2,T_B2,p_B1,p_B2,3);
        [param{indx_B5,1},param{indx_B5,2}] = reduceOrder(alpha_B3,T_B3);
    end
    % param has been updated.
end

end