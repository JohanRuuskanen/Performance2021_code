function [alpha_final,T_final,view] = calcPH(A,Service,Fork,Join,Router,dataset)

% This function is responsible for calculating the PH parameters (alpha and
% T) of a workflow network.

% Matrix A contains the link information of the workflow, the first column
% represents the start node, the second column represents the end node, the
% third column represents the transition probability.

% Service is a list containing the service nodes.

% Fork is a list containing the Fork nodes.

% Join is a list containing the Join nodes.

% Router is a list containing the Router nodes.

% dataset is a cell. Each element is a list of data.

% The output arguments are alpha_final, T_final and view, alpha_final and
% T_final are the PH parameters of the workflow network. view is a cell,
% it serves as a descriptor which can automatically describe how a workflow will be simplified.

%% Step1: initialize the param, view and calc the number of basic units in the workflow.

param = cell(length(dataset),2);
for i = 1:1:length(dataset)
    [m1,m2,m3] = getMoments(dataset{i});
    [alpha,T] = getAPH2(m1,m2,m3);
    param{i,1} = alpha;
    param{i,2} = T;
end

count = 1;
view{count} = A;

basic_s = seqDetector(A,Service);
basic_p = paraDetector(A,Service,Fork,Join);
basic_l = loopDetector(A,Service,Router);
basic_b = branchDetector(A,Service,Join);
Indicator = length(basic_s)+length(basic_p)+length(basic_l)+length(basic_b);

%% Step2: If there are basic units in the workflow network, then it will be simplified and updated by the function patternUpdate. 
          % This process will be implemented
          % repeatedly until there is no basic units in the workflow. And
          % then calc the PH parameters of this workflow.

while Indicator>0
    [new_a,new_param] = patternUpdate(A,Service,Fork,Join,Router,param);
    A = new_a;
    param = new_param;
    
    count = count+1;
    view{count} = A;
    
    basic_s = seqDetector(A,Service);
    basic_p = paraDetector(A,Service,Fork,Join);
    basic_l = loopDetector(A,Service,Router);
    basic_b = branchDetector(A,Service,Join);
    Indicator = length(basic_s)+length(basic_p)+length(basic_l)+length(basic_b);
end

% Find the final simplified service node and calc its PH parameters.
if find(Service == A(1,2))>0
    indx_final = find(Service == A(1,2));
else
    indx_final = find(Service == A(1,1));
end
C_alpha = param{indx_final,1};
C_T = param{indx_final,2};
alpha_final = real(C_alpha);
T_final = real(C_T);
end