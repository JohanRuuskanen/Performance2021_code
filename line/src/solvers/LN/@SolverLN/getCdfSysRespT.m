% This function is responsible for calculating the response time of every
% entry in a Layered Queueing Network (LQN) model.

% The input arguments are model and lqn and parent. The parameter model represents a
% certain submodel of LQN. The parameter lqn is the whole LQN model. In
% particular, lqn = lqnmodel.getStruct. parent is a list given by getStruct
% method

% The output argument is a N*4 cell. N is the number of entries in the 
% submodel. The first column represents the entry ID (the ID in lqn.graph).  
% The second and thethird column represents alpha and T respectively, they
%  are the APH parameters. The fourth column is the entry response time.

function [Output,parent] = getCdfSysRespT(model,lqn,parent)
%function Output = getCdfSysRespT(model,lqn)

%% Step1: change sparse matrix lqn.graph into matrix A
sparse_m = lqn.graph; % sparse matrix
transition_m = full(sparse_m); % transition matrix
[row,col]=find(transition_m~=0);
LENGTH = length(row);
A = zeros(LENGTH,3);
for i = 1:1:LENGTH
    A(i,1) = row(i);
    A(i,2) = col(i);
    A(i,3) = transition_m(row(i),col(i));
end
%% Step2: find all entries and activities of the submodel
Service = [];
Entry = [];
F = SolverFluid(model).getCdfRespT;
indx_data = [];
for j = 1:1:length(model.classes)
    if model.classes{j}.attribute(1) == 3
        Service = [Service,model.classes{j}.attribute(2)];
        indx_data = [indx_data,j];
    end
    if model.classes{j}.attribute(1) == 2
        Entry = [Entry,model.classes{j}.attribute(2)];
    end
end

%% Step3: classify activities according to the entries they belong to
Activity = cell(1,length(Entry));
for i = 1:1:length(Service)
    up = Service(i);
    indicator = 0;
    while indicator == 0
        up = find(lqn.graph(:,up));
        if length(up) == 1
            up = up;
        else
            up = up(1);
        end
        if ismember(up,Entry) == 1
            indicator = 1;
            indx_entry = find(Entry==up);
            Activity{indx_entry} = [Activity{indx_entry},Service(i)];
            parent(Service(i)) = up;
        end
    end
end

Output = cell(length(Entry),4);
for id = 1:1:length(Entry)
%% Step4: add fork and join nodes between activities
    % extract the activities of a certain entry
    certain_Entry = Activity{id};
    WF_1 = ismember(A,certain_Entry);
    indx_1=find(WF_1(:,1)==1 & WF_1(:,2)==1);
    if length(indx_1)>0
        WF_2 = zeros(length(indx_1),3); 
        for i = 1:1:length(indx_1)
            WF_2(i,1) = A(indx_1(i),1);
            WF_2(i,2) = A(indx_1(i),2);
            WF_2(i,3) = A(indx_1(i),3);
        end

        Fork = [];
        Join = [];
        Router = [];
        idx = max(certain_Entry);
        for i = 1:1:length(certain_Entry)
            N1=numel(find(WF_2(:,1) == certain_Entry(i)));
            if N1>1 
                indx_fork = find(WF_2(:,1) == certain_Entry(i));
                fork = [];
                fork_location = [];
                for f = 1:1:length(indx_fork)
                    fork = [fork,WF_2(indx_fork(f),3)]; % the elements of fork list are the branch probabilities
                    fork_location = [fork_location,WF_2(indx_fork(f),2)];
                end
                if length(find(fork==1)) == length(fork)
                    % recognize AND-Fork, add a fork node between activities
                    idx = idx+1;
                    Fork = [Fork,idx];

                    % update the rows of transition matrix A
                    indx_2 = find(WF_2(:,1) == certain_Entry(i));
                    WF_2(indx_2,:) = [];

                    % add new rows (containing fork node) into matrix A
                    new_r = [certain_Entry(i),idx,1];
                    WF_2 = [WF_2;new_r];
                    act = find(lqn.graph(certain_Entry(i),:));
                    for k = 1:1:length(act)
                        new_r = [idx,act(k),1];
                        WF_2 = [WF_2;new_r];
                    end
                end
                if max(fork)<1
                    % recognize OR-Fork, add a router node between activities
                    idx = idx+1;
                    Router = [Router,idx];

                    % update the rows of transition matrix A
                    indx_2 = find(WF_2(:,1) == certain_Entry(i));
                    WF_2(indx_2,:) = [];

                    % add new rows (containing router node) into matrix A
                    new_r = [certain_Entry(i),idx,1];
                    WF_2 = [WF_2;new_r];
                    act = find(lqn.graph(certain_Entry(i),:));
                    for k = 1:1:length(act)
                        or_prob = fork(find(fork_location == act(k)));
                        new_r = [idx,act(k),or_prob];
                        WF_2 = [WF_2;new_r];
                    end
                end
            end

            N2=numel(find(WF_2(:,2) == certain_Entry(i)));
            if N2>1 
                % recognize Join location, add a join node between activities
                idx = idx+1;
                Join = [Join,idx];

                % update the rows of transition matrix A
                indx_2 = find(WF_2(:,2) == certain_Entry(i));
                WF_2(indx_2,:) = [];

                % add new rows (containing join node) into matrix A
                new_r = [idx,certain_Entry(i),1];
                WF_2 = [WF_2;new_r];
                act = find(lqn.graph(:,certain_Entry(i)));
                for k = 1:1:length(act)
                    new_r = [act(k),idx,1];
                    WF_2 = [WF_2;new_r];
                end
            end
        end

        % add a start and an end node
        for i = 1:1:length(certain_Entry)
            N1=numel(find(WF_2(:,1) == certain_Entry(i)));
            N2=numel(find(WF_2(:,2) == certain_Entry(i)));
            if N1 == 1 && N2 == 0
                idx = idx+1;
                new_r = [idx,certain_Entry(i),1];
                WF_2 = [WF_2;new_r];
            end
            if N1 == 0 && N2 == 1
                idx = idx+1;
                new_r = [certain_Entry(i),idx,1];
                WF_2 = [WF_2;new_r];
            end
        end
    else
        WF_2 = [certain_Entry+1,certain_Entry,1;certain_Entry,certain_Entry+2,1];
        Fork = [];
        Join = [];
        Router = [];
    end
    
%% Step5: matach data (CDF) to its activity
    order = 1;
    data = cell(1,length(certain_Entry));
    for i = 1:1:length(certain_Entry)
        data_location = find(Service==certain_Entry(i));
        data{order} = F{2,indx_data(data_location)};
        order = order+1;       
    end
    
%% Step6: call calcPH function to get APH for entries
    [alpha,T,view] = calcPH(WF_2,certain_Entry,Fork,Join,Router,data);
    SysRespT = alpha*((-T)^(-1))*ones(2,1);
    
    Output{id,1} = Entry(id);
    Output{id,2} = alpha;
    Output{id,3} = T;
    Output{id,4} = SysRespT;
end
end