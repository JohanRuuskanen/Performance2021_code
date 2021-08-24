function features = extractFeatures(model)
% This is a helper function for tabulateData. It takes in a
% Network object (QN model in LINE) and returns a vector
% of 15 features as shown below

	qn = model.getStruct;   
    features = zeros(1, 15);
    
    % Station and scheduling information
    features(1) = sum(qn.schedid == SchedStrategy.ID_FCFS); % Num FCFS queues
    features(2) = sum(qn.schedid == SchedStrategy.ID_PS); % Num PS queues
    features(3) = sum(qn.schedid == SchedStrategy.ID_INF); % Num delays
    features(4) = qn.nnodes - qn.nstations; % Num CS nodes
    features(5) = sum(qn.nservers(~isinf(qn.nservers))); % Num queue servers
    
    % Job information
    features(6) = qn.nchains; % Num chains
    features(7) = qn.nclosedjobs; % Number of jobs in the system
    
    % Service process information
    features(8:10) = getServiceDist(model); % Num of each distribution type
    features(11) = mean(qn.rates, 'all', 'omitnan'); % Avg service rate
    features(12) = mean(qn.scv, 'all', 'omitnan'); % Avg SCV
    features(13) = mean(qn.phases, 'all', 'omitnan'); % Avg phases
    
    % Misc
    features(14) = sum(qn.nodetype == NodeType.Queue) == 1; % If only 1 Queue, special for nc.mmint
    features(15) = model.hasProductFormSolution; % Check if model has product form solution
end

function data = getServiceDist(model)
    numexp = 0;
    numhyperexp = 0;
    numerlang = 0;
    
    for i = 1 : model.getNumberOfStations
        for j = 1 : model.getNumberOfClasses
            switch model.stations{i}.serviceProcess{j}.name
               case 'Exponential'
                   numexp = numexp + 1;
               case 'HyperExp'
                   numhyperexp = numhyperexp + 1;
               case 'Erlang'
                   numerlang = numerlang + 1;
            end
        end
    end

    data = [numexp numhyperexp numerlang];
end