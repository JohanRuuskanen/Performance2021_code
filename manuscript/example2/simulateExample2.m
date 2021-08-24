function simulateExample2(suffix, logpath, mu, seed, timespan)
    addpath(genpath('../../line'));
    
    disp("Simulating")

    D = {}; 
    D{1} = Coxian.fitMeanAndSCV(mu, 1.0);
    D{2} = Coxian.fitMeanAndSCV(0.5, 0.5);
    D{3} = Coxian.fitMeanAndSCV(1.0, 10.0);
    
    [model, node, jobclass, P] = genModel(D, logpath, suffix);

    solverjmt = SolverJMT(model, 'seed', seed, 'samples', 1e9, 'timespan', timespan);
    avgJMT = solverjmt.getAvg();
    save_data_from_sim(node, jobclass, P, D, avgJMT, logpath, suffix);

end
    

function [model, node, jobclass, P] = genModel(D, logpath, suffix)
    
    model = Network('model');

    node{1} = Delay(model, ['Delay1_', suffix]); 
    node{2} = Queue(model, ['Queue2_', suffix], SchedStrategy.PS);
    node{3} = Queue(model, ['Queue3_', suffix], SchedStrategy.PS);
   
    jobclass{1} = ClosedClass(model, 'Class1', 50, node{1}, 0);
    
    node{1}.setService(jobclass{1}, D{1});
    
    node{2}.setService(jobclass{1}, D{2});
    node{2}.setNumberOfServers(4);
    
    node{3}.setService(jobclass{1}, D{3});
    node{3}.setNumberOfServers(8);
    
    P = model.initRoutingMatrix;

    P{1,1} = [  
        0.0, 1.0, 0.0;
        0.0, 0.0, 1.0;
        1.0, 0.0, 0.0;
    ];

    model.linkAndLog(P, ones(size(node)), logpath);
end

function save_data_from_sim(node, jobclass, P, D, avgJMT, logpath, suffix)
    K = {};
    Queues = {};
    Disc = {};
    for k = 1:length(node)
        Q = split(node{k}.name, "_");
        Queues{k} = Q{1};
        Disc{k} = char(node{k}.schedStrategy);
        if strcmp(Disc{k}, 'ext')
            K{k} = [1];
        else
            K{k} = [node{k}.numberOfServers];
        end
    end

    Classes = {};
    Q0 = zeros(length(Queues), length(Classes));
    for k = 1:length(jobclass)
        Classes{k} = jobclass{k}.name;

        Q = split(jobclass{k}.reference.name, "_");
        v = strcmp(Q{1}, Queues);
        if strcmp(jobclass{k}.type, 'closed')
            Q0(v, k) = jobclass{k}.population;
        else
            Q0(v, k) = 0.0;
        end
    end

    Pa = cell2mat(P);
    
    E = {};
    for i = 1:length(D)
        E{i, 1} = D{i}.params{1}.paramValue;
        E{i, 2} = D{i}.params{2}.paramValue;
    end

    save([logpath, '/params_', suffix, '.mat'], ...
        'Pa', 'Classes', 'Queues', 'Disc', 'Q0', 'K', 'E', 'avgJMT');
    
end