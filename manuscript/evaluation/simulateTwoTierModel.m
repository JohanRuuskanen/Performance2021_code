function simulateTwoTierModel(suffix, logpath, simSettings, seed, timespan)
    addpath(genpath('../../line'));
    
    disp("Simulating")

    nodeParams = getNodeParams(simSettings);

    [model, node, jobclass, P] = genTwoTierModel(nodeParams, simSettings.N, ...
        logpath, suffix);

    solverjmt = SolverJMT(model, 'seed', seed, 'samples', 1e9, ...
        'timespan', timespan);
    avgJMT = solverjmt.getAvg();
    D = getDistsFromNodeParams(nodeParams);
    save_data_from_sim(node, jobclass, P, D, avgJMT, logpath, suffix);

end
    

function [model, node, jobclass, P] = genTwoTierModel(nodeParams, N, logpath, suffix)

    nf = length(nodeParams.frontends);
    nb = length(nodeParams.backends);

    D0 = Exp.fitMean(1e9);
    
    model = Network('model');

    node{1} = Source(model, ['Source_', suffix]); 
    node{2} = Delay(model, ['Delay1_', suffix]);
    for k = 1:nf
        name =  ['Frontend', int2str(k), '_', suffix];
        node{end+1} = Queue(model, name, SchedStrategy.PS);
    end
    for k = 1:nb
        name =  ['Backend', int2str(k), '_', suffix];
        node{end+1} = Queue(model, name, SchedStrategy.PS);
    end
    node{end+1} = Sink(model, ['Sink_', suffix]);

    jobclass{1} = OpenClass(model, 'OpenClass1', 0);
    jobclass{2} = OpenClass(model, 'OpenClass2', 0);
    jobclass{3} = ClosedClass(model, 'ClosedClass3', N, node{2}, 0);
    jobclass{4} = ClosedClass(model, 'ClosedClass4', 0, node{2}, 0);
    
    node{1}.setArrival(jobclass{1}, nodeParams.source{1});

    node{2}.setService(jobclass{3}, nodeParams.delay{1});

    
    for k = 1:nf
        node{2+k}.setService(jobclass{1}, nodeParams.frontends{k}.phdist{1});
        node{2+k}.setService(jobclass{2}, nodeParams.frontends{k}.phdist{2});
        node{2+k}.setService(jobclass{3}, nodeParams.frontends{k}.phdist{3});
        node{2+k}.setService(jobclass{4}, nodeParams.frontends{k}.phdist{4});
        node{2+k}.setNumberOfServers(nodeParams.frontends{k}.K);
    end
    
    for k = 1:nb
        node{2+nf+k}.setService(jobclass{1}, nodeParams.backends{k}.phdist{1});
        node{2+nf+k}.setService(jobclass{3}, nodeParams.backends{k}.phdist{3});
        node{2+nf+k}.setNumberOfServers(nodeParams.backends{k}.K);
    end
    
    P = model.initRoutingMatrix;
    
    for i = 1:length(jobclass)
        for j = 1:length(jobclass)
            P{i,j} = zeros(length(node), length(node));
        end
    end

    P{1, 1}(1, 2+1:2+nf) = 1/nf;
    P{1, 1}(2+1:2+nf, 2+nf+1:2+nf+nb) = 1/nb;

    P{1, 2}(2+nf+1:2+nf+nb, 2+1:2+nf) = 1/nf;
    P{2, 1}(2+1:2+nf, end) = 1.0;
    
    P{3, 3}(2, 2+1:2+nf) = 1/nf;
    P{3, 3}(2+1:2+nf, 2+nf+1:2+nf+nb) = 1/nb;
    P{3, 4}(2+nf+1:2+nf+nb, 2+1:2+nf) = 1/nf;
    P{4, 3}(2+1:2+nf, 2) = 1.0;

    model.linkAndLog(P, [0, ones(1, nf+nb+1), 0], logpath);
end

function D = getDistsFromNodeParams(nodeParams)
    D = {};
    
    D{1} = nodeParams.source{1};
    
    D{2} = nodeParams.delay{1};
    
    for k = 1:length(nodeParams.frontends)
        for i = 1:4
            D{end+1} = nodeParams.frontends{k}.phdist{i};
        end
    end
    
    for k = 1:length(nodeParams.backends)
        D{end+1} = nodeParams.backends{k}.phdist{1};
        D{end+1} = nodeParams.backends{k}.phdist{3};
    end
    
end

function nodeParams = getNodeParams(s)
    nodeParams.source = {Coxian.fitMeanAndSCV(s.ma, s.scva)};
    
    nodeParams.delay = {Coxian.fitMeanAndSCV(s.md, s.scvd)};
    
    for k = 1:s.nf
        i = (k-1)*4;
        nodeParams.frontends{k}.K = s.Kf(k);
        nodeParams.frontends{k}.phdist = ...
            {Coxian.fitMeanAndSCV(s.mf(i+1), s.scvf(i+1)), ...
             Coxian.fitMeanAndSCV(s.mf(i+2), s.scvf(i+2)), ...
             Coxian.fitMeanAndSCV(s.mf(i+3), s.scvf(i+3)), ...
             Coxian.fitMeanAndSCV(s.mf(i+4), s.scvf(i+4))};
    end
    
    for k = 1:s.nb
        i = (k-1)*2;
        nodeParams.backends{k}.K = s.Kb(k);
        nodeParams.backends{k}.phdist = ...
            {Coxian.fitMeanAndSCV(s.mb(i+1), s.scvb(i+1)), ...
             [], ...
             Coxian.fitMeanAndSCV(s.mb(i+2), s.scvb(i+2)), ...
             []};
    end
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