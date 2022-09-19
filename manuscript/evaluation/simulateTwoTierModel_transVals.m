function simulateTwoTierModel(logpath, suffix, simSettings, itrs, timespan)
    addpath(genpath('../../line'));
    
    disp("Simulating")

    nodeParams = getNodeParams(simSettings);

    [model, node, jobclass, P] = genTwoTierModel(nodeParams, simSettings.N, ...
        logpath, suffix);

    t = {};
    q = {};
    
    for k = 1:itrs
        seed = 1 + 123*k;
        solver = SolverJMT(model, 'seed', seed, 'timespan', timespan);
        res = solver.sampleSysAggr();
        t{k} = res.t;
        q{k} = cell2mat({res.state{2:end}});

        if any(q{k} < 0)
            t{k} = [];
            q{k} = [];
        end

         qnan = isnan(q{k});
         if sum(sum(qnan) > 0)
             [cols, ~] = find(qnan);
             t{k}(cols) = [];
             q{k}(cols, :) = [];
             disp(["Found NaN values, removing ", int2str(length(cols)), " rows"])
         end

    end

    
    save([logpath, '/transientVals_', suffix, '.mat'], 't', 'q');

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

    model.link(P);
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