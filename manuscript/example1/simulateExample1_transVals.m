function simulateExample1_transVals(logpath, suffix, mu, itrs, timespan)
    addpath(genpath('../../line'));
    
    disp("Simulating")
    
    % Create the model
    D = {}; 
    D{1} = Coxian.fitMeanAndSCV(mu, 1.0);
    D{2} = Coxian.fitMeanAndSCV(1.0, 1.0);
    
    model = Network('model');

    node{1} = Source(model, ['Source']);
    node{2} = Queue(model, ['Queue'], SchedStrategy.PS);
    node{3} = Sink(model, ['Sink']);
   
    jobclass{1} = OpenClass(model, 'Class1');
    
    node{1}.setArrival(jobclass{1}, D{1});
    
    node{2}.setService(jobclass{1}, D{2});
    node{2}.setNumberOfServers(1);
    
    P = model.initRoutingMatrix;

    P{1,1} = [  
        0.0, 1.0, 0.0;
        0.0, 0.0, 1.0;
        0.0, 0.0, 0.0;
    ];

    model.link(P);
    
    % Draw transient values from the model using JMT
    t = {};
    q = {};
    
    for k = 1:itrs
        seed = 1 + 123*k;
        solver = SolverJMT(model, 'seed', seed, 'timespan', timespan);
        res = solver.sampleSysAggr();
        t{k} = res.t;
        q{k} = res.state{2};
        
        if any(q{k} < 0)
            t{k} = [0];
            q{k} = [0];
        end
        
    end

    % Export data
    save([logpath, '/transientVals_', suffix, '.mat'], 't', 'q');
end