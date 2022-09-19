function simulateExample2_transVals(logpath, suffix, mu, itrs, timespan)
    addpath(genpath('../../line'));
    
    disp("Simulating")

    D = {}; 
    D{1} = Coxian.fitMeanAndSCV(mu, 1.0);
    D{2} = Coxian.fitMeanAndSCV(0.5, 0.5);
    D{3} = Coxian.fitMeanAndSCV(1.0, 10.0);
    
    model = Network('model');

    node{1} = Delay(model, 'Delay1'); 
    node{2} = Queue(model, 'Queue2', SchedStrategy.PS);
    node{3} = Queue(model, 'Queue3', SchedStrategy.PS);
   
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

    model.link(P);
    
    t = {};
    q = {};
    
    for k = 1:itrs
        seed = 1 + 123*k;
        solver = SolverJMT(model, 'seed', seed, 'timespan', timespan);
        res = solver.sampleSysAggr();
        t{k} = res.t;
        q{k} = cell2mat(res.state);
        
        if any(q{k} < 0)
            t{k} = [];
            q{k} = [];
        end
        
        qnan = isnan(q{k});
        if sum(sum(qnan) > 0)
            disp("Found NaN values")
            if any(sum(qnan, 2) > 1)
                error("Too many nans on rows, cannot interpolate!")
            end
            [cols, rows] = find(qnan);
            disp("Interpolating...")
            for i = 1:length(cols)
                q{k}(cols(i), rows(i)) = 50 - nansum(q{k}(cols(i), :));
            end
        end
        
    end

    % Export data
    save([logpath, '/transientVals_', suffix, '.mat'], 't', 'q');

    
    
end