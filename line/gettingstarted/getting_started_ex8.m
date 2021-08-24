function p_opt = getting_started_ex8()
model = Network('LoadBalCQN');
% Block 1: nodes
delay = Delay(model,'Think');
queue1 = Queue(model, 'Queue1', SchedStrategy.PS);
queue2 = Queue(model, 'Queue2', SchedStrategy.PS);
% Block 2: classes
cclass = ClosedClass(model, 'Job1', 16, delay);
delay.setService(cclass, Exp(1));
queue1.setService(cclass, Exp(0.75));
queue2.setService(cclass, Exp(0.50));
% Block 3: topology
P = model.initRoutingMatrix;
P{cclass}(queue1, delay) = 1.0;
P{cclass}(queue2, delay) = 1.0;

% Block 4: solution
    function R = objFun(p)
        P{cclass}(delay, queue1) = p;
        P{cclass}(delay, queue2) = 1-p;
        model.link(P);
        R = SolverMVA(model,'exact','verbose',false).getAvgSysRespT;
    end
p_opt = fminbnd(@(p) objFun(p), 0,1)
end