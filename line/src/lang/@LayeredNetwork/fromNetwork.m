function myLN = fromNetwork(model)

% myLN = LayeredNetwork('LQN1');
% qn = model.getStruct;
% 
% P={}; T={}; E={};
% for i=1:qn.nstations
%     P{i} = Processor(myLN,['P',num2str(i)], 1, SchedStrategy.INF);
%     T{i} = Task(myLN,['T',num2str(i)], 1, qn.sched{i}).on(P{i});
%     for r=1:qn.nclasses
%         E{i,r} = Entry(myLN,['E',num2str(i),num2str()]).on(T{i});
%     end
% end
% 
% for i=1:M
%     for r=1:R
%         % definition of processors, tasks and entries
%         T0 = Task(myLN, 'T0', N1, SchedStrategy.REF).on(P0);
%         E0 = Entry(myLN, 'E0').on(T0);
%         T1 = Task(myLN, 'T1', N2, SchedStrategy.REF).on(P0);
%         E1 = Entry(myLN, 'E1').on(T1);
%         
%         P2 = Processor(myLN, 'P2', 1, );
%         T2 = Task(myLN, 'T2', scale, SchedStrategy.FCFS).on(P2);
%         E2 = Entry(myLN, 'E2').on(T2);
%         E3 = Entry(myLN, 'E3').on(T2);
%         
%         % definition of activities
%         T0.setThinkTime(Exp.fitMean(Z1));
%         T1.setThinkTime(Exp.fitMean(Z2));
%         
%         A0 = Activity(myLN, 'A0', Immediate()).on(T0).boundTo(E0).synchCall(E2,1);
%         A1 = Activity(myLN, 'A1', Immediate()).on(T1).boundTo(E1).synchCall(E3,1);
%         A2 = Activity(myLN, 'A2', APH.fitMeanAndSCV(c/50,SCV)).on(T2).boundTo(E2).repliesTo(E2);
%         A3 = Activity(myLN, 'A3', APH.fitMeanAndSCV(1,SCV)).on(T2).boundTo(E3).repliesTo(E3);
%         
%         % instantiate solvers
%         options = SolverLQNS.defaultOptions;
%         options.keep = true;
%         options.verbose = 1;
%         %options.method = 'lqsim';
%         %options.samples = 1e4;
%         lqnssolver = SolverLQNS(myLN, options);
%         AvgTableLQNS = lqnssolver.getAvgTable;
%         lqn(c)= (N1-AvgTableLQNS.Tput(3)*Z1)  + (N2-AvgTableLQNS.Tput(4)*Z2)
%     end
% end

end