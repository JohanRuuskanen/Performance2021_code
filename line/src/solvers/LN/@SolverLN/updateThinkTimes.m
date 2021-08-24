function updateThinkTimes(self, it)
if size(self.lqn.iscaller,2) > 0 % ignore task models if no callers
    %torder = randperm(self.lqn.ntasks); % randomize order in which the tasks are updated
    torder = 1:(self.lqn.ntasks); % randomize order in which the tasks are updated
    % solve all task models
    for t = torder
        tidx = self.lqn.tshift + t;
        tidx_thinkTime = self.lqn.think{tidx}.getMean;
        if ~isnan(self.idxhash(tidx)) % this skips all REF tasks
            self.tput(tidx) = self.lqn.repl(tidx)*sum(self.results{end,self.idxhash(tidx)}.TN(self.ensemble{self.idxhash(tidx)}.attribute.serverIdx,:),2); % throughput of t as a server, sum is correct since servoice only entries
            if self.lqn.sched(tidx) == SchedStrategy.INF
                self.util(tidx) = sum(self.results{end,self.idxhash(tidx)}.UN(self.ensemble{self.idxhash(tidx)}.attribute.serverIdx,:),2);
                njobs = self.lqn.mult(tidx)*self.lqn.repl(tidx);
                if isinf(njobs)
                    % this section correct lqn.mult for an infinite
                    % server by replacing it with the number of jobs
                    callers_of_tidx = find(self.lqn.taskgraph(:,tidx));
                    njobs = sum(self.lqn.mult(callers_of_tidx) .* self.lqn.repl(callers_of_tidx)); %#ok<FNDSB>
                    if isinf(njobs)
                        % if also the callers of tidx_caller are inf servers, then use
                        % an heuristic
                        njobs = sum(self.lqn.mult(isfinite(self.lqn.mult)) .* self.lqn.repl(isfinite(self.lqn.mult)));
                    end
                end
                self.idlet(tidx) = (njobs-self.util(tidx)) / self.tput(tidx) - tidx_thinkTime;
            else
                self.util(tidx) = sum(self.results{end,self.idxhash(tidx)}.UN(self.ensemble{self.idxhash(tidx)}.attribute.serverIdx,:),2); % utilization of t as a server
                njobs = self.lqn.mult(tidx)*self.lqn.repl(tidx);
                self.idlet(tidx) = njobs*abs(1-self.util(tidx)) / self.tput(tidx) - tidx_thinkTime;
            end            
            self.idletproc{tidx} = Exp.fitMean(self.idlet(tidx) + tidx_thinkTime);            
        end
    end
end
end