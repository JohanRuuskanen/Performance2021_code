function [init_sol, state] = solver_fluid_initsol(qn, options) %#ok<INUSD>
% [INIT_SOL, STATE] = SOLVER_FLUID_INITSOL(QN, OPTIONS) %#OK<INUSD>

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

if ~exist('options','var')
    options = Solver.defaultOptions; %#ok<NASGU>
end

init_sol = [];
for ind=1:qn.nnodes
    if qn.isstateful(ind)
        isf = qn.nodeToStateful(ind);
        ist = qn.nodeToStation(ind);
        state_i = [];
        init_sol_i = [];  % compared to state_i, this does not track disabled classes and removes Inf entries in the Sources
        [~, nir, ~, kir_i] = State.toMarginal(qn, ind, qn.state{isf});
        switch qn.sched(ist)
            case {SchedStrategy.EXT}
                state_i(:,1) = Inf; % fluid does not model infinite buffer?
                for r=1:size(kir_i,2)
                    for k=1:length(qn.mu{ist,r})
                        state_i(:,end+1) = kir_i(:,r,k);
                        if ~isnan(qn.rates(ist,r))
                            init_sol_i(:,end+1) = kir_i(:,r,k);
                        end
                    end
                end
            case {SchedStrategy.FCFS, SchedStrategy.PS, SchedStrategy.INF, SchedStrategy.DPS, SchedStrategy.HOL}
                for r=1:size(kir_i,2)
                    for k=1:length(qn.mu{ist,r})
                        if k==1
                            state_i(:,end+1) = nir(:,r) - sum(kir_i(:,r,2:end),3); % jobs in waiting buffer are re-started phase 1
                            if ~isnan(qn.rates(ist,r))
                                init_sol_i(:,end+1) = nir(:,r) - sum(kir_i(:,r,2:end),3); % jobs in waiting buffer are re-started phase 1
                            end
                        else
                            state_i(:,end+1) = kir_i(:,r,k);
                            if ~isnan(qn.rates(ist,r))
                                init_sol_i(:,end+1) = kir_i(:,r,k);
                            end
                        end
                    end
                end
            otherwise
                line_error(mfilename,'Unsupported scheduling policy at station %d',ist);
                return
        end
        init_sol = [init_sol, init_sol_i];
        state{isf} = state_i;
    end
end
end
