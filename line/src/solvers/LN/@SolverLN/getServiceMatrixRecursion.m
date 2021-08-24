function U = getServiceMatrixRecursion(self, lqn, aidx, U)
% auxiliary function to getServiceMatrix
nextaidxs = find(lqn.graph(aidx,:));
for nextaidx = nextaidxs
    if ~isempty(nextaidx)
        if ~(lqn.parent(aidx) == lqn.parent(nextaidx))
            %% if the successor activity is a call
            for cidx = lqn.callsof{aidx}
                switch lqn.calltype(cidx)
                    case CallType.ID_SYNC
                        U(aidx,lqn.nidx+cidx) = lqn.callproc{cidx}.getMean;
                    case CallType.ID_ASYNC
                        % nop - doesn't contribute to respt
                end
            end
        end
        % here we have processed all calls, let us do the activities now
        %% if the successor activity is not a call
        if (lqn.parent(aidx) == lqn.parent(nextaidx))
            U(aidx,nextaidx) = full(lqn.graph(aidx,nextaidx));
            %% now recursively build the rest of the routing matrix graph
            U = self.getServiceMatrixRecursion(lqn, nextaidx, U);
        end
    end
end % nextaidx
end
