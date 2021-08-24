function [Pnir,runtime,fname] = solver_ctmc_jointaggr(qn, options)
% [PNIR,RUNTIME,FNAME] = SOLVER_CTMC_JOINTAGGR(QN, OPTIONS)
%
% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

M = qn.nstations;    %number of stations
K = qn.nclasses;    %number of classes
fname = '';
rt = qn.rt;
Tstart = tic;

myP = cell(K,K);
for k = 1:K
    for c = 1:K
        myP{k,c} = zeros(qn.nstations);
    end
end

for i=1:qn.nstations
    for j=1:qn.nstations
        for k = 1:K
            for c = 1:K
                % routing table for each class
                myP{k,c}(i,j) = rt((i-1)*K+k,(j-1)*K+c);
            end
        end
    end
end

[Q,~,SSq,~,~,~,qn] = solver_ctmc(qn, options);
% SSq is an aggregate state space
if options.keep
    fname = tempname;
    save([fname,'.mat'],'Q','SSq')
    line_printf('\nCTMC generator and aggregate state space saved in: ');
    line_printf([fname, '.mat'])
end
pi = ctmc_solve(Q);
pi(pi<1e-14)=0;

state = qn.state;
nvec = [];
for i=1:qn.nstations
    if qn.isstateful(i)
        isf = qn.stationToStateful(i);
        [~,nir,~,~] = State.toMarginal(qn, isf, state{isf}, options);
        nvec = [nvec, nir(:)'];
    end
end
Pnir = sum(pi(findrows(SSq,nvec)));

runtime = toc(Tstart);

if options.verbose > 0
    line_printf('\nCTMC analysis completed. Runtime: %f seconds.\n',runtime);
end
end
