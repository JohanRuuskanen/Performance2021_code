function [RD] = solver_mam_passage_time(qn, PH, options)
% [RD] = SOLVER_MAM_PASSAGE_TIME(QN, PH, OPTIONS)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

global BuToolsVerbose;
global BuToolsCheckInput;
global BuToolsCheckPrecision;
%% generate local state spaces
M = qn.nstations;
K = qn.nclasses;
N = qn.njobs';
rt = qn.rt;
V = qn.visits;

Q = zeros(M,K);
U = zeros(M,K);
R = zeros(M,K);
T = zeros(M,K);
C = zeros(1,K);
X = zeros(1,K);

if M==2 && all(isinf(N))
    % open queueing system (one node is the external world)
    BuToolsVerbose = false;
    BuToolsCheckInput = true;
    BuToolsCheckPrecision = 1e-12;
    pie = {};
    S = {};
    for i=1:M
        switch qn.sched(i)
            case SchedStrategy.EXT
                na = cellfun(@(x) length(x{1}),{PH{i,:}});
                A = {PH{i,1}{1},PH{i,1}{2},PH{i,1}{2}};
                for k=2:K
                    A = mmap_super(A,{PH{i,k}{1},PH{i,k}{2},PH{i,k}{2}});
                end
                idx_arv = i;
            case {SchedStrategy.FCFS, SchedStrategy.HOL}
                row = size(S,1) + 1;
                for k=1:K
                    PH{i,k} = map_scale(PH{i,k}, map_mean(PH{i,k})/qn.nservers(i));
                    pie{k} = map_pie(PH{i,k});
                    S{k} = PH{i,k}{1};
                end
                idx_q = i;
            otherwise
                line_error(mfilename,'Unsupported scheduling strategy');
        end
    end
    
    if any(qn.classprio ~= qn.classprio(1)) % if priorities are not identical
        [uK,iK] = unique(qn.classprio);
        if length(uK) == length(qn.classprio) % if all priorities are different
            %            [Ret{1:2*K}] = MMAPPH1NPPR({A{[1,3:end]}}, {pie{:}}, {S{:}}, 'stDistrPH');
            line_error(mfilename,'Response time distribution in priority models not yet supported.');
        else
            line_error(mfilename,'SolverMAM requires either identical priorities or all distinct priorities');
        end
    else
        RD = cell(K,1);
        [Ret{1:2*K}] = MMAPPH1FCFS({A{[1,3:end]}}, {pie{:}}, {S{:}}, 'stDistrPH');
    end
    for k=1:K
        alpha{k} = Ret{(k-1)*2+1};
        D0{k} = Ret{(k-1)*2+2};
        RDph{k}={D0{k},(-D0{k})*ones(length(alpha{k}),1)*alpha{k}(:)'};
        % now estimate the range of the CDF
        sigma(k) = sqrt(map_var(RDph{k}));
        mean(k) = map_mean(RDph{k});
        n = 5;
        while map_cdf(RDph{k},mean(k)+n*sigma(k)) < 1-Distrib.Zero
            n = n+1;
        end
        % generate 10000 CDF points
        X = linspace(0,(mean(k)+n*sigma(k)),10000);
        F = map_cdf(RDph{k},X);
        RD{idx_arv,k} = [];
        RD{idx_q,k} = [F(:),X(:)];
    end
else
    line_warning(mfilename,'This model is not supported by SolverMAM yet. Returning with no result.');
end

end
