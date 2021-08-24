function [lGn,X,Q] = pfqn_nc(L,N,Z,varargin)
% [LGN] = PFQN_NC(L,N,Z,VARARGIN)

options = Solver.parseOptions(varargin, SolverNC.defaultOptions);
% backup initial parameters
Rin = length(N);

X=[]; Q=[];

% first remove empty classes
nnzClasses = find(N);
L = L(:,nnzClasses);
N = N(:,nnzClasses);
Z = Z(:,nnzClasses);

% then scale demands in [0,1], importat that stays before the other
% simplications in case both D and Z are all very small or very large in a
% given class, in which case the may look to filter but not if all of them
% are at the same scale
R = length(N);
scalevec = ones(1,R);
for r=1:R
    scalevec(r) = max([L(:,r);Z(:,r)]);
end
L = L ./ repmat(scalevec,size(L,1),1);
Z = Z ./ scalevec;

% remove stations with no demand
Lsum = sum(L,2);
Lmax = max(L,[],2);
demStations = find((Lmax./Lsum)>Distrib.Zero);
noDemStations = setdiff(1:size(L,1), demStations);
L = L(demStations,:);
if any(N((sum(L,1) + sum(Z,1)) == 0)>0) % if there is a class with jobs but L and Z all zero
    line_warning(mfilename,'The model has no positive demands in any class.');
    if isempty(Z) || sum(Z(:))<options.tol
        lGn = 0;
    else
        lGn = - sum(factln(N)) + sum(N.*log(sum(Z,1))) + N*log(scalevec)';
    end
    return
end

% update M and R
[M,R]=size(L);

% return immediately if degenerate case
if isempty(L) || sum(L(:))<options.tol % all demands are zero
    if isempty(Z) || sum(Z(:))<options.tol
        lGn = 0;
    else
        lGn = - sum(factln(N)) + sum(N.*log(sum(Z,1))) + N*log(scalevec)';
    end
    return
elseif M==1 && (isempty(Z) || sum(Z(:))<options.tol) % single node and no think time
    lGn = factln(sum(N)) - sum(factln(N)) + sum(N.*log(sum(L,1))) + N*log(scalevec)';
    return
end

% determine contribution from jobs that permanently loop at delay
zeroDemandClasses = find(sum(L,1)<options.tol); % all jobs in delay
nonzeroDemandClasses = setdiff(1:R, zeroDemandClasses);

if isempty(sum(Z,1)) || all(sum(Z(:,zeroDemandClasses),1)<options.tol)
    lGzdem = 0;
    Nz = 0;
else
    if isempty(zeroDemandClasses) % for old MATLAB release compatibility
        lGzdem = 0;
        Nz = 0;
    else
        Nz = N(zeroDemandClasses);
        lGzdem = - sum(factln(Nz)) + sum(Nz.*log(sum(Z(:,zeroDemandClasses),1))) + Nz*log(scalevec(zeroDemandClasses))';
    end
end
L = L(:,nonzeroDemandClasses);
N = N(nonzeroDemandClasses);
Zz = Z(:,zeroDemandClasses);
Z = Z(:,nonzeroDemandClasses);
scalevecz = scalevec(nonzeroDemandClasses);
% compute G for classes No with non-zero demand
[lGnnzdem,Xnnzdem,Qnnzdem] = sub_method(L, N, Z, options);

if isempty(Xnnzdem)
    X = [];
else
    zClasses = setdiff(1:Rin, nnzClasses);
    Xz = zeros(1,length(zClasses));
    Xnnz = zeros(1,length(nnzClasses));
    Xnnz(zeroDemandClasses) = Nz./ sum(Zz,1)./ scalevec(zeroDemandClasses);
    Xnnz(nonzeroDemandClasses) = Xnnzdem./ scalevec(nonzeroDemandClasses);
    X(1,[zClasses, nnzClasses]) = [Xz, Xnnz];
    
    Qz = zeros(size(Qnnzdem,1),length(zClasses));
    Qnnz = zeros(size(Qnnzdem,1),length(nnzClasses));
    Qnnz(:,zeroDemandClasses) = 0; % they are all in the delay
    Qnnz(:,nonzeroDemandClasses) = Qnnzdem; % Q does not require scaling
    Q(noDemStations,:) = 0;
    Q(demStations,[zClasses, nnzClasses]) = [Qz, Qnnz];
end
% scale back to original demands
lGn = lGnnzdem + lGzdem + N*log(scalevecz)';
end



function [lG,X,Q] = sub_method(L,N,Z,options)
% LG = SUB_METHOD(L,N,Z,OPTIONS)
[M,R] = size(L);
X=[];Q=[];
switch options.method
    case {'ca'}
        [~,lG] = pfqn_ca(L,N,sum(Z,1));
    case {'exact'}
        if M>=R || sum(N)>20 || sum(Z)>0
            [~,lG] = pfqn_ca(L,N,sum(Z,1));
        else
            [~,lG] = pfqn_recal(L,N,sum(Z,1));% implemented with Z=0
        end
    case {'default','adaptive'}
        if M>1
            if R <= 3 && sum(N)<50
                [~,~,~,~,lG] = pfqn_mva(L,N,sum(Z,1));
            else
                if M>R
                    [~,lG] = pfqn_kt(L,N,sum(Z,1));
                else
                    [~,lG] = pfqn_le(L,N,sum(Z,1));
                end
            end
        elseif sum(Z,1)==0 % single queue, no delay
            lG = -N*log(L)';
        else % repairman model
            if N < 10
                 [~,~,~,~,lG] = pfqn_mva(L,N,sum(Z,1));
            elseif N < 50 % otherwise numerical issues
                [~,lG] = pfqn_mmint2(L,N,sum(Z,1));
                if isnan(lG)
                    [~,lG] = pfqn_le(L,N,sum(Z,1));
                end
            else
                [~,lG] = pfqn_le(L,N,sum(Z,1));
            end
        end
    case {'sampling'}
        if M==1
            [~,lG] = pfqn_grm(L,N,sum(Z,1),options.samples);
        elseif M>R
            [~,lG] = pfqn_mci(L,N,sum(Z,1),options.samples,'imci');
        else
            [~,lG] = pfqn_ls(L,N,sum(Z,1),options.samples);
        end
    case {'mmint','pnc2'}
        if size(L,1)>1
            line_error(mfilename,'The %s method requires a model with a delay and a single queueing station.',options.method);
        end
        [~,lG] = pfqn_mmint2(L,N,sum(Z,1));
    case {'grm'}
        if size(L,1)>1
            line_error(mfilename,'The %s method requires a model with a delay and a single queueing station.',options.method);
        end
        [~,lG] = pfqn_grm(L,N,sum(Z,1),options.samples);
    case {'pana','panacea','pnc'}
        [~,lG] = pfqn_panacea(L,N,sum(Z,1));
        if isnan(lG)
            line_warning(mfilename,'Model is not in normal usage, panacea cannot continue.');
        end
    case 'gm'
        [~,lG] = pfqn_gm(L,N,sum(Z,1));        
    case 'le'
        [~,lG] = pfqn_le(L,N,sum(Z,1));
    case 'kt'
        [~,lG] = pfqn_kt(L,N,sum(Z,1));
    case 'ls'
        [~,lG] = pfqn_ls(L,N,sum(Z,1),options.samples);
    case {'mci','imci'}
        [~,lG] = pfqn_mci(L,N,sum(Z,1),options.samples,'imci');
    case {'mva'}
        [~,~,~,~,lG] = pfqn_mva(L,N,sum(Z,1));
    case 'mom'
        if length(N)>1
            try
                [~,lG,X,Q] = pfqn_mom(L,N,Z);
            catch
                % java exception, probably singular linear system
                line_warning(mfilename,'Numerical problems.');
                lG = NaN;
            end
        else
            [X,Q,~,~,lG] = pfqn_mva(L,N,Z);
        end
    case 'comom'
        if R>1
            try
                [~,lG,X,Q] = pfqn_comombtf(L,N,Z);
            catch
                % java exception, probably singular linear system
                line_warning(mfilename,'Numerical problems.');
                lG = NaN;
            end
        else
            [X,Q,~,~,lG] = pfqn_mva(L,N,Z);
        end        
    case 'propfair'
        [~,lG] = pfqn_propfair(L,N,sum(Z,1));
    case {'recal'}
        if sum(Z)>0
            line_error(mfilename,'RECAL is currently available only for models with non-zero think times.');
        end
        [~,lG] = pfqn_recal(L,N,sum(Z,1));
    case 'rgf'
        if sum(Z)>0
            line_error(mfilename,'RGF is defined only for models with non-zero think times.');
        end
        [~,lG] = pfqn_rgf(L,N);
    otherwise
        line_error(mfilename,'Unrecognized method: %s',options.method);
end
return
end