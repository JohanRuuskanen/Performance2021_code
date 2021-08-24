function [QNclass,UNclass,RNclass,TNclass] = getAvg(self,Q,U,R,T)
% [QNCLASS,UNCLASS,RNCLASS,TNCLASS] = GETAVG(SELF,Q,U,R,T)

% Return average station metrics at steady-state
%
% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

if nargin == 1 % no parameter
    if isempty(self.model.handles) || ~isfield(self.model.handles,'Q') || ~isfield(self.model.handles,'U') || ~isfield(self.model.handles,'R') || ~isfield(self.model.handles,'T') || ~isfield(self.model.handles,'A')
        self.resetResults(); % reset in case there are partial results saved
    end
    [Q,U,R,T] = self.model.getAvgHandles;
elseif nargin == 2
    handlers = Q;
    Q=handlers{1};
    U=handlers{2};
    R=handlers{3};
    T=handlers{4};
end

%if isfield(self.options,'timespan')
%    if isfinite(self.options.timespan(2))
%        line_error(mfilename,'The getAvg method does not support the timespan option, use the getTranAvg method instead.');
%   end
%else
%    self.options.timespan = [0,Inf];
%end

if ~self.hasAvgResults || ~self.options.cache
    try
        self.runAnalysis();
    catch ME
        switch ME.identifier
            case {'Line:FeatureNotSupportedBySolver', 'Line:ModelTooLargeToSolve', 'Line:UnspecifiedOption'}
                if self.options.verbose
                    line_printf('\n%s',ME.message);
                end
                QNclass=[];
                UNclass=[];
                RNclass=[];
                TNclass=[];
                ANclass=[];
                return
            otherwise
                rethrow(ME)
        end
    end
    if ~self.hasAvgResults
        line_error(mfilename,'Line is unable to return results for this model.');
    end
end % else return cached value

M = self.model.getNumberOfStations();
K = self.model.getNumberOfClasses();

QNclass = [];
UNclass = [];
RNclass = [];
TNclass = [];

if ~isempty(Q)
    QNclass = zeros(M,K);
    for k=1:K
        for i=1:M          
            QNclass(i,k) = Q{i,k}.get(self.result, self.model);
        end
    end
end
if ~isempty(U)
    UNclass = zeros(M,K);
    for k=1:K
        for i=1:M
            UNclass(i,k) = U{i,k}.get(self.result, self.model);
        end
    end
end
if ~isempty(R)
    RNclass = zeros(M,K);
    for k=1:K
        for i=1:M
            RNclass(i,k) = R{i,k}.get(self.result, self.model);
        end
    end
end

if ~isempty(T)
    TNclass = zeros(M,K);
    for k=1:K
        for i=1:M
            TNclass(i,k) = T{i,k}.get(self.result, self.model);
        end
    end
end


%% nan values indicate that a metric is disabled
QNclass(isnan(QNclass))=0;
UNclass(isnan(UNclass))=0;
RNclass(isnan(RNclass))=0;
TNclass(isnan(TNclass))=0;

%% set to zero entries associated to immediate transitions
QNclass(RNclass < 10/Distrib.InfRate)=0;
UNclass(RNclass < 10/Distrib.InfRate)=0;
RNclass(RNclass < 10/Distrib.InfRate)=0;

%% round to zero numerical perturbations
QNclass(QNclass < Distrib.Zero)=0;
UNclass(UNclass < Distrib.Zero)=0;
RNclass(RNclass < Distrib.Zero)=0;
TNclass(TNclass < Distrib.Zero)=0;

end
