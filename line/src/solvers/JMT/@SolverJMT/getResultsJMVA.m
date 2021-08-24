function [result, parsed] = getResultsJMVA(self)
% [RESULT, PARSED] = GETRESULTSJMVA()

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

try
    fileName = strcat(self.getFilePath(),'jmva',filesep,self.getFileName(),'.jmva-result.jmva');
    if exist(fileName,'file')
        Pref.Str2Num = 'always';
        parsed = xml_read(fileName,Pref);
    else
        line_error(mfilename,'JMT did not output a result file, the analysis has likely failed.');
    end
catch me
    line_error(mfilename,'Unknown error upon parsing JMT result file. ');
end
self.result.('solver') = self.getName();
self.result.('model') = parsed.ATTRIBUTE;
self.result.('metric') = {};
self.result.('Prob') = struct();
try
    % older JMVA versions do not have the logValue field and will throw an
    % exception
    switch class(parsed.solutions.algorithm.normconst.ATTRIBUTE.logValue)
        case 'double'
            self.result.Prob.logNormConstAggr = parsed.solutions.algorithm.normconst.ATTRIBUTE.logValue;
        otherwise
            self.result.Prob.logNormConstAggr = NaN;
    end
catch
    self.result.Prob.logNormConstAggr = NaN;
end

qn = self.model.getStruct;
%%%
M = qn.nstations;    %number of stations
S = qn.nservers;
NK = qn.njobs';  % initial population per class
C = qn.nchains;
SCV = qn.scv;

% determine service times
ST = 1./qn.rates;
ST(isnan(qn.rates))=0;
SCV(isnan(SCV))=1;

alpha = zeros(qn.nstations,qn.nclasses);
Vchain = zeros(qn.nstations,qn.nchains);
for c=1:qn.nchains
    inchain = find(qn.chains(c,:));
    for i=1:qn.nstations
        Vchain(i,c) = sum(qn.visits{c}(i,inchain)) / sum(qn.visits{c}(qn.refstat(inchain(1)),inchain));
        for k=inchain
            alpha(i,k) = alpha(i,k) + qn.visits{c}(i,k) / sum(qn.visits{c}(i,inchain));
        end
    end
end
Vchain(~isfinite(Vchain))=0;
alpha(~isfinite(alpha))=0;
alpha(alpha<1e-12)=0;

STchain = zeros(M,C);

SCVchain = zeros(M,C);
Nchain = zeros(1,C);
refstatchain = zeros(C,1);
for c=1:qn.nchains
    inchain = find(qn.chains(c,:));
    isOpenChain = any(isinf(qn.njobs(inchain)));
    for i=1:qn.nstations
        % we assume that the visits in L(i,inchain) are equal to 1
        STchain(i,c) = ST(i,inchain) * alpha(i,inchain)';
        if isOpenChain && i == qn.refstat(inchain(1)) % if this is a source ST = 1 / arrival rates
            STchain(i,c) = 1 / sumfinite(qn.rates(i,inchain)); % ignore degenerate classes with zero arrival rates
        else
            STchain(i,c) = ST(i,inchain) * alpha(i,inchain)';
        end
        SCVchain(i,c) = SCV(i,inchain) * alpha(i,inchain)';
    end
    Nchain(c) = sum(NK(inchain));
    refstatchain(c) = qn.refstat(inchain(1));
    if any((qn.refstat(inchain(1))-refstatchain(c))~=0)
        line_error(sprintf('Classes in chain %d have different reference station.',c));
    end
end
STchain(~isfinite(STchain))=0;
%%%
statres = parsed.solutions.algorithm.stationresults;

for k=1:qn.nclasses
    for i=1:qn.nstations
        switch qn.nodetype(self.getStruct.stationToNode(i))
            case NodeType.Source
                s = struct();
                s.('alfa') = NaN;
                s.('analyzedSamples') = Inf;
                s.('class') = qn.classnames{k};
                s.('discardedSamples') = 0;
                s.('lowerLimit') = qn.rates(i,k);
                s.('maxSamples') = Inf;
                s.('meanValue') = qn.rates(i,k);
                s.('measureType') = Metric.Tput;
                s.('nodeType') = 'station';
                s.('precision') = Inf;
                s.('station') = qn.nodenames{self.getStruct.stationToNode(i)};
                s.('successful') = 'true';
                s.('upperLimit') = qn.rates(i,k);
                self.result.metric{end+1} = s;
                
                s = struct();
                s.('alfa') = NaN;
                s.('analyzedSamples') = Inf;
                s.('class') = qn.classnames{k};
                s.('discardedSamples') = 0;
                s.('lowerLimit') = 0;
                s.('maxSamples') = Inf;
                s.('meanValue') = 0;
                s.('measureType') = Metric.QLen;
                s.('nodeType') = 'station';
                s.('precision') = Inf;
                s.('station') = qn.nodenames{self.getStruct.stationToNode(i)};
                s.('successful') = 'true';
                s.('upperLimit') = 0;
                self.result.metric{end+1} = s;
                
                s = struct();
                s.('alfa') = NaN;
                s.('analyzedSamples') = Inf;
                s.('class') = qn.classnames{k};
                s.('discardedSamples') = 0;
                s.('lowerLimit') = 0;
                s.('maxSamples') = Inf;
                s.('meanValue') = 0;
                s.('measureType') = Metric.RespT;
                s.('nodeType') = 'station';
                s.('precision') = Inf;
                s.('station') = qn.nodenames{self.getStruct.stationToNode(i)};
                s.('successful') = 'true';
                s.('upperLimit') = 0;
                self.result.metric{end+1} = s;
                
                s = struct();
                s.('alfa') = NaN;
                s.('analyzedSamples') = Inf;
                s.('class') = qn.classnames{k};
                s.('discardedSamples') = 0;
                s.('lowerLimit') = 0;
                s.('maxSamples') = Inf;
                s.('meanValue') = 0;
                s.('measureType') = Metric.Util;
                s.('nodeType') = 'station';
                s.('precision') = Inf;
                s.('station') = qn.nodenames{self.getStruct.stationToNode(i)};
                s.('successful') = 'true';
                s.('upperLimit') = 0;
                self.result.metric{end+1} = s;
        end
    end
end

%%
%Rchain
%Xchain
%Tchain
% for c=1:qn.nchains
%     inchain = find(qn.chains(c,:));
%     for k=inchain(:)'
%         X(k) = Xchain(c) * alpha(qn.refstat(k),k);
%         for i=1:qn.nstations
%             if isinf(S(i))
%                 U(i,k) = ST(i,k) * (Xchain(c) * Vchain(i,c) / Vchain(qn.refstat(k),c)) * alpha(i,k);
%             else
%                 U(i,k) = ST(i,k) * (Xchain(c) * Vchain(i,c) / Vchain(qn.refstat(k),c)) * alpha(i,k) / S(i);
%             end
%             if Lchain(i,c) > 0
%                 Q(i,k) = Rchain(i,c) * ST(i,k) / STchain(i,c) * Xchain(c) * Vchain(i,c) / Vchain(qn.refstat(k),c) * alpha(i,k);
%                 T(i,k) = Tchain(i,c) * alpha(i,k);
%                 R(i,k) = Q(i,k) / T(i,k);
%             else
%                 T(i,k) = 0;
%                 R(i,k)=0;
%                 Q(i,k)=0;
%             end
%         end
%         C(k) = qn.njobs(k) / X(k);
%     end
% end
% Q=abs(Q); R=abs(R); X=abs(X); U=abs(U); T=abs(T); C=abs(C);
% T(~isfinite(T))=0; U(~isfinite(U))=0; Q(~isfinite(Q))=0; R(~isfinite(R))=0; X(~isfinite(X))=0; C(~isfinite(C))=0;
%%

for i=1:length(statres)
    classres = statres(i).classresults;
    for c=1:length(classres)
        inchain = find(qn.chains(c,:));
        for m=1:length(classres(c).measure)
            for k=inchain(:)'
                s = struct();
                s.('alfa') = NaN;
                s.('analyzedSamples') = Inf;
                s.('class') = qn.classnames{k};
                s.('discardedSamples') = 0;
                s.('meanValue') = classres(c).measure(m).ATTRIBUTE.meanValue;
                if strcmp(s.meanValue,'NaN')
                    s.meanValue = NaN;
                end
                s.('maxSamples') = Inf;
                s.('measureType') = classres(c).measure(m).ATTRIBUTE.measureType;
                switch classres(c).measure(m).ATTRIBUTE.measureType
                    case 'Utilization'
                        if isinf(qn.nservers(i))
                            s.meanValue = ST(i,k) * (s.meanValue / STchain(i,c)) * Vchain(i,c) / Vchain(qn.refstat(k),c) * alpha(i,k);
                        else
                            s.meanValue = ST(i,k) * (s.meanValue / STchain(i,c)) / Vchain(qn.refstat(k),c) * alpha(i,k) * min(sum(NK(isfinite(NK))), qn.nservers(i)) / qn.nservers(i);
                        end
                        s.('measureType') = classres(c).measure(m).ATTRIBUTE.measureType;
                    case 'Throughput'
                        s.meanValue = s.meanValue * alpha(i,k);
                    case 'Number of Customers'
                        % Q(i,k) = Rchain(i,c) * ST(i,k) / STchain(i,c) * Xchain(c) * Vchain(i,c) / Vchain(qn.refstat(k),c) * alpha(i,k);
                        s.meanValue = s.meanValue * ST(i,k) / STchain(i,c) / Vchain(qn.refstat(k),c) * alpha(i,k);
                    case 'Residence time'
                        s.('measureType') = 'Response Time';
                        s.meanValue = s.meanValue / qn.visits{c}(i,k); % this is to convert from JMVA's residence into LINE's response time per visit
                        if isinf(qn.nservers(i))
                            s.meanValue = s.meanValue * ST(i,k) / STchain(i,c) / Vchain(qn.refstat(k),c) * alpha(i,k);
                        else
                            s.meanValue = s.meanValue * ST(i,k) / STchain(i,c) / Vchain(qn.refstat(k),c) * alpha(i,k);
                        end
                    otherwise
                        s.('measureType') = classres(c).measure(m).ATTRIBUTE.measureType;
                end
                s.('lowerLimit') = s.meanValue;
                s.('upperLimit') = s.meanValue;
                s.('nodeType') = 'station';
                s.('precision') = Inf;
                s.('station') = statres(i).ATTRIBUTE.station;
                s.('successful') = classres(c).measure(m).ATTRIBUTE.successful;
                self.result.metric{end+1} = s;
            end
        end
    end
end
result = self.result;
end
