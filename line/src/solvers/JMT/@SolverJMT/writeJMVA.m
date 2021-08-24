function [outputFileName] = writeJMVA(self, outputFileName)
% [OUTPUTFILENAME] = WRITEJMVA(OUTPUTFILENAME)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

if ~self.model.hasProductFormSolution
    line_error(mfilename,'JMVA requires the model to have a product-form solution.');
end

if self.model.hasClassSwitch
    %    line_error(mfilename,'JMVA does not support class switching.');
end

if isoctave
    mvaDoc = javaObject('org.apache.xerces.dom.DocumentImpl');
    mvaElem = mvaDoc.createElement('sim');
    mvaDoc.appendChild(mvaElem);
else
    mvaDoc = com.mathworks.xml.XMLUtils.createDocument('model');
    mvaElem = mvaDoc.getDocumentElement;
end

mvaElem.setAttribute('xmlns:xsi', self.xmlnsXsi);
mvaElem.setAttribute('xsi:noNamespaceSchemaLocation', 'JMTmodel.xsd');
if ~exist('outFileName','var')
    outputFileName = self.getJMVATempPath();
end

qn = self.model.getStruct();

algTypeElement = mvaDoc.createElement('algType');
switch self.options.method
    case {'jmva.recal'}
        if max(qn.nservers(isfinite(qn.nservers))) > 1
            line_error(sprintf('%s does not support multi-server stations.',self.options.method));
        end
        algTypeElement.setAttribute('name','RECAL');
    case {'jmva.comom'}
        if max(qn.nservers(isfinite(qn.nservers))) > 1
            line_error(sprintf('%s does not support multi-server stations.',self.options.method));
        end
        algTypeElement.setAttribute('name','CoMoM');
    case {'jmva.chow'}
        if max(qn.nservers(isfinite(qn.nservers))) > 1
            line_error(sprintf('%s does not support multi-server stations.',self.options.method));
        end
        algTypeElement.setAttribute('name','Chow');
    case {'jmva.bs','jmva.amva'}
        if max(qn.nservers(isfinite(qn.nservers))) > 1
            line_error(sprintf('%s does not support multi-server stations.',self.options.method));
        end
        algTypeElement.setAttribute('name','Bard-Schweitzer');
    case {'jmva.aql'}
        if max(qn.nservers(isfinite(qn.nservers))) > 1
            line_error(sprintf('%s does not support multi-server stations.',self.options.method));
        end
        algTypeElement.setAttribute('name','AQL');
    case {'jmva.lin'}
        if max(qn.nservers(isfinite(qn.nservers))) > 1
            line_error(sprintf('%s does not support multi-server stations.',self.options.method));
        end
        algTypeElement.setAttribute('name','Linearizer');
    case {'jmva.dmlin'}
        if max(qn.nservers(isfinite(qn.nservers))) > 1
            line_error(sprintf('%s does not support multi-server stations.',self.options.method));
        end
        algTypeElement.setAttribute('name','De Souza-Muntz Linearizer');
    case {'jmva.ls'}
        algTypeElement.setAttribute('name','Logistic Sampling');
    otherwise
        algTypeElement.setAttribute('name','MVA');
end
algTypeElement.setAttribute('tolerance','1.0E-7');
algTypeElement.setAttribute('maxSamples',num2str(self.options.samples));

%%%%%%%%%%
M = qn.nstations;    %number of stations
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

Lchain = zeros(M,C);
STchain = zeros(M,C);

SCVchain = zeros(M,C);
Nchain = zeros(1,C);
refstatchain = zeros(C,1);
for c=1:qn.nchains
    inchain = find(qn.chains(c,:));
    isOpenChain = any(isinf(qn.njobs(inchain)));
    for i=1:qn.nstations
        % we assume that the visits in L(i,inchain) are equal to 1
        Lchain(i,c) = Vchain(i,c) * ST(i,inchain) * alpha(i,inchain)';
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
Lchain(~isfinite(Lchain))=0;
%%%%%%%%%%
parametersElem = mvaDoc.createElement('parameters');
classesElem = mvaDoc.createElement('classes');
classesElem.setAttribute('number',num2str(qn.nchains));
stationsElem = mvaDoc.createElement('stations');
stationsElem.setAttribute('number',num2str(qn.nstations - sum(self.getStruct.nodetype == NodeType.Source)));
refStationsElem = mvaDoc.createElement('ReferenceStation');
refStationsElem.setAttribute('number',num2str(qn.nchains));
algParamsElem = mvaDoc.createElement('algParams');

sourceid = self.getStruct.nodetype == NodeType.Source;
for c=1:qn.nchains
    if isfinite(sum(qn.njobs(qn.chains(c,:))))
        classElem = mvaDoc.createElement('closedclass');
        classElem.setAttribute('population',num2str(Nchain(c)));
        classElem.setAttribute('name',sprintf('Chain%02d',c));
    else
        classElem = mvaDoc.createElement('openclass');
        classElem.setAttribute('rate',num2str(sum(qn.rates(sourceid,qn.chains(c,:)))));
        classElem.setAttribute('name',sprintf('Chain%02d',c));
    end
    classesElem.appendChild(classElem);
end

isLoadDep = false(1,qn.nstations);
for i=1:qn.nstations
    switch self.getStruct.nodetype(self.getStruct.stationToNode(i))
        case NodeType.Delay
            statElem = mvaDoc.createElement('delaystation');
            statElem.setAttribute('name',qn.nodenames{self.getStruct.stationToNode(i)});
        case NodeType.Queue
            if qn.nservers(i) == 1
                isLoadDep(i) = false;
                statElem = mvaDoc.createElement('listation');
                statElem.setAttribute('name',qn.nodenames{self.getStruct.stationToNode(i)});
                statElem.setAttribute('servers',num2str(1));
            else
                isLoadDep(i) = true;
                statElem = mvaDoc.createElement('ldstation');
                statElem.setAttribute('name',qn.nodenames{self.getStruct.stationToNode(i)});
                statElem.setAttribute('servers',num2str(1));
            end
        otherwise
            continue
    end
    srvTimesElem = mvaDoc.createElement('servicetimes');
    for c=1:qn.nchains
        if isLoadDep(i)
            statSrvTimeElem = mvaDoc.createElement('servicetimes');
            statSrvTimeElem.setAttribute('customerclass',sprintf('Chain%02d',c));
            ldSrvString = num2str(STchain(i,c));
            if any(isinf(NK))
                line_error(mfilename,'JMVA does not support open classes in load-dependent models;');
            end
            
            for n=2:sum(NK)
                ldSrvString = sprintf('%s;%s',ldSrvString,num2str(STchain(i,c)/min( n, qn.nservers(i) )));
            end
            statSrvTimeElem.appendChild(mvaDoc.createTextNode(ldSrvString));
            srvTimesElem.appendChild(statSrvTimeElem);
        else
            statSrvTimeElem = mvaDoc.createElement('servicetime');
            statSrvTimeElem.setAttribute('customerclass',sprintf('Chain%02d',c));
            statSrvTimeElem.appendChild(mvaDoc.createTextNode(num2str(STchain(i,c))));
            srvTimesElem.appendChild(statSrvTimeElem);
        end
    end
    statElem.appendChild(srvTimesElem);
    visitsElem = mvaDoc.createElement('visits');
    for c=1:qn.nchains
        statVisitElem = mvaDoc.createElement('visit');
        statVisitElem.setAttribute('customerclass',sprintf('Chain%02d',c));
        if STchain(i,c) > 0
            val = Lchain(i,c) ./ STchain(i,c);
        else
            val = 0;
        end
        statVisitElem.appendChild(mvaDoc.createTextNode(num2str(val)));
        visitsElem.appendChild(statVisitElem);
    end
    statElem.appendChild(visitsElem);
    
    stationsElem.appendChild(statElem);
end

for c=1:qn.nchains
    classRefElem = mvaDoc.createElement('Class');
    classRefElem.setAttribute('name',sprintf('Chain%d',c));
    classRefElem.setAttribute('refStation',qn.nodenames{qn.stationToNode(refstatchain(c))});
    refStationsElem.appendChild(classRefElem);
end

compareAlgsElem = mvaDoc.createElement('compareAlgs');
compareAlgsElem.setAttribute('value','false');
algParamsElem.appendChild(algTypeElement);
algParamsElem.appendChild(compareAlgsElem);

parametersElem.appendChild(classesElem);
parametersElem.appendChild(stationsElem);
parametersElem.appendChild(refStationsElem);
mvaElem.appendChild(parametersElem);
mvaElem.appendChild(algParamsElem);

xmlwrite(outputFileName, mvaDoc);
end
