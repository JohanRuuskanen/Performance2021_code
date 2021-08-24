function writeXML(self,filename)
% WRITEXML(SELF,FILENAME)
% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.DocumentBuilder;
import org.w3c.dom.Document;
import org.w3c.dom.NodeList;
import org.w3c.dom.Node;
import org.w3c.dom.Element;
import java.io.File;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import javax.xml.transform.OutputKeys;

precision = '%10.15e'; %precision for doubles
docFactory = DocumentBuilderFactory.newInstance();
docBuilder = docFactory.newDocumentBuilder();
doc = docBuilder.newDocument();

%Root Element
rootElement = doc.createElement('lqn-model');
doc.appendChild(rootElement);
rootElement.setAttribute('xmlns:xsi', 'http://www.w3.org/2001/XMLSchema-instance');
rootElement.setAttribute('xsi:noNamespaceSchemaLocation', 'lqn.xsd');
rootElement.setAttribute('name', self.getName());

for p = 1:length(self.hosts)
    %processor
    curProc = self.hosts{p};
    procElement = doc.createElement('processor');
    rootElement.appendChild(procElement);
    procElement.setAttribute('name', curProc.name);        
    procElement.setAttribute('scheduling', SchedStrategy.toText(curProc.scheduling));    
    if curProc.replication>1
        procElement.setAttribute('replication', num2str(curProc.replication));
    end
    if curProc.scheduling ~= SchedStrategy.INF
        mult = num2str(curProc.multiplicity);
        if isinf(mult), mult=1; end
        procElement.setAttribute('multiplicity', mult);
    end
    if curProc.scheduling== SchedStrategy.PS
        procElement.setAttribute('quantum', num2str(curProc.quantum));
    end
    procElement.setAttribute('speed-factor', num2str(curProc.speedFactor));
    for t=1:length(curProc.tasks)
        curTask = curProc.tasks(t);
        taskElement = doc.createElement('task');
        procElement.appendChild(taskElement);
        taskElement.setAttribute('name', curTask.name);
        taskElement.setAttribute('scheduling', SchedStrategy.toText(curTask.scheduling));
        if curTask.replication>1
        taskElement.setAttribute('replication',  num2str(curTask.replication));
        end
        if curTask.scheduling ~= SchedStrategy.INF
            taskElement.setAttribute('multiplicity', num2str(curTask.multiplicity));
        end
        if curTask.scheduling == SchedStrategy.REF
            taskElement.setAttribute('think-time', num2str(curTask.thinkTimeMean));
        end
        for e=1:length(curTask.entries)
            curEntry = curTask.entries(e);
            entryElement = doc.createElement('entry');
            taskElement.appendChild(entryElement);
            entryElement.setAttribute('name', curEntry.name);
            entryElement.setAttribute('type', 'NONE');
        end
        taskActsElement = doc.createElement('task-activities');
        taskElement.appendChild(taskActsElement);
        for a=1:length(curTask.activities)
            curAct = curTask.activities(a);
            actElement = doc.createElement('activity');
            taskActsElement.appendChild(actElement);
            actElement.setAttribute('host-demand-mean', num2str(curAct.hostDemandMean));
            actElement.setAttribute('host-demand-cvsq', num2str(curAct.hostDemandSCV));
            if ~isempty(curAct.boundToEntry)
                actElement.setAttribute('bound-to-entry', curAct.boundToEntry);
            end
            actElement.setAttribute('call-order', curAct.callOrder);
            actElement.setAttribute('name', curAct.name);
            for sd=1:length(curAct.syncCallDests)
                syncCallElement = doc.createElement('synch-call');
                actElement.appendChild(syncCallElement);
                syncCallElement.setAttribute('dest', curAct.syncCallDests(sd));
                syncCallElement.setAttribute('calls-mean', num2str(curAct.syncCallMeans(sd)));
            end
            for asd=1:length(curAct.asyncCallDests)
                asyncCallElement = doc.createElement('asynch-call');
                actElement.appendChild(asyncCallElement);
                asyncCallElement.setAttribute('dest', curAct.asyncCallDests(asd));
                asyncCallElement.setAttribute('calls-mean', num2str(curAct.asyncCallMeans(asd)));
            end
        end
        for ap=1:length(curTask.precedences)
            curActPrec = curTask.precedences(ap);
            actPrecElement = doc.createElement('precedence');
            taskActsElement.appendChild(actPrecElement);
            
            preElement = doc.createElement(curActPrec.preType);
            actPrecElement.appendChild(preElement);
            if strcmp(curActPrec.preType, ActivityPrecedence.PRE_AND) && ~isempty(curActPrec.preParams)
                preElement.setAttribute('quorum', num2str(curActPrec.preParams(1)));
            end
            for pra = 1:length(curActPrec.preActs)
                preActElement = doc.createElement('activity');
                preElement.appendChild(preActElement);
                preActElement.setAttribute('name', curActPrec.preActs{pra});
            end
            
            postElement = doc.createElement(curActPrec.postType);
            actPrecElement.appendChild(postElement);
            if strcmp(curActPrec.postType, ActivityPrecedence.POST_OR)
                for poa = 1:length(curActPrec.postActs)
                    postActElement = doc.createElement('activity');
                    postElement.appendChild(postActElement);
                    postActElement.setAttribute('name', curActPrec.postActs{poa});
                    postActElement.setAttribute('prob', num2str(curActPrec.postParams(poa)));
                end
            elseif strcmp(curActPrec.postType, ActivityPrecedence.POST_LOOP)
                for poa = 1:length(curActPrec.postActs)-1
                    postActElement = doc.createElement('activity');
                    postElement.appendChild(postActElement);
                    postActElement.setAttribute('name', curActPrec.postActs{poa});
                    postActElement.setAttribute('count', num2str(curActPrec.postParams(poa)));
                end
                postElement.setAttribute('end', curActPrec.postActs{end});
            else
                for poa = 1:length(curActPrec.postActs)
                    postActElement = doc.createElement('activity');
                    postElement.appendChild(postActElement);
                    postActElement.setAttribute('name', curActPrec.postActs{poa});
                end
            end
        end
        for e=1:length(curTask.entries)
            curEntry = curTask.entries(e);
            if ~isempty(curEntry.replyActivity)
                entryReplyElement = doc.createElement('reply-entry');
                taskActsElement.appendChild(entryReplyElement);
                entryReplyElement.setAttribute('name', curEntry.name);
                for r=1:length(curEntry.replyActivity)
                    entryReplyActElement = doc.createElement('reply-activity');
                    entryReplyElement.appendChild(entryReplyActElement);
                    entryReplyActElement.setAttribute('name', curEntry.replyActivity{r});
                end
            end
        end
    end
end

%write the content into xml file
transformerFactory = TransformerFactory.newInstance();
transformer = transformerFactory.newTransformer();
transformer.setOutputProperty(OutputKeys.INDENT, 'yes');
source = DOMSource(doc);
result = StreamResult(File(filename));
transformer.transform(source, result);
end
