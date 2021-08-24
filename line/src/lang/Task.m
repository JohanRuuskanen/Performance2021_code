classdef Task < LayeredNetworkElement
    % A software server in a LayeredNetwork.
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    properties
        multiplicity;       %int
        replication;       %int
        scheduling;         %string
        thinkTime;
        thinkTimeMean;      %double
        thinkTimeSCV;       %double
        entries = [];
        activities = [];
        precedences = [];
        replyEntry;
    end
    
    
    methods
        %public methods, including constructor
        
        %constructor
        function self = Task(model, name, multiplicity, scheduling, thinkTime)
            % self = TASK(MODEL, NAME, MULTIPLICITY, SCHEDULING, THINKTIME)
            
            if ~exist('name','var')
                line_error(mfilename,'Constructor requires to specify at least a name.');
            end
            self@LayeredNetworkElement(name);
            
            if ~exist('multiplicity','var')
                multiplicity = 1;
            end
            if ~exist('scheduling','var')
                scheduling = SchedStrategy.INF;
            end
            if ~exist('thinkTime','var')
                thinkTime = Distrib.Zero;
            end
            self.replication = 1;
            self.multiplicity = multiplicity;
            self.scheduling = scheduling;
            self.setThinkTime(thinkTime);
            model.tasks{end+1} = self;
            switch scheduling
                case 'ref'
                    model.reftasks{end+1} = self;
            end            
        end
        
        function self = setReplication(self, replication)
            self.replication = replication;
        end
                
        function self = on(self, parent)
            % self = ON(self, PARENT)
            
            parent.addTask(self);
        end
        
        function self = setAsReferenceTask(self)
            % self = SETASREFERENCETASK(self)
            
            self.scheduling = SchedStrategy.REF;
        end
        
        function self = setThinkTime(self, thinkTime)
            % self = SETTHINKTIME(self, THINKTIME)
            
            if isnumeric(thinkTime)
                if thinkTime <= Distrib.Zero
                    self.thinkTime = Immediate();
                    self.thinkTimeMean = Distrib.Zero;
                    self.thinkTimeSCV = Distrib.Zero;
                else
                    self.thinkTime = Exp(1/thinkTime);
                    self.thinkTimeMean = thinkTime;
                    self.thinkTimeSCV = 1.0;
                end
            elseif isa(thinkTime,'Distrib')
                self.thinkTime = thinkTime;
                self.thinkTimeMean = thinkTime.getMean();
                self.thinkTimeSCV = thinkTime.getSCV();
            end
        end
        
        %addEntry
        function self = addEntry(self, newEntry)
            % self = ADDENTRY(self, NEWENTRY)
            
            self.entries = [self.entries; newEntry];
        end
        
        %addActivity
        function self = addActivity(self, newAct)
            % self = ADDACTIVITY(self, NEWACT)
            
            newAct.setParent(self.name);
            self.activities = [self.activities; newAct];
        end
        
        %setActivity
        function self = setActivity(self, newAct, index)
            % self = SETACTIVITY(self, NEWACT, INDEX)
            
            self.activities(index,1) = newAct;
        end
        
        %removeActivity
        function self = removeActivity(self, index)
            % self = REMOVEACTIVITY(self, INDEX)
            
            idxToKeep = [1:index-1,index+1:length(self.activities)];
            self.activities = self.activities(idxToKeep);
            self.actNames = self.actNames(idxToKeep);
        end
        
        %addPrecedence
        function self = addPrecedence(self, newPrec)
            % self = ADDPRECEDENCE(self, NEWPREC)
            
            if iscell(newPrec)
                for m=1:length(newPrec)
                    self.precedences = [self.precedences; newPrec{m}];
                end
            else
                self.precedences = [self.precedences; newPrec];
            end
        end
        
        %setReplyEntry
        function self = setReplyEntry(self, newReplyEntry)
            % self = SETREPLYENTRY(self, NEWREPLYENTRY)
            
            self.replyEntry = newReplyEntry;
        end
        
        function meanHostDemand = getMeanHostDemand(self, entryName)
            % MEANHOSTDEMAND = GETMEANHOSTDEMAND(self, ENTRYNAME)
            
            % determines the demand posed by the entry entryName
            % the demand is located in the activity of the corresponding entry
            
            meanHostDemand = -1;
            for j = 1:length(self.entries)
                if strcmp(self.entries(j).name, entryName)
                    meanHostDemand = self.entries(j).activities(1).hostDemandMean;
                    break;
                end
            end
        end
        
    end
    
end