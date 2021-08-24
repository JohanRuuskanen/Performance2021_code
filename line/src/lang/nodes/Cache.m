classdef Cache < StatefulNode
    % A class switch node based on cache hits or misses
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    properties
        cap;
        schedPolicy;
        schedStrategy;
        replacementPolicy;
        popularity;
        nLevels;
        itemLevelCap;
        items;
        accessProb;
    end
    
    methods
        %Constructor
        function self = Cache(model, name, nitems, itemLevelCap, replPolicy)
            % SELF = CACHE(MODEL, NAME, NITEMS, ITEMLEVELCAP, REPLPOLICY)
            
            self@StatefulNode(name);
            if ~exist('levels','var')
                levels = 1;
            end
            classes = model.classes;
            self.input = Buffer(classes);
            self.output = Dispatcher(classes);
            self.schedPolicy = SchedStrategyType.NP;
            self.schedStrategy = SchedStrategy.FCFS;
            self.items = ItemSet(model, [name,'_','Items'], nitems, self);
            self.nLevels = nnz(itemLevelCap);
            self.cap = Inf; % job capacity
            self.accessProb = {};
            self.itemLevelCap = itemLevelCap; % item capacity
            if sum(itemLevelCap) > nitems
                line_error(mfilename,'The number of items is smaller than the capacity of %s.',name);
            end
            self.replacementPolicy = ReplacementStrategy.toId(replPolicy);
            self.server =  CacheClassSwitcher(classes, levels, itemLevelCap); % replace Server created by Queue
            self.popularity = {};
            self.setModel(model);
            self.model.addNode(self);
        end
        
        %         function setMissTime(self, distribution)
        % SETMISSTIME(DISTRIBUTION)
        
        %             itemclass = self.items;
        %             self.server.serviceProcess{1, itemclass.index} = {[], ServiceStrategy.ID_SD, distribution};
        %         end
        %
        %         function setHitTime(self, distribution, level)
        % SETHITTIME(DISTRIBUTION, LEVEL)
        
        %             itemclass = self.items;
        %             if ~exist('level','var')
        %                 levels = 2:self.nLevels;
        %             else
        %                 levels = level;
        %             end
        %             for level = levels
        %                 self.server.serviceProcess{1+level, itemclass.index} = {[], ServiceStrategy.ID_SD, distribution};
        %             end
        %         end

        function self = reset(self)
            % SELF = RESET()
            %
            % Reset internal data structures when the network model is
            % reset
            self.server.actualHitProb = sparse([]);
            self.server.actualMissProb = sparse([]);
        end
        
        function self = setScheduling(self, class, strategy)
            %noop
        end
        
        
        function setHitClass(self, jobinclass, joboutclass)
            % SETHITCLASS(JOBINCLASS, JOBOUTCLASS)
            
            self.server.hitClass(jobinclass.index) = joboutclass.index;
        end
        
        function setMissClass(self, jobinclass, joboutclass)
            % SETMISSCLASS(JOBINCLASS, JOBOUTCLASS)
            
            self.server.missClass(jobinclass.index) = joboutclass.index;
        end
        
        function setRead(self, jobclass, distribution)
            % SETREAD(JOBCLASS, DISTRIBUTION)
            
            itemclass = self.items;
            if distribution.isDiscrete
                self.popularity{itemclass.index, jobclass.index} = distribution.copy;
                if self.popularity{itemclass.index, jobclass.index}.support(2) ~= itemclass.nitems
                    line_error(mfilename,'The reference model is defined on a number of items different from the ones used to instantiate %s.',self.name);
                end
                switch class(distribution)
                    case 'Zipf'
                        self.popularity{itemclass.index, jobclass.index}.setParam(2, 'n', itemclass.nitems, 'java.lang.Integer');
                end
                %                self.probselect(itemclass.index, jobclass.index) = probselect;
            else
                line_error(mfilename,'A discrete popularity distribution is required.');
            end
        end
        
        function setAccessProb(self, R)
            % SETACCESSCOSTS(R)
            
            self.accessProb = R;
        end
        
        
        function setProbRouting(self, class, destination, probability)
            % SETPROBROUTING(CLASS, DESTINATION, PROBABILITY)
            
            setRouting(self, class, RoutingStrategy.PROB, destination, probability);
        end
        
        function sections = getSections(self)
            % SECTIONS = GETSECTIONS()
            
            sections = {self.input, self.server, self.output};
        end
    end
end
