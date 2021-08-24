classdef Node < NetworkElement
    % An abstract for a node in a Network model
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    properties
        model;
        input;
        server;
        output;
    end
    
    methods(Hidden)
        %Constructor
        function self = Node(name)
            % SELF = NODE(NAME)
            
            self@NetworkElement(name);
        end
        
        function self = setModel(self, model)
            % SELF = SETMODEL(MODEL)
            
            self.model = model;
        end
        
        function self = link(self, nodeTo)
            % SELF = LINK(NODETO)
            
            self.model.addLink(self,nodeTo);
        end
        
        function self = reset(self)
            % SELF = RESET()
            %
            % Reset internal data structures when the network model is
            % reset
            
        end
    end
    
    methods
        
        function sections = getSections(self)
            % SECTIONS = GETSECTIONS()
            
            sections = {self.input, self.server, self.output};
        end
        
        function setProbRouting(self, class, destination, probability)
            % SETPROBROUTING(CLASS, DESTINATION, PROBABILITY)
            
            self.setRouting(class, RoutingStrategy.PROB, destination, probability);
        end
        
        function self = setScheduling(self, class, strategy)
            % SETSCHEDULING(CLASS, STRATEGY)
            
            self.input.inputJobClasses{class.index}{2} = strategy;
        end
        
        function setRouting(self, class, strategy, destination, probability)
            % SETROUTING(CLASS, STRATEGY, DESTINATION, PROBABILITY)
            
            switch nargin
                case 3
                    self.output.outputStrategy{1, class.index}{2} = RoutingStrategy.toType(strategy);
                case 5
                    self.output.outputStrategy{1, class.index}{2} = RoutingStrategy.toType(strategy);
                    if length(self.output.outputStrategy{1, class.index})<3
                        self.output.outputStrategy{1, class.index}{3}{1} = {destination, probability};
                    else
                        self.output.outputStrategy{1, class.index}{3}{end+1} = {destination, probability};
                    end
            end
        end
        
        function bool = hasClassSwitch(self)
            % BOOL = HASCLASSSWITCH()
            
            bool = isa(self.server,'ClassSwitcher');
        end
        
        function bool = isStateful(self)
            % BOOL = ISSTATEFUL()
            
            bool = isa(self,'StatefulNode');
        end
        
        function bool = isStation(self)
            % BOOL = ISSTATION()
            
            bool = isa(self,'Station');
        end
    end
    
    methods(Access = protected)
        % Override copyElement method:
        function clone = copyElement(self)
            % CLONE = COPYELEMENT()
            
            % Make a shallow copy of all properties
            clone = copyElement@Copyable(self);
            % Make a deep copy of each object
            clone.input = self.input.copy;
            clone.server = self.server.copy;
            clone.output = self.output.copy;
        end
    end
    
    methods (Access = public)
        function ind = subsindex(self)
            % IND = SUBSINDEX()
            
            ind = double(self.model.getNodeIndex(self.name))-1; % 0 based
        end
        
        function V = horzcat(self, varargin)
            % V = HORZCAT(VARARGIN)
            
            V = zeros(1,length(varargin));
            V(1) = 1+ self.subsindex;
            for v=1:length(varargin)
                V(1+v) = 1+varargin{v}.subsindex;
            end
        end
        
        function V = vertcat(self, varargin)
            % V = VERTCAT(VARARGIN)
            
            V = zeros(length(varargin),1);
            V(1) = 1+ self.subsindex;
            for v=1:length(varargin)
                V(1+v) = 1+varargin{v}.subsindex;
            end
        end
        
        function summary(self)
            % SUMMARY()
            
            line_printf('\nNode: <strong>%s</strong>',self.getName);
            %            self.input.summary;
            %            self.server.summary;
            %            self.output.summary;
        end
    end
end
