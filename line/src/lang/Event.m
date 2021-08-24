classdef Event 
    % A generic event occurring in a Network.
    %
    % Object of the Event class are not passed by handle.
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
        
    properties
        node;
        event;
        class;
        prob;
        state; % state information when the event occurs (optional)
        t; % timestamp when the event occurs (optional)
    end
    
    methods
        function self = Event(event, node, class, prob, state, t)
            % SELF = EVENT(EVENT, NODE, CLASS, PROB, STATE, TIMESTAMP)
            
            self.node = node;
            self.event = event;
            self.class = class;
            if ~exist('prob','var')
                prob = NaN;
            end
            self.prob = prob;
            if ~exist('state','var')
                state = []; % local state of the node or environment transition
            end
            self.state = state;
            
            if ~exist('ts','var')
                t = NaN; % timestamp
            end
            self.t = t;
        end
        
        function print(self)
            % PRINT()
            
            if isnan(self.t)
                line_printf('\n(%s: node: %d, class: %d)',EventType.toText(self.event),self.node,self.class);
            else
                line_printf('\n(%s: node: %d, class: %d, time: %d)',EventType.toText(self.event),self.node,self.class,self.t);
            end
        end
    end
    
    
end
