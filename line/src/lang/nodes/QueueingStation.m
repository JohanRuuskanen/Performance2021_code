classdef QueueingStation < Queue
    % Alias for the Queue class
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    methods
        %Constructor
        function self = QueueingStation(model, name, schedStrategy)
            % SELF = QUEUEINGSTATION(MODEL, NAME, SCHEDSTRATEGY)
            
            self@Queue(model, name, schedStrategy);
        end
    end
end
