classdef StatelessClassSwitcher < ClassSwitcher
    % A class switch node basic on class switch probabilities
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    methods
        %Constructor
        function self = StatelessClassSwitcher(classes, csMatrix)
            % SELF = STATELESSCLASSSWITCHER(CLASSES, CSMATRIX)
            
            self@ClassSwitcher(classes, 'StatelessClassSwitcher');
            % this is slower than indexing the matrix, but it is a small
            % matrix anyway
            self.csFun = @(r,s,state,statep) csMatrix(r,s); % state parameter if present is ignored
        end
        
        function self = updateClassSwitch(self, csMatrix)
            self.csFun = @(r,s,state,statep) csMatrix(r,s); % state parameter if present is ignored
        end
    end
    
end
