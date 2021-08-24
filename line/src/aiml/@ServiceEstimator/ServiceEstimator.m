classdef ServiceEstimator < handle
    % Abstract class for distribution
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    properties (Access = public)
        options; % Data structure with engine options
        model; % Model to be analyzed
        samples; % Input dataset
        samplesAggr; % Input dataset
    end
    
    methods (Hidden)
        %Constructor
        function self = ServiceEstimator(model, options)
            % SELF = SOLVER(MODEL, NAME, OPTIONS)
            if ~exist('options','var')
                options = self.defaultOptions();
            end
            self.model = model;
            self.options = options;
            self.samples = cell(model.getNumberOfNodes, model.getNumberOfClasses);
            self.samplesAggr = cell(model.getNumberOfNodes, 1);
        end
    end
    
    methods
        function self = addSamples(self, sampleData)
            i = self.model.getNodeIndex(sampleData.node);
            r = 0;
            if sampleData.isAggregate
                %jobclass = NaN;
                self.samplesAggr{i}{end+1} = sampleData;
            else
                r = self.model.getClassIndex(sampleData.class);
                % store
                self.samples{i,r}{end+1} = sampleData;
            end
        end
        
        function data = getData(self)
            data = self.samples;
        end
        
        function data = getDataAggr(self)
            data = self.samplesAggr;
        end
        
        % look into the data available for that node and class for the
        % first ArvR dataset
        function data = getArvR(self, node, jobclass)
            data = [];
            i = self.model.getNodeIndex(node);
            r = self.model.getClassIndex(jobclass);
            nodeData = self.samples{i,r};
            for d=1:length(nodeData)
                if strcmp(nodeData{d}.type, Metric.ArvR)
                    data = nodeData{d};
                    return
                end
            end
        end
        
        % look into the data available for that node and class for the
        % first ArvR dataset
        function data = getUtil(self, node, jobclass)
            data = [];
            i = self.model.getNodeIndex(node);
            r = self.model.getClassIndex(jobclass);
            nodeData = self.samples{i,r};
            for d=1:length(nodeData)
                if strcmp(nodeData{d}.type, Metric.Util)
                    data = nodeData{d};
                    return
                end
            end
        end
        
        % look into the data available for that node and class for the
        % first ArvR dataset
        function data = getRespT(self, node, jobclass)
            data = [];
            i = self.model.getNodeIndex(node);
            r = self.model.getClassIndex(jobclass);
            nodeData = self.samples{i,r};
            for d=1:length(nodeData)
                if strcmp(nodeData{d}.type, Metric.RespT)
                    data = nodeData{d};
                    return
                end
            end
        end
        
        % look into the data available for that node and class for the
        % first ArvR dataset
        function data = getAggrUtil(self, node)
            data = [];
            i = self.model.getNodeIndex(node);
            nodeAggrData = self.samplesAggr{i};
            for d=1:length(nodeAggrData)
                if strcmp(nodeAggrData{d}.type, Metric.Util)
                    data = nodeAggrData{d};
                    return
                end
            end
        end
        
        % look into the data available for that node and class for the
        % first QLen dataset
        function data = getQLen(self, node, jobclass, ev)
            data = [];
            i = self.model.getNodeIndex(node);
            r = self.model.getClassIndex(jobclass);
            nodeData = self.samples{i,r};
            if isempty(ev)
                for d=1:length(nodeData)
                    if strcmp(nodeData{d}.type, Metric.QLen)
                        data = nodeData{d};
                        return
                    end
                end
            else
                for d=1:length(nodeData)
                    % e.g., arrival queue-length
                    if strcmp(nodeData{d}.type, Metric.QLen) && (nodeData{d}.event.node == ev.node) && (nodeData{d}.event.class == ev.class) && (nodeData{d}.event.node == ev.event)
                        data = nodeData{d};
                        return
                    end
                end
            end
        end
        
        % look into the data available for that node and class for the
        % first QLen dataset
        function data = getAggrQLen(self, node, ev)
            data = [];
            i = self.model.getNodeIndex(node);
            nodeData = self.samplesAggr{i};
            if isempty(ev)
                for d=1:length(nodeData)
                    if strcmp(nodeData{d}.type, Metric.QLen)
                        data = nodeData{d};
                        return
                    end
                end
            else
                for d=1:length(nodeData)
                    % e.g., arrival queue-length
                    if strcmp(nodeData{d}.type, Metric.QLen) && (nodeData{d}.cond.node == ev.node) && (nodeData{d}.cond.class == ev.class) && (nodeData{d}.cond.event == ev.event)
                        data = nodeData{d};
                        return
                    end
                end
            end
        end
        
        interpolate(self);
        estVal = estimateAt(self, node);
        estVal = estimatorUBO(self, nodes);
        estVal = estimatorUBR(self, node);
    end
    
    methods (Static)
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            % Return default options
            options = struct();
            options.verbose = 1;
            options.method = 'ubr';
            options.iter_max = 1000;
        end
    end
end