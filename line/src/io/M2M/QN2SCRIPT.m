function QN2SCRIPT(model, modelName, fid)
% QN2SCRIPT(MODEL, MODELNAME, FID)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
if ~exist('modelName','var')
    modelName='myModel';
end
if ~exist('fid','var')
    fid=1;
end
if ischar(fid)
    fid = fopen(fid,'w');
end
qn = model.getStruct;
%% initialization
fprintf(fid,'model = Network(''%s'');\n',modelName);
fprintf(fid,'\n%%%% Block 1: nodes');
rt = qn.rt;
rtnodes = qn.rtnodes;
hasSink = 0;
sourceID = 0;
PH = qn.proc;

fprintf(fid,'\n');
%% write nodes
for i= 1:qn.nnodes
    switch qn.nodetype(i)
        case NodeType.Source
            sourceID = i;
            fprintf(fid,'node{%d} = Source(model, ''%s'');\n',i,qn.nodenames{i});
            hasSink = 1;
        case NodeType.Delay
            fprintf(fid,'node{%d} = DelayStation(model, ''%s'');\n',i,qn.nodenames{i});
        case NodeType.Queue
            fprintf(fid,'node{%d} = Queue(model, ''%s'', SchedStrategy.%s);\n', i, qn.nodenames{i}, SchedStrategy.toProperty(qn.sched(qn.nodeToStation(i))));
            if qn.nservers(qn.nodeToStation(i))>1
                fprintf(fid,'node{%d}.setNumServers(%d);\n', i, qn.nservers(qn.nodeToStation(i)));
            end
        case NodeType.Router
            fprintf(fid,'node{%d} = Router(model, ''%s'');\n',i,qn.nodenames{i});
        case NodeType.Fork
            fprintf(fid,'node{%d} = Fork(model, ''%s'');\n',i,qn.nodenames{i});
        case NodeType.Join
            fprintf(fid,'node{%d} = Join(model, ''%s'');\n',i,qn.nodenames{i});
        case NodeType.Sink
            fprintf(fid,'node{%d} = Sink(model, ''%s'');\n',i,qn.nodenames{i});
        case NodeType.ClassSwitch
%            csMatrix = eye(qn.nclasses);
%            fprintf(fid,'\ncsMatrix%d = zeros(%d);\n',i,qn.nclasses);
%             for k = 1:qn.nclasses
%                 for c = 1:qn.nclasses
%                     for m=1:qn.nnodes
%                         % routing matrix for each class
%                         csMatrix(k,c) = csMatrix(k,c) + rtnodes((i-1)*qn.nclasses+k,(m-1)*qn.nclasses+c);
%                     end
%                 end
%             end
%            for k = 1:qn.nclasses
%                for c = 1:qn.nclasses
%                    if csMatrix(k,c)>0
%                        fprintf(fid,'csMatrix%d(%d,%d) = %f; %% %s -> %s\n',i,k,c,csMatrix(k,c),qn.classnames{k},qn.classnames{c});
%                    end
%                end
%            end
            %fprintf(fid,'node{%d} = ClassSwitch(model, ''%s'', csMatrix%d);\n',i,qn.nodenames{i},i);
            fprintf(fid,'node{%d} = ClassSwitch(model, ''%s'', eye(%d)); %% Class switching is embedded in the routing matrix P \n',i,qn.nodenames{i},qn.nclasses);
    end
end
%% write classes
fprintf(fid,'\n%%%% Block 2: classes\n');
for k = 1:qn.nclasses
    if qn.njobs(k)>0
        if isinf(qn.njobs(k))
            fprintf(fid,'jobclass{%d} = OpenClass(model, ''%s'', %d);\n',k,qn.classnames{k},qn.classprio(k));
        else
            fprintf(fid,'jobclass{%d} = ClosedClass(model, ''%s'', %d, node{%d}, %d);\n',k,qn.classnames{k},qn.njobs(k),qn.stationToNode(qn.refstat(k)),qn.classprio(k));
        end
    else
        % if the reference node is unspecified, as in artificial classes,
        % set it to the first node where the rate for this class is
        % non-null
        iref = 0;
        for i=1:qn.nstations
            if sum(nnz(qn.proc{i,k}{1}))>0
                iref = i;
                break
            end
        end
        if isinf(qn.njobs(k))
            fprintf(fid,'jobclass{%d} = OpenClass(model, ''%s'', %d);\n',k,qn.classnames{k},qn.classprio(k));
        else
            fprintf(fid,'jobclass{%d} = ClosedClass(model, ''%s'', %d, node{%d}, %d);\n',k,qn.classnames{k},qn.njobs(k),iref,qn.classprio(k));
        end
    end
end
fprintf(fid,'\n');
%% arrival and service processes
for k=1:qn.nclasses
    for i=1:qn.nstations
        if qn.nodetype(qn.stationToNode(i)) ~= NodeType.Join
            if isprop(model.stations{i},'serviceProcess') && strcmp(class(model.stations{i}.serviceProcess{k}),'Replayer')
                switch qn.sched(i)
                    case SchedStrategy.EXT
                        fprintf(fid,'node{%d}.setArrival(jobclass{%d}, Replayer(''%s'')); %% (%s,%s)\n',qn.stationToNode(i),k,model.stations{i}.serviceProcess{k}.params{1}.paramValue,qn.nodenames{qn.stationToNode(i)},qn.classnames{k});
                    otherwise
                        fprintf(fid,'node{%d}.setService(jobclass{%d}, Replayer(''%s'')); %% (%s,%s)\n',qn.stationToNode(i),k,model.stations{i}.serviceProcess{k}.params{1}.paramValue,qn.nodenames{qn.stationToNode(i)},qn.classnames{k});
                end
            else
                SCVik = map_scv(PH{i,k});
                if SCVik >= 0.5
                    switch qn.sched(i)
                        case SchedStrategy.EXT
                            if SCVik == 1
                                fprintf(fid,'node{%d}.setArrival(jobclass{%d}, Exp.fitMean(%f)); %% (%s,%s)\n',qn.stationToNode(i),k,map_mean(PH{i,k}),qn.nodenames{qn.stationToNode(i)},qn.classnames{k});
                            else
                                fprintf(fid,'node{%d}.setArrival(jobclass{%d}, Cox2.fitMeanAndSCV(%f,%f)); %% (%s,%s)\n',qn.stationToNode(i),k,map_mean(PH{i,k}),SCVik,qn.nodenames{qn.stationToNode(i)},qn.classnames{k});
                            end
                        otherwise
                            if SCVik == 1
                                fprintf(fid,'node{%d}.setService(jobclass{%d}, Exp.fitMean(%f)); %% (%s,%s)\n',qn.stationToNode(i),k,map_mean(PH{i,k}),qn.nodenames{qn.stationToNode(i)},qn.classnames{k});
                            else
                                fprintf(fid,'node{%d}.setService(jobclass{%d}, Cox2.fitMeanAndSCV(%f,%f)); %% (%s,%s)\n',qn.stationToNode(i),k,map_mean(PH{i,k}),SCVik,qn.nodenames{qn.stationToNode(i)},qn.classnames{k});
                            end
                    end
                else
                    % this could be made more precised by fitting into a 2-state
                    % APH, especially if SCV in [0.5,0.1]
                    nPhases = max(1,round(1/SCVik));
                    switch qn.sched(i)
                        case SchedStrategy.EXT
                            if isnan(PH{i,k}{1})
                                fprintf(fid,'node{%d}.setArrival(jobclass{%d}, Disabled()); %% (%s,%s)\n',qn.stationToNode(i),k,qn.nodenames{qn.stationToNode(i)},qn.classnames{k});
                            else
                                fprintf(fid,'node{%d}.setArrival(jobclass{%d}, Erlang(%f,%f)); %% (%s,%s)\n',qn.stationToNode(i),k,nPhases/map_mean(PH{i,k}),nPhases,qn.nodenames{qn.stationToNode(i)},qn.classnames{k});
                            end
                        otherwise
                            if isnan(PH{i,k}{1})
                                fprintf(fid,'node{%d}.setService(jobclass{%d}, Disabled()); %% (%s,%s)\n',qn.stationToNode(i),k,qn.nodenames{qn.stationToNode(i)},qn.classnames{k});
                            else
                                fprintf(fid,'node{%d}.setService(jobclass{%d}, Erlang(%f,%f)); %% (%s,%s)\n',qn.stationToNode(i),k,nPhases/map_mean(PH{i,k}),nPhases,qn.nodenames{qn.stationToNode(i)},qn.classnames{k});
                            end
                    end
                end
            end
        end        
    end
end

fprintf(fid,'\n%%%% Block 3: topology');
if hasSink
    rt(qn.nstations*qn.nclasses+(1:qn.nclasses),qn.nstations*qn.nclasses+(1:qn.nclasses)) = zeros(qn.nclasses);
    for k=find(isinf(qn.njobs))' % for all open classes
        for i=1:qn.nstations
            % all open class transitions to ext station are re-routed to sink
            rt((i-1)*qn.nclasses+k, qn.nstations*qn.nclasses+k) = rt((i-1)*qn.nclasses+k, (sourceID-1)*qn.nclasses+k);
            rt((i-1)*qn.nclasses+k, (sourceID-1)*qn.nclasses+k) = 0;
        end
    end
end

fprintf(fid,'\n');
fprintf(fid,'P = model.initRoutingMatrix(); %% initialize routing matrix \n',qn.nclasses,qn.nclasses,qn.nnodes,qn.nnodes);
for k = 1:qn.nclasses
    for c = 1:qn.nclasses
        for i=1:qn.nnodes
            for m=1:qn.nnodes
                % routing matrix for each class
                myP{k,c}(i,m) = rtnodes((i-1)*qn.nclasses+k,(m-1)*qn.nclasses+c);
                if myP{k,c}(i,m) > 0 && qn.nodetype(i) ~= NodeType.Sink
                    % do not change %d into %f to avoid round-off errors in
                    % the total probability
                    fprintf(fid,'P{%d,%d}(%d,%d) = %d; %% (%s,%s) -> (%s,%s)\n',k,c,i,m,myP{k,c}(i,m),qn.nodenames{i},qn.classnames{k},qn.nodenames{m},qn.classnames{c});
                end
            end
        end
        %fprintf(fid,'P{%d,%d} = %s;\n',k,c,mat2str(myP{k,c}));
    end
end

fprintf(fid,'model.link(P);\n');
%if fid~=1
%    fclose(fid);
%end
end
