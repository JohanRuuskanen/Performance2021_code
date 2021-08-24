function ret = tget(AvgTable,station,class)
if ~isstr(station) % inputs are objects
    if nargin==2
        if isa(station,'JobClass')
            class = station;
            station=[];
        else
            class=[];
        end
    end
    if isempty(station)
        ret = AvgTable(AvgTable.JobClass == class.name,:);
    elseif isempty(class)
        ret = AvgTable(AvgTable.Station == station.name,:);
    else
        ret = AvgTable(AvgTable.Station == station.name & AvgTable.JobClass == class.name,:);
    end
else % inputs are strings
    inputstring = station;
    if nargin==2
        ret = AvgTable(AvgTable.Station == inputstring,:);
        if isempty(ret)
            ret = AvgTable( AvgTable.JobClass == inputstring,:);
        end
    else        
        ret = AvgTable(AvgTable.Station == station & AvgTable.JobClass == class,:);
    end
end
end