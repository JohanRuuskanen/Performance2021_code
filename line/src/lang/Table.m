function t = Table(varargin) 
% T = TABLE(VARARGIN)

if isoctave
    n = length(varargin{1});
    m = length(varargin);
    v = cell(m,n);
    for j=1:m
        v{j,1} = inputname(j);
        for i=1:length(varargin{j})
            if iscell(varargin{j}(i))
                v{j,i+1} = varargin{j}(i);
                v{j,i+1} = v{j,i+1}{:};
            else
                v{j,i+1} = varargin{j}(i);
            end
        end
    end
    try
        t = dataframe(v');
    catch
        try
            pkg load dataframe
        catch
            pkg install -forge dataframe
            pkg load dataframe
        end
        t = dataframe(v);
    end
else
    t = table(varargin{:});
    for i=1:length(varargin)
        t.Properties.VariableNames{i} = inputname(i);
    end
end
end
