function SS = spaceCache(n, m)
% SS = SPACECACHE(N, M)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

SS = [];
HS = nchoosek(1:sum(n),sum(m)); % hotset
for i=1:size(HS,1)
    SS = [SS; perms(HS(i,:))];
end

end
