function MMAP = mmap_hide(MMAP, types)
% MMAP_HIDE. Makes a subset of the MMAP type hidden.
% MMAP = MMAP_HIDE(MMAP, TYPES) takes a MMAP with K arrival types and
% returns a new MMAP for the same stochastic process with those arrivals
% hidden from observation.

D0 = MMAP{1};
D1 = MMAP{2};
for k=types(:)'
    D0 = D0 + MMAP{2+k};
    D1 = D1 + MMAP{2+k};
    MMAP{2+k} = 0 * MMAP{2+k};
end
MMAP = mmap_normalize(MMAP);
end