%  re = EmpiricalRelativeEntropy(f1, f2, intBounds)
%  
%  Returns the relative entropy (aka Kullback–Leibler
%  divergence) of two continuous functions given by samples
%  and the bounds of the corresponding intervalls.
%  
%  This function can be used to characterize the distance
%  between two density functions, distribution functions, 
%  etc.
%  
%  Parameters
%  ----------
%  f1 : vector, length M
%      Samples of the first continuous function
%  f2 : vector, length M
%      Samples of the second continuous function
%  intBounds : vector, length M+1
%      The bounds of the intervals. The ith sample
%      corresponds to the interval 
%      (intbounds(i),intbounds(i+1))
%  
%  Returns
%  -------
%  re : double
%      The relative entropy

function re = EmpiricalRelativeEntropy (f1, f2, intBounds)

    intlens = intBounds(2:end) - intBounds(1:end-1);
    intlens = reshape(intlens, size(f1));
    p1 = f1 .* intlens;
    p2 = f2 .* intlens;
    re = RelativeEntropy (p1, p2);
end

