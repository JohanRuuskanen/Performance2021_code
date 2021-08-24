%  moms = MarginalMomentsFromMMAP(D, K, precision)
%  
%  Returns the moments of the marginal distribution of a 
%  marked Markovian arrival process.
%  
%  Parameters
%  ----------
%  D : list/cell of matrices of shape(M,M), length(N)
%      The D0...DN matrices of the MMAP
%  K : int, optional
%      Number of moments to compute. If K=0, 2*M-1 moments
%      are computed. The default value is K=0.
%  precision : double, optional
%      Numerical precision for checking if the input is valid.
%      The default value is 1e-14
%  
%  Returns
%  -------
%  moms : row vector of doubles, length K
%      The vector of moments.

function moms = MarginalMomentsFromMMAP (D,K)

    if ~exist('K','var') || K==0
        K = 2*size(D{1},1)-1;
    end

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckMMAPRepresentation(D)
        error('MarginalMomentsFromMMAP: Input isn''t a valid MMAP representation');
    end

    [alpha,A] = MarginalDistributionFromMMAP(D);
    moms = MomentsFromPH(alpha,A,K);
end
