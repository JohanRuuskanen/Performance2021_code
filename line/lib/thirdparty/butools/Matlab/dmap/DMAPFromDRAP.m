%  [D0, D1] = DMAPFromDRAP(H0, H1, precision)
%  
%  Obtains a Markovian representation of a discrete rational
%  arrival process of the same size, if possible, using the
%  procedure published in [1]_.
%  
%  Parameters
%  ----------
%  H0 : matrix, shape (M,M)
%      The H0 matrix of the discrete rational arrival process
%  H1 : matrix, shape (M,M)
%      The H1 matrix of the discrete rational arrival process
%  precision : double, optional
%      A representation is considered to be a Markovian one
%      if it is closer to it than this precision
%  
%  Returns
%  -------
%  D0 : matrix, shape (M,M)
%      The D0 matrix of the discrete Markovian arrival process
%  D1 : matrix, shape (M,M)
%      The D1 matrix of the discrete Markovian arrival process
%  
%  References
%  ----------
%  .. [1] G Horvath, M Telek, "A minimal representation of 
%         Markov arrival processes and a moments matching 
%         method," Performance Evaluation 64:(9-12) pp. 
%         1153-1168. (2007)       

function [D0, D1] = DMAPFromDRAP (H0, H1, prec)

    if ~exist('prec','var')
        prec = 1e-14;
    end

    global BuToolsCheckInput;

    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckDRAPRepresentation(H0,H1)
        error('DMAPFromDRAP: Input isn''t a valid DRAP representation!');
    end

    H{1}=H0;
    H{2}=H1;
    Y=DMMAPFromDMRAP(H, prec);
    D0=Y{1};
    D1=Y{2};
end