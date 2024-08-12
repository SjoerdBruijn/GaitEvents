function M = standing(M)
% function makes sure the vector or matrix is a standing vector
% (a vector/matrix of which the n_rows > n_cols).
%
% FUNCTION: 
%           M = standing(M)
%    INPUT:
%           M : any matrix of vector
%   OUTPUT: 
%           M : standing vector/matrix;
% 
% ----V-----U------A-----M-----S-----T-----E-----r-----D-----A-----M-------

if ~iscolumn(M)
    M = transpose(M);
end

end
