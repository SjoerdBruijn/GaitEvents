function [max_value,index] = mindist(target,dat)
% find minimal dististance to target
% returns max_value and the index of the point closest to target
% function 
%           mindist(target,dat)
% 
% input     target: Target value from where distance in computed(scalar or standing vector)
%           dat   : The data from which the distance is measured (matrix, calculates min distance by collumns).
% output    max_value: value of data point closest to target
%           index: index of closest point in dat structure (per collumn)
% Nick Kluft

target = repmat(target,size(dat));
[~,index] = min(abs(target-dat));
if isequal(sum(ismember(size(dat),1)),0)
max_value = diag(dat(index,1:size(dat,2)))';
else
    max_value = dat(index);
end
