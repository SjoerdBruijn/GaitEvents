function [maxout, minout] = peakfind(x,ws)
%
% function [maxout, minout] = peakfind(x,ws)
%
% Finds local maxima and minima of the vector x in a window of size ws (in
% samples).
%
% Any extremum is a local maximum (or minimum) in the sense that the preceeding ws/2 samples and
% the following ws/2 samples had lower (or higher) values.
%
% The minimum window size is 3, even window sizes are increased by 1.
% The speed of the algorithm is inversely proportional to ws 
%   
% Data processing routine
% Arne Ridderikhoff 9/22/2002

% x is a column vector
x = x(:);

% if the window size is not odd, add one sample
ws = max(3,2*fix(ws/2)+1);

% the number of bins
nbins = ceil(length(x)/ws);

% the centers of each bin are evaluated, all centers are elements on one
% row in the matrix y (see below)
center = ceil(ws/2);

% padding parameters
npad = fix(ws/2);

% locations of extrema
locmax = [];
locmin = [];

% make a column vector and pad it with nan
% every element k in the vector z corresponds to the element k-npad in the
% vector x
z= [NaN*ones(npad,1); x; NaN*ones(2*ws,1)];

% the window shift
for i = 1:ws
    % the vector for analysis
    % every element k in the vector y corresponds to the element k + i - 1
    % in the vector z, and hence to the element k + i - npad - 1 in x
    y = z(i:ws*nbins+i-1);
    
    % reshape y in [nbins] columns of length [ws] (nbins x ws matrix)
    % the key idea is that the center of each bin is going to be evaluated
    % relative to its neighbours
    % ==> this exploits the matrix calculation of MATLAB
    y = reshape(y,ws,nbins);
    
    % evaluate center of each bin, and if the center is a local extremum,
    % keep its index (in the original vector y)
    % if the center is not a local extremum, ignore it
    locmax_in_y = find(max(y) == y(center,:)).*ws - center + 1;
    locmin_in_y = find(min(y) == y(center,:)).*ws - center + 1;    
      
    % transform the indices of the extrema in y to indices of the extrema
    % in x
    locmax = [locmax locmax_in_y + i - npad - 1];
    locmin = [locmin locmin_in_y + i - npad - 1];
end    

% remove 1st and last sample of timeseries (artificial extremes)
locmax(find(locmax== 1 | locmax == length(x))) = [];
locmin(find(locmin== 1 | locmin == length(x))) = [];

% sort data and create output matrices
locmax = sort(locmax);
locmin = sort(locmin);
maxout = [locmax(:) x(locmax)];
minout = [locmin(:) x(locmin)];

% check (if no output is specified)
if nargout == 0
    plot(x),hold on
    plot(maxout(:,1),maxout(:,2),'ro')
    plot(minout(:,1),minout(:,2),'r+')
end    



