function  [data_new, ones_gaps_limit_all,ones_gaps_all]=spline_interp_find_gaps(data,gap_limit, t_old,t_new)

% [data_new, ones_gaps_limit_all,ones_gaps_all]=spline_interp_find_gaps(data,gap_limit, t_old,t_new)
% This function interpolates (spline) all the gaps in the data
%
% input
% data       : data you want to have interpolated
% gap_limit  : maximum gap size
% optional inputs
% t_old     : uneven spaced time vector belonging to the data, which may
% have with missing and or double samples
% t_new     : new, evenly spaced time vector where data should be resampled
% to. It is advisable to use this as:
% fs_intended = round(1/nanmedian(diff(t_old)));
% t_new=t_old(1):1/fs_intended:t_old(end);
% so that fs is the intended sample frequency of the data, which may not always have been
% met
%
% output
% data_new : interoploated data
% gap_matrix_limit  : Matrix with ones indicating all gaps exeeding the gap limit
% gap_matrix_all    : Matrix with ones indicating all gaps
reshapeflag=0;
n=[];
m=[];
if ndims(data)==2
    if size(data,2)>size(data,1)
        data=data';
        reshapeflag=1;
        disp('found more colums than rows in input data, assuming time is in colums')
    end
elseif ndims(data)==3
    [n,m,k]=size(data);
    data=reshape(data,n*m,k)';
    reshapeflag=2;
    disp('found 3d input data, assuming time is the third dimension')
end


if nargin ==4 % timeline corrections
    % remove nansamples (samples which have nan as timestamp
    data(isnan(t_old),:)    = [];
    t_old(isnan(t_old),:)   = [];

    % remove time samples that occur double for some reason, or jumps back in
    % time (process needs to be repeated until no such values are present
    % anymore
    t_diff_vlag=1;
    while t_diff_vlag==1
        data(diff(t_old)<=0,:)  = [];
        t_old(diff(t_old)<=0)   = [];
        if sum(diff(t_old)<=0)==0
            t_diff_vlag=0;
        end
    end    % interpolate to fixed timebase
    data   = interp1(t_old,data,t_new);
end %end of correction of time of data

% make zeros nan, because real zeros are hardly ever present in the data
data(data==0)=nan;

ones_gaps_all           = ones(size(data))*nan;
ones_gaps_limit_all     = ones(size(data))*nan;
data_new                = ones(size(data))*nan;

for i_col=1:size(data,2)
    data_col=data(:,i_col);
    ones_gaps=isnan(data_col);
    if sum(~isnan(data_col))>1  %more than 2 visible data points needed
        %% find indexes without gaps
        x_old_oke=find(~isnan(data_col));
        y_old_oke=data_col(x_old_oke);

        %% spline over all data_col
        x_new=1:(size(data_col,1)*size(data_col,2));%whole time serie
        y_new=interp1(x_old_oke,y_old_oke,x_new,'spline');

        % find start and end of gaps, and calculate gap size
        ones_gaps2=[0;ones_gaps;0];
        start_gaps=find((diff(ones_gaps2)==1));
        end_gaps=find((diff(ones_gaps2)==-1))-1;
        gap_size=(end_gaps-start_gaps)+1;

        ones_gaps_limit=ones_gaps*0;
        if ~isempty(gap_size)%
            %% make gap matrix with gaps larger than gap limit and with gaps at the
            %% start and the end
            large_gaps=[find(gap_size>gap_limit)]';%find gaps larger than gap_limit
            for i_gap=large_gaps% gaps larger than gap limit
                ones_gaps_limit(start_gaps(i_gap):end_gaps(i_gap))=1;
            end
            if start_gaps(1)==1%gaps at the start
                ones_gaps_limit(start_gaps(1):end_gaps(1))=1;
            end
            if end_gaps(end)==size(data_col,1)%gaps at the end
                ones_gaps_limit(start_gaps(end):end_gaps(end))=1;
            end
        end
        ones_gaps_all(:,i_col)=ones_gaps; %#ok<AGROW>
        ones_gaps_limit_all(:,i_col)=ones_gaps_limit; %#ok<AGROW>
        data_new(:,i_col)=y_new; %#ok<AGROW>
    else %% less than 2 visible data points found
        ones_gaps_all(:,i_col)=ones_gaps; %#ok<AGROW>
        ones_gaps_limit_all(:,i_col)=ones_gaps; %#ok<AGROW>
        data_new(:,i_col)=data_col; %#ok<AGROW>
    end
    keep ones_gaps_all ones_gaps_limit_all data_new data gap_limit reshapeflag ones_gaps_limit n m

end
if reshapeflag==1
    data_new=data_new';
    ones_gaps_all=ones_gaps_all';
    ones_gaps_limit_all=ones_gaps_limit';
elseif reshapeflag==2
    data_new=reshape(data_new',n,m,length(data_new));
    ones_gaps_all=reshape(ones_gaps_all',n,m,length(data_new));
    ones_gaps_limit_all=reshape(ones_gaps_limit_all',n,m,length(data_new));
end


