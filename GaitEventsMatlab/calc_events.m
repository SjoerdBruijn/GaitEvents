function [events,signal]=calc_events(signal,fs,varargin)
% Function to calculate heel strikes and toe offs
% FUNCTION:
%       [events] = calc_events(signal,fs,tmpcor)
% INTPUT:
%       signal          =   Center of Pressure data in collumns (or traj
%       structure), if CoP; the columns should be: 1.ML(+ = right) the  2.AP(+ = forward)
%       fs              =   Sample frequency of your data
%       tmpcor          =   Correction for the estimation of the step length
% OUTPUT:
%       events          =   Structure containing the indices of the gait
%                           events of the CoP data
%
% Edit: Nick Kluft, March 2015
%       Changed: A` new approach to estimate the main frequency, Now ground
%       reaction forces not needed as input for this function. The main
%       frequency is now determined in the first 2 minutes of the data.
%
% ----V-----U------A-----M-----S-----T-----E-----r-----D-----A-----M-------
if nargin==2
    tmpcor = 1;
    fc = nan;
elseif nargin==3
    tmpcor = varargin{1};
    fc = nan;
elseif nargin ==4
    tmpcor = varargin{1};
    fc = varargin{2};
end


if isstruct(signal)
    traj = signal;
    sigistraj = true;
    signal=[];
    for i_seg = 1:size(traj.segment,2)
        if contains(lower(traj.segment(i_seg).name{1}),'foot')
            if ~isempty(traj.segment(i_seg).joint_name)
                iheel = find(~cellfun(@isempty,strfind(lower(traj.segment(i_seg).joint_name),'heel')),1,'first');
                if isempty(iheel)
                    iheel=3;
                end
                z = squeeze(traj.segment(i_seg).joint(3,iheel,:));
            else
                z = squeeze(traj.segment(i_seg).origin(3,1,:));
            end

            [z,tempz]   = spline_interp_find_gaps(z,30);

            if ~isnan(fc)
                [b,a] = butter(4, fc/(fs/2),'low');
                z     = filtfilt(b,a,z);
            end


            z(logical(tempz))=nan;

            %% differentiate

            zp = gradient(z);
            signal=[signal z zp];
            %% Using peakfind to get highest peaks in 2nd derivative of z
            [~,mDown] = peakfind(z,fs*tmpcor);
            [aUp,~] = peakfind(zp,fs*tmpcor);

            %% make structured output
            if contains(lower(traj.segment(i_seg).name{1}),'r')
                events.rto=aUp(:,1);
                events.rhs=mDown(:,1);
            else
                events.lto=aUp(:,1);
                events.lhs=mDown(:,1);
            end
        end
    end
else
    %% get some temporary heelstrikes, to get an idea of main frequency
    if ~isnan(fc)
        [b,a] = butter(4, fc/(fs/2),'low');
        signal= filtfilt(b,a,signal);
    end
    % Check the length of the data
    if length(signal(:,2))>120*fs
        nEst = (120*fs);
    else
        nEst = length(signal(:,2));
    end
    % Check power spectrum for first peak
    [pxx,f] = pwelch(signal(1:nEst,2),2048,25,0.1:.05:3,fs);
    try
        [M,~]   = peakfind(pxx,(.5));
        tmp     = (f(M(1))*fs)/2;
    catch
        rhs=find(diff(signal(1:nEst,2)>0)==1);
        tmp=(nanmean(diff(rhs)));
    end
    % If specified use the template correction
    tmp     = tmp*tmpcor;
    %% find heel contacts and toe contacts, using main period
    [maxout,minout]=peakfind(signal(:,2),tmp);
    hc=minout(:,1);
    to=maxout(:,1);
    %% create y-signal with no drift
    [B,A] = butter(4, 0.5/(fs/2), 'low');
    y=filtfilt(B,A,signal(:,1));
    y=signal(:,1)-y;

    %% split heel and toe to left and right
    events.lhs=hc(y(hc)> 0);
    events.rhs=hc(y(hc)< 0);
    events.lto=to(y(to)> 0);
    events.rto=to(y(to)< 0);
end
%% try ordering the evenst (may fail)
try
    events      = order_events(events);
catch
end
events.fs   = fs;
