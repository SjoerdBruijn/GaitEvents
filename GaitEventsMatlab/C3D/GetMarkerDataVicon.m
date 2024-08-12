function [data]=GetMarkerData(POINTdat,ParameterGroup,labels)
% get markerdata from POINTdat, returned from readC3D. 
% Example usage: [data]=GetMarkerData(POINTdat,ParameterGroup,{'RANK','LANK'});
% April 2012, SMB

[n,m,k]=size(POINTdat);
[n1,m1]=size(labels);
data=nan*ones(3,m1,k);
parameternames = [ParameterGroup(strcmp([ParameterGroup.name],'POINT')).Parameter.name];
pointlabels = ParameterGroup(strcmp([ParameterGroup.name],'POINT')).Parameter(strcmp(parameternames,'LABELS')).data;

for i=1:length(pointlabels)
    tmp=pointlabels{i};
    if ~isempty(tmp)
        tmp=strsplit(tmp,':');
        if size(tmp,2)>1
            pointlabels(i)=tmp(2);
        else
            pointlabels(i)=tmp;
        end
    else
        pointlabels{i}=tmp;
    end
end
%% RB Test to fix. remove subj/platform name from 'labels'
for i=1:length(labels)
    tmp=labels{i};
    if ~isempty(tmp)
        tmp=strsplit(tmp,':');
        if size(tmp,2)>1
            labels(i)=tmp(2);
        else
            labels(i)=tmp;
        end
    else
        labels{i}=tmp;
    end
end
%
if iscellstr(labels)
    for i=1:m1
        if sum(strcmp(pointlabels,char(labels(i)))) == 1
            data(1:3,i,:)= POINTdat(:,strcmp(pointlabels,char(labels(i))),:);
        else
            disp(['Label :' char(labels(i)) ' not found, returning without data (NaNs instead) for this marker'])
        end
            
    end
else
    disp('input variable "labels" should be a cell containing strings. Returning without results')
end