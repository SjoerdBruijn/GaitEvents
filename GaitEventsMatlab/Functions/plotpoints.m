function      plotpoints(lhs,lto,rhs,rto,signal,i_plot,maxsubplots,maxwindow)
subplot(maxsubplots,1,1);
% make axis of first plot as long as the signal for easy control
xas = [1:(length(signal)/length(diff(lhs))):length(signal)];
% create plots of the first panel
plot([1:(length(signal)/length(diff(lto))):length(signal)],diff(lto),'Color',[0 0 1]);hold on
plot([1:(length(signal)/length(diff(lhs))):length(signal)],diff(lhs),'Color',[1 0 1])
plot([1:(length(signal)/length(diff(rto))):length(signal)],diff(rto),'Color',[0 1 0])
plot([1:(length(signal)/length(diff(rhs))):length(signal)],diff(rhs),'Color',[0 0 0])
% draw to lines to see whether it's under the std
line([xas(1) xas(end)],[nanmean(diff(lhs))-std(diff(lhs)) nanmean(diff(lhs))-std(diff(lhs))])
line([xas(1) xas(end)],[nanmean(diff(lhs))+std(diff(lhs)) nanmean(diff(lhs))+std(diff(lhs))])

Meve = [[lhs;rto;rhs;lto],...
    [zeros(size(lhs));ones(size(rto));...
    ones(size(rhs))*2;ones(size(lto))*3]];
idxwrong=[];
if length(Meve)>32
    % sort now
    Meve = sortrows(Meve,1);
    Morder = [rem(1:length(Meve),4);rem(2:length(Meve)+1,4);...
        rem(3:length(Meve)+2,4);rem(4:length(Meve)+3,4)]';
    dMeve = Morder(:,:)-Meve(:,2);
    [~,eveordercol] = min(sum(dMeve~=0));
    ieve_wrong = find(abs(Morder(:,eveordercol)-Meve(:,2))>0,1,'first');
    idxwrong = Meve(ieve_wrong,1);
    evewrong = Meve(ieve_wrong,2)+1;
    vline(idxwrong,':r',[],[],[],2);
end


lijnkleur = ['k','b','k','b','c','k'];

evenames = {'lhs','rto','rhs','lto'};
% get domain you are interested in
domein = [round((i_plot*maxwindow)-(maxwindow-1)):round(i_plot*maxwindow)];
% When the selected area does not contain any data points, don't go any
% further:
if domein(1)<0
    domein = (domein+abs(domein(1)))+1;
end
if domein(end)> length(signal)
    domein = domein-(domein(end)-length(signal));
end
% set axis in a nice manner
axis tight
set(gca,'YLim',[nanmean(diff(lhs))-500 nanmean(diff(lhs))+500])
h=get(gca);
% draw lines: borders of domain
line([domein(1) domein(1)],[h.YLim(1) h.YLim(2)],'LineWidth',2);
line([domein(end) domein(end)],[h.YLim(1) h.YLim(2)],'LineWidth',2);hold off
% draw the subplots
for j=1:maxsubplots-1
    % divide domein in the number of subplots
    eval(['domein',num2str(j),' = domein(round(length(domein)*(1/(maxsubplots-1))*(j-1)+1):'....,
        'round(length(domein)*(j/(maxsubplots-1))));'])
    for en = 1:length(evenames)
        % Find indices which fall inside boundries given by domein
        eval(['i_',evenames{en},'=find((',evenames{en},' > domein',...
            num2str(j),'(1)) & (',evenames{en},' < domein',num2str(j),'(end)));']);
    end
    % create first subplot
    subplot(maxsubplots,1,j+1);
    hold off;
    for k = 1:size(signal,2)
        % plot the signal input:
        eval(['hl=plot((domein',num2str(j),'),signal(domein',num2str(j),',k));']);
        % Change colors of the lines
        set(hl,'Color',lijnkleur(k));
        hold on;
    end
    % plot the events on the curves
    subplot(maxsubplots,1,j+1); plot(lto(i_lto),signal(lto(i_lto),:),'o','Color',[0 0 1]);
    subplot(maxsubplots,1,j+1); plot(lhs(i_lhs),signal(lhs(i_lhs),:), 'o','Color',[1 0 1]);
    subplot(maxsubplots,1,j+1); plot(rto(i_rto),signal(rto(i_rto),:), 'o','Color',[0 1 0]);
    subplot(maxsubplots,1,j+1); plot(rhs(i_rhs),signal(rhs(i_rhs),:), 'o','Color',[0 0 0]);
    if ~isempty(idxwrong)
       if and(idxwrong>domein(1),idxwrong<domein(end))
            subplot(maxsubplots,1,j+1); 
            plot(idxwrong,signal(idxwrong,:),...
                'o','Color',[1 0 0],'MarkerSize',10);
       end
    end
    hold off
    title('points check correct using options left')
    xlim([eval(['domein',num2str(j),'(1)']) eval(['domein',num2str(j),'(end)'])])
end
