function [events,notgood]=check_events(events, rawsignal, fs, varargin)
% Function to evaluate and edit calculated gait events
%   FUNCTION:
%       [events,notgood]=
%           check_signal(events, rawsignal, fs, ws, maxsubplots)
%   INPUT:
%       events (with fields)
%             .lhs     =   Indices of left heel strike
%             .rhs     =   Indices of right heel strike
%             .lto     =   Indices of left toe off
%             .rto     =   Indices of right toe off
%       signal  =   Collumn matrix Center of Pressure (CoP) data (x,y,z) or
%       fs      =   Sample Frequency of the signal
%       maxsubplots = Number of panels in figure
%
%   OUTPUT:
%       events (with fields)
%             .lhs     =   Indices of left heel strike
%             .rhs     =   Indices of right heel strike
%             .lto     =   Indices of left toe off
%             .rto     =   Indices of right toe off
%       notgood =   Gives a one when the data should not be saved
%
%   OTHER FUNCTIONS:
%       With the keyboard arrows you can scroll to the data. Using the up
%       and down arrow you can change the scope of the data presentation.
%       With a left button mouse press in the first panel you can display
%       the data surrounding the selected area.
%       Keyboard buttons a) and d) respectivelly adds the next on the most
%       nearby peak or removes the nearest heel strike or toe off. Using
%       the 'Calc Eve' button you can redo the events determination.
%       Specifying the cut-off frequency and a correction for the window
%       size used to calculate the events. When your done you can either
%       press 'done, exit' or 'exit, not good' to exit the figure. The
%       events can always be seen in the output of these function
%
%   Edit: Nick Kluft, 17 March 2015
%           added: Button Press in first panel, change scope with arrow
%           keys, scroll through the data using arrow keys, add/remove
%           points using keyboard buttons a) or d), redo determination of
%           events with different cut-off or other tmp.
%
%           removed: delete all points before or after ginput, Drawing
%           based on length heelstrike-> now based on total length signal.
%
% ----V-----U------A-----M-----S-----T-----E-----r-----D-----A-----M-------
[~,signal] = calc_events(rawsignal,fs);
%% try ordering the evenst (may fail)
try
    events      = order_events(events);
catch
end

ws = 10;%in seconds
maxwindow   = ws*fs;
if maxwindow > size(signal,1)
    ws = round((size(signal,1)/fs)/2);
    maxwindow   = ws*fs;
end

maxplots    = ceil(length(signal)/maxwindow)+1;
i_plot      = 1;
scrsz       = get(0,'ScreenSize');

notgood     = true;
if nargin>3
    maxsubplots = varargin{1};
else
    maxsubplots = 2;
end


warning off
figure('keyReleaseFcn',@KeyDownListener,'Tag','event_check');
set(gcf,'position',scrsz)
%% first set of buttons
uicontrol('Callback',@done_fnct, ...
    'String','Done, exit', ...
    'position',[5 470 140 30]);
uicontrol('Callback',@exit_fnct, ...
    'String','No good data',...
    'position',[5 440 140 30]);
%% second set of buttons
uicontrol('Style','text',...
    'String','fc =',...
    'FontSize',12,...
    'Position',[5 370 55 30]);
uicontrol('Style','text',...
    'String','win = % gait-cycle',...
    'FontSize',12,...
    'Position',[65 380 80 30]);

hfc = uicontrol('style','edit',...
    'String','10',...
    'position',[10 350 55 30]);
htmp= uicontrol('style','edit',...
    'String','1',...
    'position',[80 350 55 30]);
uicontrol('Callback',@det_events,...
    'String','Calc Events',...
    'position',[5 310 140 30]);
%% third set of buttons
uicontrol('Callback',@nxt_window_fnct, ...
    'String','OK, next window',...
    'position',[5 240 140 30]);
uicontrol('Callback',@prv_window_fnct, ...
    'String','One window back',...
    'position',[5 210 140 30]);
%% fourth set of buttons
uicontrol('Callback',@del_points_select_fnct, ...
    'String','Delete selected point(d)', ...
    'position',[5 130 140 30]);
uicontrol('Callback',@add_point_fnct2, ...
    'String','Add selected point(a)', ...
    'position',[5 100 140 30]);
%% keydownlistener
set(gcf,'keyReleaseFcn',@KeyDownListener)

%lhs==>rto==>rhs==>lto
%% buttons for direct control (alligned with x axis)
hlhs = uicontrol('Callback',@add_point_fnct, ...
    'Tag','lhs',...
    'String','Insert LHS',...
    'position',[0 30 100 30],...
    'ForegroundColor',  [1 0 1],...
    'style','togglebutton',...
    'value',0);
hlhs.BackgroundColor =[1 0 1];

uicontrol('Callback',@del_point_fnct, ...
    'Tag','lhs',...
    'String','Delete LHS ',...
    'position',[0 0 100 30],...
    'ForegroundColor',  [1 0 1],...
    'style','togglebutton',...
    'value',0);

uicontrol('Callback',@add_point_fnct, ...
    'Tag','rto',...
    'String','Insert RTO',...
    'position',[100 30 100 30],...
    'ForegroundColor',  [0 1 0],...
    'style','togglebutton',...
    'value',0);
%     hrto.BackgroundColor = [0 1 0];
uicontrol('Callback',@del_point_fnct, ...
    'Tag','rto',...
    'String','Delete RTO ',...
    'position',[100 0 100 30],...
    'ForegroundColor',  [0 1 0],...
    'style','togglebutton',...
    'value',0);

hrhs = uicontrol('Callback',@add_point_fnct, ...
    'Tag','rhs',...
    'String','Insert RHS',...
    'position',[200 30 100 30],...
    'ForegroundColor',  [0 0 0],...
    'style','togglebutton',...
    'value',0);
hrhs.BackgroundColor = [0 0 0];
uicontrol('Callback',@del_point_fnct, ...
    'Tag','rhs',...
    'String','Delete RHS ',...
    'position',[200 0 100 30],...
    'ForegroundColor',  [0 0 0],...
    'style','togglebutton',...
    'value',0);

hlto = uicontrol('Callback',@add_point_fnct, ...
    'Tag','lto',...
    'String','Insert LTO',...
    'position',[300 30 100 30],...
    'ForegroundColor',  [0 0 1],...
    'style','togglebutton',...
    'value',0);
hlto.BackgroundColor = [0 0 1];
uicontrol('Callback',@del_point_fnct, ...
    'Tag','lto',...
    'String','Delete LTO ',...
    'position',[300 0 100 30],...
    'ForegroundColor',  [0 0 1],...
    'style','togglebutton',...
    'value',0);

uicontrol('Style', 'text',...
    'String', 'No max found!',...
    'Units','normalized',...
    'Position', [0.9 0.7 0.1 0.1],...
    'ForegroundColor',[1 0 0],...
    'HandleVisibility','on',...
    'Visible','off');

%% draw plots
drawmain

%% Callback functions
    function drawmain
        plotpoints(events.lhs,events.lto,events.rhs,events.rto,signal,i_plot,maxsubplots,maxwindow);
        h_ax = findobj(gcf,'Type','axes');
        uph=zeros(1,length(h_ax));
        for k = 1:length(h_ax)
            uph(k) = h_ax(k).Position(2);
        end
        [~,upG] = max(uph);
        set(h_ax(upG),'ButtonDownFcn',@PressinAx);

        subplot(maxsubplots,1,1);
        if ~(length(events.lto)==length(events.rto)&&length(events.rhs)==length(events.rto)&&length(events.lhs)==length(events.rto))
            title(['Not all points are in the order lhs==>rto==>rhs==>lto, keep correcting: LTO: ',num2str(length(events.lto)),' LHS: ',num2str(length(events.lhs)),' RTO: ',num2str(length(events.rto)),' RHS: ',num2str(length(events.rhs))])
        elseif any((events.rto-events.lhs)<0)||any((events.rhs-events.rto)<0)||any((events.lto-events.rhs)<0)||any((events.lhs(2:end)-events.lto(1:end-1))<0)
            title(['Not all points are in the order lhs==>rto==>rhs==>lto, keep correcting: LTO: ',num2str(length(events.lto)),' LHS: ',num2str(length(events.lhs)),' RTO: ',num2str(length(events.rto)),' RHS: ',num2str(length(events.rhs))])
        else
            title(['All points are in the order lhs==>rto==>rhs==>lto :) : LTO: ',num2str(length(events.lto)),' LHS: ',num2str(length(events.lhs)),' RTO: ',num2str(length(events.rto)),' RHS: ',num2str(length(events.rhs))])
        end

    end

    function PressinAx(o,e)
        %function to proceed in windows
        i_plot = round(maxplots *(e.IntersectionPoint(1)/o.XLim(2)));
        gcf;
        drawmain
    end

    function det_events(~,~)
        % read events from text edit space
        fc = str2num(get(hfc,'String'));
        tmp = str2num(get(htmp,'String'));
        % redo the determination of gait events with new parameters

        events = [];
        [events,signal] = calc_events(rawsignal,fs,tmp,fc);
        clear cop
        drawmain
    end

    function nxt_window_fnct(~,~)
        %function to advance one window
        if i_plot<ceil(maxplots)
            i_plot=i_plot+1;
            drawmain
        end
    end

    function prv_window_fnct(~,~)
        %function to go to previous window within file
        if i_plot>1
            i_plot=i_plot-1;
            drawmain
        end
    end

    function KeyDownListener(handles,blah)
        switch blah.Key
            case 'd'
                del_points_select_fnct(handles,blah);
            case 'a'
                add_point_fnct2(handles,blah);
            case 'rightarrow'
                nxt_window_fnct(handles,blah)
            case 'leftarrow'
                prv_window_fnct(handles,blah)
            case 'uparrow'
                ratP = i_plot/maxplots;
                if ws< floor(length(signal)/2)-1
                    ws = ws+1;
                    maxwindow   = ws*fs;
                    maxplots    = ceil(length(signal)/maxwindow)+1;
                end
                i_plot = ratP*maxplots;
                drawmain
            case 'downarrow'
                ratP = i_plot/maxplots;
                if ws>1
                    ws = ws-1;
                    maxwindow   = ws*fs;
                    maxplots    = ceil(length(signal)/maxwindow)+1;
                end
                i_plot = ratP*maxplots;
                drawmain
        end
    end

    function add_point_fnct(h,~)
        %function to add points
        variable=get(h,'Tag');
        [x,~] = ginput(1);
        new =[events.(variable);round(x)];
        events.(variable) = sort(new);
        clear x new
        set(h,'Value',0)
        drawmain
    end


    function add_point_fnct2(~,~)
        % get users input
        [x,y] = ginput(1);
        % find closest event-point to input x
        [~,e1(1)] = mindist(x,events.lhs);
        [~,e1(2)] = mindist(x,events.rto);
        [~,e1(3)] = mindist(x,events.rhs);
        [~,e1(4)] = mindist(x,events.lto);
        [val_e,e2] = mindist(x,[events.lhs(e1(1));events.rto(e1(2));events.rhs(e1(3));events.lto(e1(4))]);
        if x-val_e > 0
            vars ={'rto','rhs','lto','lhs'}; % if
        else
            vars ={'lto','lhs','rto','rhs'};
        end
        % point to add= variable
        variable = vars{e2};
        % round x to calculate minimal distance to top
        x = round(x);
        % find signal closest to selected point
        [~,sign] = mindist(y,signal(x,:));
        % let the peak find window size be a function of part of the signal
        % shown
        R = round(maxwindow/20);
        % don't survey comlete signal only 1 second surrounding selected
        % point
        range       = (x-(R)):(x+(R));
        % set boundaries to range
        range(range<=0) =[];range(range>length(signal))=[];
        % select the window for the peakfind
        Peak_ws     = round((.1*length(range)));
        % find closest peak to the selected point
        [Ma,Mi] = peakfind(signal(range,sign),Peak_ws);
        if ~isempty(Ma) && ~isempty(Mi)
            % if a Maxima and Minima is found check for nearest peak,
            % mindist finds index closest. Translate it with range to the
            % index relative to the signal
            [new,~] = mindist(x,[range(Ma(:,1)),range(Mi(:,1))]);
        elseif isempty(Ma)
            new = range(min(Mi(:,1)));
        elseif isempty(Mi)
            new = range(max(Ma(:,1)));
        else
            % if none is found use the selected point as new point
            new = x;
        end
        % add new point to the existing points
        events.(variable) = sort([events.(variable);new]);
        clear x y new vec win index
        drawmain
    end

    function del_point_fnct(h,~)
        %function to delete points
        variable=get(h,'Tag');
        [x,~]=ginput(1);
        x=sort(x);
        ind = find(x<events.(variable));
        events.(variable)(ind(1))=[];
        clear x y ind
        set(h,'Value',0)
        drawmain
    end

    function del_points_select_fnct(~,~)
        %function to delete nearest point
        [x,~]=ginput(1);
        [a,i(1)] = min(abs(events.lto-x));
        [b,i(2)] = min(abs(events.rto-x));
        [c,i(3)] = min(abs(events.lhs-x));
        [d,i(4)] = min(abs(events.rhs-x));
        [~,e]     = min([a;b;c;d]);
        variable= {'lto','rto','lhs','rhs'};
        events.(variable{e})(i(e))=[];
        clear x y ind a b c d e i variable
        drawmain
    end

    function exit_fnct(~,~)
        notgood=true;
        events.tol=[];
        events.tor=[];
        events.hcl=[];
        events.hcr=[];
        uiresume
    end

    function done_fnct(~,~)
        notgood = false;
        i_plot  = ceil(maxplots);
        uiresume
    end

uiwait
figuur = findobj('Type','figure','Tag','event_check');
if ishandle(figuur)
    close(figuur)
end

end






