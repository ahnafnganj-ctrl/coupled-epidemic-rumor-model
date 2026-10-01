function Case_studies()
clc; clear; close all;

dataDir = fullfile(pwd,'matlab_output','data');

co = [ 31 119 180;
       214 39 40;
        44 160 44;
       148 103 189;
       255 127 14]/255;

set(groot,'defaultAxesFontName','Arial','defaultAxesFontSize',12,...
    'defaultAxesFontWeight','normal','defaultLineLineWidth',2.0,...
    'defaultTextFontName','Arial','defaultTextFontSize',12,...
    'defaultTextInterpreter','tex','defaultLegendInterpreter','tex',...
    'defaultAxesTickLabelInterpreter','tex');

% Figure 3 only

% Figure 3.
f=newfig([1193 481]);
cases={'baseline','low_rumour','high_rumour','high_vaccination','high_leakiness'};
labels={'Baseline','\lambda = 0.05','\lambda = 0.9','\psi = 0.7','\sigma = 0.7'};
sty={'-','--','-.',':','-'}; states={'I','A','V'};
titles={'(a) Infection','(b) Active rumour','(c) Vaccination'};

for j=1:3
    ax=subplot(1,3,j,'Parent',f);
    hold(ax,'on');

    for k=1:5
        lw=2.0;
        if k==4
            lw=2.5;
        end

        d=readcsv(dataDir,['trajectory_',cases{k}]);

        % Remove t = 0 because logarithmic axes cannot display zero
        mask = d.t > 0 & d.t <= 1000;

        semilogx(ax,d.t(mask),d.(states{j})(mask),...
            'Color',co(k,:),...
            'LineStyle',sty{k},...
            'LineWidth',lw);
    end

    % Explicitly use logarithmic x-axis
    set(ax,'XScale','log');

    xlabel(ax,'Time','FontSize',15);
    ylabel(ax,['\it',states{j},'(t)'],...
        'FontSize',15,'Interpreter','tex');

    % Same lower limit as Figure 2
    xlim(ax,[0.01 1000]);

    clean(ax);
    tit(ax,titles{j},15);

    if j==1
        legend(ax,labels,...
            'Box','off',...
            'FontName','Arial',...
            'FontSize',10.5,...
            'Location','best');
    end
end
showfig(f);



end

%% ================= HELPER FUNCTIONS =================

function f = newfig(sz)

scr = get(groot,'ScreenSize');

left = max(40, round((scr(3)-sz(1))/2));
bottom = max(40, round((scr(4)-sz(2))/2));

f = figure( ...
    'Color','w',...
    'Units','pixels',...
    'Position',[left bottom sz(1) sz(2)],...
    'Visible','on',...
    'Renderer','painters');

end


function showfig(fh)

figure(fh);
drawnow;

end


function clean(ax)

grid(ax,'on');

set(ax,...
    'FontName','Arial',...
    'FontSize',12,...
    'FontWeight','normal',...
    'TickLabelInterpreter','tex',...
    'Box','off',...
    'TickDir','out',...
    'GridColor',[0.8 0.8 0.8],...
    'GridAlpha',0.5,...
    'LineWidth',1.0);

end


function tit(ax,s,fs)

title(ax,s,...
    'FontName','Arial',...
    'FontSize',fs,...
    'FontWeight','bold',...
    'Interpreter','tex');

end


function paneltit(ax,letter,rest,fs)

title(ax,...
    ['\bf(',letter,')\rm ',rest],...
    'FontName','Arial',...
    'FontSize',fs,...
    'FontWeight','normal',...
    'Interpreter','tex');

end


function d = readcsv(dataDir,name)

d = readtable(fullfile(dataDir,[name,'.csv']));

end