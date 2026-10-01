function Heatmap_analysis()
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

% Figures 6-10 only

% Figure 6.
f=newfig([1144 473]);
ax1=subplot(1,2,1,'Parent',f);
heat(ax1,'map_psi_sigma','I','(a) Vaccination and leakiness',[],...
    '\psi','\sigma','\itI(T)');
ax2=subplot(1,2,2,'Parent',f);
heat(ax2,'map_beta_sigma','I','(b) Transmission and leakiness',[],...
    '\beta','\sigma','\itI(T)');
showfig(f);

% Figure 7.
f=newfig([712 489]);
heat(axes('Parent',f),'map_psi_delta','V','Vaccination and rumour cessation',[],...
    '\psi','\delta','\itV(T)');
showfig(f);

% Figure 8.
f=newfig([1144 473]);
ax1=subplot(1,2,1,'Parent',f);
heat(ax1,'map_alpha_psi','V','(a) Rumour effect and uptake',[],...
    '\alpha','\psi','\itV(T)');
ax2=subplot(1,2,2,'Parent',f);
heat(ax2,'map_alpha_sigma','V','(b) Rumour effect and leakiness',[],...
    '\alpha','\sigma','\itV(T)');
showfig(f);

% Figure 9.
f=newfig([1224 809]);
maps={'map_alpha_phi','map_sigma_phi'}; states={'I','V','A'};
xlabs={'\alpha','\sigma'}; panelLetters={'a','b','c','d','e','f'};
for row=1:2
    for j=1:3
        k=(row-1)*3+j; ax=subplot(2,3,k,'Parent',f);
        cblab=['\it',states{j},'(T)'];
        heat(ax,maps{row},states{j},'',[],xlabs{row},'\phi',cblab);
        paneltit(ax,panelLetters{k},cblab,15);
    end
end
showfig(f);

% Figure 10.
f=newfig([712 489]);
heat(axes('Parent',f),'map_psi_phi','A',...
    'Vaccination and infection-induced activation',[],...
    '\psi','\phi','\itA(T)');
showfig(f);



end


%% ================= HELPER FUNCTIONS =================

function f = newfig(sz)

scr=get(groot,'ScreenSize');
left=max(40,round((scr(3)-sz(1))/2));
bottom=max(40,round((scr(4)-sz(2))/2));

f=figure('Color','w',...
    'Units','pixels',...
    'Position',[left bottom sz(1) sz(2)],...
    'Visible','on',...
    'Renderer','painters');

end


function showfig(fh)
figure(fh);
drawnow;
end


function d=readcsv(dataDir,name)
d=readtable(fullfile(dataDir,[name,'.csv']));
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


function paneltit(ax,letter,rest,fs)

title(ax,['\bf(',letter,')\rm ',rest],...
    'FontName','Arial',...
    'FontSize',fs,...
    'FontWeight','normal',...
    'Interpreter','tex');

end


function heat(ax,name,state,titleText,limits,xlab,ylab,cblab)

dataDir = fullfile(pwd,'matlab_output','data');

d = readtable(fullfile(dataDir,[name,'.csv']));

vars = d.Properties.VariableNames;

xv = unique(d.(vars{1}));
yv = unique(d.(vars{2}));

zz = reshape(d.(state),numel(xv),numel(yv))';


contourf(ax,xv,yv,zz,100,'LineStyle','none');

colormap(ax,parula);


if ~isempty(limits)
    caxis(ax,limits);
else
    caxis(ax,'auto');
end


set(ax,...
    'YDir','normal',...
    'Box','on',...
    'TickDir','in',...
    'Layer','top',...
    'LineWidth',1.0,...
    'FontName','Arial',...
    'FontSize',11.5);


xlim(ax,[min(xv) max(xv)]);
ylim(ax,[min(yv) max(yv)]);

axis(ax,'square');


xlabel(ax,xlab,...
    'Interpreter','tex',...
    'FontName','Arial',...
    'FontSize',15);

ylabel(ax,ylab,...
    'Interpreter','tex',...
    'FontName','Arial',...
    'FontSize',15);


cb=colorbar(ax);

cb.FontName='Arial';
cb.FontSize=11.5;

cb.Label.String=cblab;
cb.Label.Interpreter='tex';
cb.Label.FontName='Arial';
cb.Label.FontSize=15;


if ~isempty(titleText)

    title(ax,titleText,...
        'FontName','Arial',...
        'FontSize',15,...
        'FontWeight','bold',...
        'Interpreter','tex');

end

end
