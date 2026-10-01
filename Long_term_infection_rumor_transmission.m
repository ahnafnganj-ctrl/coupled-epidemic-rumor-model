function Long_term_infection_rumor_transmission()
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

% Figures 4 and 5 only

% Figure 4.
f=newfig([1193 825]); ks={'psi','alpha'}; sv={[0.2 0.5 0.8],[0.05 0.3 0.9]};
for row=1:2
    for j=1:3
        ax=subplot(2,3,(row-1)*3+j,'Parent',f); hold(ax,'on');
        d=readcsv(dataDir,sprintf('response_%s_%d',ks{row},j-1));
        vals=unique(d.(ks{row})); labs=cell(1,numel(vals));
        for k=1:numel(vals)
            mask=d.(ks{row})==vals(k);
            plot(ax,d.beta(mask),d.I(mask),'Color',co(k,:));

            % Analytical invasion thresholds at E0 and, when it exists, E1.
            p4=baseline_parameters();
            p4.sigma=sv{row}(j); p4.(ks{row})=vals(k);
            Nbar4=p4.Lambda/p4.mu; Mbar4=p4.eta/p4.mu;
            beta0=(p4.gamma+p4.mu)*(p4.psi+p4.omega+p4.mu)/...
                (Nbar4*(p4.omega+p4.mu+p4.sigma*p4.psi));
            xline(ax,beta0,':','Color',co(k,:),'LineWidth',1.25,...
                'HandleVisibility','off');
            if p4.lambda*Mbar4 > p4.delta+p4.mu
                A1=(Mbar4-(p4.delta+p4.mu)/p4.lambda)/...
                    (1+p4.delta/(p4.theta+p4.mu));
                q1=p4.psi*(1-p4.alpha*A1);
                beta1=(p4.gamma+p4.mu)*(q1+p4.omega+p4.mu)/...
                    (Nbar4*(p4.omega+p4.mu+p4.sigma*q1));
                xline(ax,beta1,'--','Color',co(k,:),'LineWidth',1.25,...
                    'HandleVisibility','off');
            end
            if strcmp(ks{row},'psi')
                labs{k}=sprintf('\\psi = %g',vals(k));
            else
                labs{k}=sprintf('\\alpha = %g',vals(k));
            end
        end
        xlabel(ax,'\beta','FontSize',15,'Interpreter','tex');
        ylabel(ax,'\itI(T)','FontSize',15,'Interpreter','tex');
        axis(ax,[0 3 0 0.07]); clean(ax);
        paneltit(ax,char(96+(row-1)*3+j),sprintf('\\sigma = %g',sv{row}(j)),15);
        if j==1
            legend(ax,labs,'Box','off','FontName','Arial','FontSize',10,...
                'Location','best');
        end
    end
end
annotation(f,'textbox',[0.08 0.005 0.84 0.035],...
    'String',['Threshold markers: dotted R_{0E} = 1; ',...
    'dashed R_{E|1} = 1 (colours match curves).'],...
    'EdgeColor','none','HorizontalAlignment','center',...
    'Interpreter','tex','FontName','Arial','FontSize',10);
showfig(f);

% Figure 5.
f=newfig([1097 793]); bv=[0.2 0.7 1.4 3];

for j=1:4
    ax=subplot(2,2,j,'Parent',f); hold(ax,'on');

    d=readcsv(dataDir,sprintf('response_phi_%d',j-1));
    vals=unique(d.phi); labs=cell(1,4);

    for k=1:4
        mask=d.phi==vals(k);

        plot(ax,d.lambda(mask),d.A(mask),...
            'Color',co(k,:));

        labs{k}=sprintf('\\phi = %g',vals(k));
    end

    xlabel(ax,'\lambda','FontSize',15,'Interpreter','tex');
    ylabel(ax,'\itA(T)','FontSize',15,'Interpreter','tex');

    xlim(ax,[0 2]);
    clean(ax);

    paneltit(ax,char(96+j),sprintf('\\beta = %g',bv(j)),15);

    if j==1
        legend(ax,labs,...
            'Box','off',...
            'FontName','Arial',...
            'FontSize',10,...
            'Location','best',...
            'NumColumns',2);
    end
end

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


function d = readcsv(dataDir,name)

d = readtable(fullfile(dataDir,[name,'.csv']));

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


function heat(ax,name,state,dataDir)

d=readcsv(dataDir,name);

vars=d.Properties.VariableNames;

xv=unique(d.(vars{1}));
yv=unique(d.(vars{2}));

zz=reshape(d.(state),numel(xv),numel(yv))';

contourf(ax,xv,yv,zz,100,'LineStyle','none');

colormap(ax,parula);

set(ax,'YDir','normal',...
    'Box','on',...
    'TickDir','in',...
    'Layer','top',...
    'LineWidth',1,...
    'FontName','Arial',...
    'FontSize',11.5);

xlabel(ax,vars{1},'FontSize',15);
ylabel(ax,vars{2},'FontSize',15);

cb=colorbar(ax);
cb.FontName='Arial';

axis(ax,'square');

end

function p = baseline_parameters()

p = struct('Lambda',0.02,...
    'eta',0.02,...
    'mu',0.02,...
    'beta',0.70,...
    'gamma',0.25,...
    'psi',0.30,...
    'alpha',0.50,...
    'omega',0.05,...
    'sigma',0.40,...
    'lambda',0.25,...
    'phi',0.60,...
    'delta',0.15,...
    'theta',0.10);

end