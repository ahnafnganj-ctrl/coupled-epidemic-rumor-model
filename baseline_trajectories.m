function baseline_trajectories()

clc;
clear;
close all;

%% ================================================================
% BASELINE TRAJECTORIES (FIGURE 2 ONLY)
% Panels:
% (a) Baseline disease states
% (b) Baseline rumour states
% (c) Joint extinction
% (d) Rumour persistence
% (e) Coexistence
% ================================================================


%% Global formatting

set(groot,'defaultAxesFontName','Arial',...
    'defaultAxesFontSize',12,...
    'defaultAxesFontWeight','normal',...
    'defaultLineLineWidth',2.0,...
    'defaultTextFontName','Arial',...
    'defaultTextFontSize',12,...
    'defaultTextInterpreter','tex',...
    'defaultLegendInterpreter','tex',...
    'defaultAxesTickLabelInterpreter','tex');


%% Colour palette

co = [31 119 180;       % Blue
      214 39 40;        % Red
      44 160 44;        % Green
      148 103 189;      % Purple
      255 127 14]/255;  % Orange


%% Baseline parameters

p0 = struct('Lambda',0.02,...
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


%% Initial condition
% State order:
% S I V R U A C

y0 = [0.89;
      0.10;
      0.01;
      0.00;
      0.90;
      0.10;
      0.00];


%% Solver options

opts = odeset('RelTol',1e-9,...
              'AbsTol',1e-11,...
              'NonNegative',1:7);


%% Time

T = 10000;

tEval = unique([0,logspace(-2,log10(T),850)]);



%% ================================================================
% FIGURE 2
% ================================================================

fig = figure('Color','w',...
    'Position',[100 100 1016 837],...
    'Renderer','painters');



%% Baseline trajectory

[~,Y_base] = ode45( ...
    @(t,y) coupled_rhs(t,y,p0),...
    tEval,...
    y0,...
    opts);



%% (a) Baseline disease states

ax = subplot(2,6,1:3,'Parent',fig);

trajectory(ax,tEval,Y_base,...
    [1 2 3 4],...
    {'S','I','V','R'},...
    co(1:4,:),T);

title(ax,...
    '(a) Baseline disease states',...
    'FontName','Arial',...
    'FontSize',15,...
    'FontWeight','bold');

legend(ax,{'S','I','V','R'},...
    'Location','best',...
    'Orientation','horizontal',...
    'Box','off',...
    'FontName','Arial',...
    'FontSize',11);



%% (b) Baseline rumour states

ax = subplot(2,6,4:6,'Parent',fig);

trajectory(ax,tEval,Y_base,...
    [5 6 7],...
    {'U','A','C'},...
    co(1:3,:),T);

title(ax,...
    '(b) Baseline rumour states',...
    'FontName','Arial',...
    'FontSize',15,...
    'FontWeight','bold');

legend(ax,{'U','A','C'},...
    'Location','best',...
    'Orientation','horizontal',...
    'Box','off',...
    'FontName','Arial',...
    'FontSize',11);



%% Panels (c), (d), (e)

cases = {'extinction','rumour persistence','coexistence'};

titles_cases = { ...
    '(c) Joint extinction',...
    '(d) Rumour persistence',...
    '(e) Coexistence'};


p_cases = repmat(p0,1,3);


% Joint extinction
p_cases(1).beta = 0.4;
p_cases(1).lambda = 0.1;


% Rumour persistence
p_cases(2).beta = 0.4;


% Coexistence = baseline



for j = 1:3

    ax = subplot(2,6,6+(2*j-1:2*j),'Parent',fig);

    hold(ax,'on');


    [~,Y_case] = ode45( ...
        @(t,y) coupled_rhs(t,y,p_cases(j)),...
        tEval,...
        y0,...
        opts);


    mask = tEval>0 & tEval<=T;


    % Infection
    semilogx(ax,...
        tEval(mask),...
        Y_case(mask,2),...
        'Color',co(2,:),...
        'LineWidth',2);


    % Active rumour
    semilogx(ax,...
        tEval(mask),...
        Y_case(mask,6),...
        'Color',co(1,:),...
        'LineWidth',2);


    % Vaccination
    semilogx(ax,...
        tEval(mask),...
        Y_case(mask,3),...
        '--',...
        'Color',co(3,:),...
        'LineWidth',2);



    set(ax,'XScale','log');

    xlim(ax,[0.01 T]);

    ylim(ax,[0 1]);


    xlabel(ax,'Time','FontSize',14);

    ylabel(ax,'Active fraction','FontSize',14);


    title(ax,...
        titles_cases{j},...
        'FontName','Arial',...
        'FontSize',15,...
        'FontWeight','bold');


    legend(ax,{'I','A','V'},...
        'Location','best',...
        'Orientation','horizontal',...
        'Box','off',...
        'FontName','Arial',...
        'FontSize',10.5);


    clean_axis(ax);

end



disp('Baseline trajectory Figure 2 generated successfully.');

end



%% ================================================================
% MODEL
% ================================================================

function dy = coupled_rhs(~,y,p)

S=y(1);
I=y(2);
V=y(3);
R=y(4);

U=y(5);
A=y(6);
C=y(7);


vaccination = p.psi*(1-p.alpha*A)*S;

infection = p.beta*S*I;

breakthrough = p.sigma*p.beta*V*I;

activation = (p.lambda*A+p.phi*I)*U;


dy=[
p.Lambda-infection-vaccination+p.omega*V-p.mu*S;

infection+breakthrough-(p.gamma+p.mu)*I;

vaccination-breakthrough-(p.omega+p.mu)*V;

p.gamma*I-p.mu*R;

p.eta-activation+p.theta*C-p.mu*U;

activation-(p.delta+p.mu)*A;

p.delta*A-(p.theta+p.mu)*C
];

end



%% ================================================================
% TRAJECTORY FUNCTION
% ================================================================

function trajectory(ax,t,Y,states,labels,colors,T)

hold(ax,'on');

mask=t>0 & t<=T;


for k=1:length(states)

    semilogx(ax,...
        t(mask),...
        Y(mask,states(k)),...
        'Color',colors(k,:),...
        'LineWidth',2);

end


set(ax,'XScale','log');

xlim(ax,[0.01 T]);

ylim(ax,[0 1]);


xlabel(ax,'Time','FontSize',14);

ylabel(ax,'Population fraction','FontSize',14);


clean_axis(ax);

end



%% ================================================================
% AXIS FORMAT
% ================================================================

function clean_axis(ax)

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
    'LineWidth',1);


if strcmp(ax.XScale,'log')

    ticks=[1e-2 1e-1 1e0 1e1 1e2 1e3 1e4];

    ticks=ticks(ticks>=ax.XLim(1) & ticks<=ax.XLim(2));

    set(ax,...
        'XTick',ticks,...
        'XTickLabel',arrayfun(@(x) ...
        sprintf('10^{%d}',round(log10(x))),...
        ticks,...
        'UniformOutput',false));

end

end