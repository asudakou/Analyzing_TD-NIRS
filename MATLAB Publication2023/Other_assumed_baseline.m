% Written by Aleh Sudakou (https://www.researchgate.net/profile/Aleh-Sudakou)
% Nałęcz Institute of Biocybernetics and Biomedical Engineering, Polish Academy of Sciences 
% This is one of the codes shared on  https://github.com/asudakou/Analyzing_TD-NIRS
% Last updated: 20 May 2023


%% Description:
% Generate two DTOFs (before and after change in absorption) and retreive 
% changes in optical properties using varying values of baseline optical properties
% 
% The results are presented in Fig. 3 in the 2023 Publication.


%% Main

clear; % clc

cut_lim_mom = [0.25 0.03];

rho = 30;

n = 1.33;

L_known = 15; 

opt_prop_base = [0.01 0.01 0.01 1 1 1 L_known 40]; % [Mua1 Mua2 Mua3 Musp1 Musp2 Musp3 L1 L2]

opt_prop_after = [0.01 0.015 0.015 1 1 1 L_known 40]; % Incrase second layer's Mua

[R_base, ~]        = DTOF_generate_Liemert(opt_prop_base,  n, rho, -1);
[R_after, time_ns] = DTOF_generate_Liemert(opt_prop_after, n, rho, -1);

[ R_base, cut_ind_mom ] = DTOF_filter( R_base, time_ns, 0, 0, cut_lim_mom, 1);
[ R_after, ~ ] = DTOF_filter( R_after, time_ns, 0, 0, cut_ind_mom, 3);

scal_fact = 10^6 / sum(R_base); % This step makes no difference to calculations, but good for plotting DTOFs
R_base = R_base * scal_fact;
R_after = R_after * scal_fact;

Mom_base = DTOF_CentralMom(time_ns, R_base);
Mom_after = DTOF_CentralMom(time_ns, R_after);

delMom = DTOF_DelMom(Mom_base, Mom_after);


% Define ranges of varying baseline parameters:
ranges = [];
ranges.L = 1:0.5:20;

ranges.musp1 = 0.7:0.05:1.3;

ranges.musp2 = ranges.musp1;

ranges.mua1 = 0.005:0.0005:0.015;

ranges.mua2 = ranges.mua1;

ranges.n1 = sort(1.13:0.025:1.53);

ranges.n2 = ranges.n1;

Results = '';


%% L
OptProp_base = opt_prop_base(1:6); % [Mua1 Mua2 Mua3 Musp1 Musp2 Musp3]

delOptProp_all = [-7 -7 -8 0 0 0 L_known 40]; % [delMua1 delMua2 delMua3 delMusp1 delMusp2 delMusp3 L1 L2]

for j = 1:length(ranges.L)
    delOptProp_all(7) = ranges.L(j); % Vary the value of L1
    
    % Use LMA #2 to retrieve delMua1 and delMua2
    [X_found, ~] = LMA_2_FittingDelMom(delOptProp_all, OptProp_base, n, rho, delMom, Mom_base, time_ns, cut_ind_mom, -1);

    Results.L(j,:) = X_found(1:2);
end

disp([char(datetime('now','Format','HH:mm:ss')) '  Finished for L'])


%% Musp1  and   Musp2
OptProp_base = opt_prop_base(1:6); % [Mua1 Mua2 Mua3 Musp1 Musp2 Musp3]

delOptProp_all = [-7 -7 -8 0 0 0 L_known 40]; % [delMua1 delMua2 delMua3 delMusp1 delMusp2 delMusp3 L1 L2]

for j = 1:length(ranges.musp1)
    OptProp_base(4) = ranges.musp1(j); % Vary the value of baseline Musp1
    
    % Use LMA #2 to retrieve delMua1 and delMua2
    [X_found, ~] = LMA_2_FittingDelMom(delOptProp_all, OptProp_base, n, rho, delMom, Mom_base, time_ns, cut_ind_mom, -1);

    Results.Musp1(j,:) = X_found(1:2);
end

% Repeating for Musp2:
OptProp_base = opt_prop_base(1:6); % [Mua1 Mua2 Mua3 Musp1 Musp2 Musp3]

for j = 1:length(ranges.musp2)
    OptProp_base(5:6) = ranges.musp2(j); % Vary the value of baseline Musp2
    
    % Use LMA #2 to retrieve delMua1 and delMua2
    [X_found, ~] = LMA_2_FittingDelMom(delOptProp_all, OptProp_base, n, rho, delMom, Mom_base, time_ns, cut_ind_mom, -1);

    Results.Musp2(j,:) = X_found(1:2);
end

disp([char(datetime('now','Format','HH:mm:ss')) '  Finished for Musp1 and Musp2'])


%% n1   and    n2
OptProp_base = opt_prop_base(1:6); % [Mua1 Mua2 Mua3 Musp1 Musp2 Musp3]

delOptProp_all = [-7 -7 -8 0 0 0 L_known 40]; % [delMua1 delMua2 delMua3 delMusp1 delMusp2 delMusp3 L1 L2]

n_temp = [n n n];
for j = 1:length(ranges.n1)
    n_temp(1) = ranges.n1(j); % Vary the value of n1
    
    % Use LMA #2 to retrieve delMua1 and delMua2
    [X_found, ~] = LMA_2_FittingDelMom(delOptProp_all, OptProp_base, n_temp, rho, delMom, Mom_base, time_ns, cut_ind_mom, -1);

    Results.n1(j,:) = X_found(1:2);
end

% Repeating for n2:
n_temp = [n n n];
for j = 1:length(ranges.n2)
    n_temp(2:3) = ranges.n2(j); % Vary the value of n2
    
    % Use LMA #2 to retrieve delMua1 and delMua2
    [X_found, ~] = LMA_2_FittingDelMom(delOptProp_all, OptProp_base, n_temp, rho, delMom, Mom_base, time_ns, cut_ind_mom, -1);

    Results.n2(j,:) = X_found(1:2);
end

disp([char(datetime('now','Format','HH:mm:ss')) '  Finished for n1 and n2'])


%% Mua1  and   Mua2
OptProp_base = opt_prop_base(1:6); % [Mua1 Mua2 Mua3 Musp1 Musp2 Musp3]

delOptProp_all = [-7 -7 -8 0 0 0 L_known 40]; % [delMua1 delMua2 delMua3 delMusp1 delMusp2 delMusp3 L1 L2]

for j = 1:length(ranges.mua1)
    OptProp_base(1) = ranges.mua1(j); % Vary the value of baseline Mua1
    
    % Use LMA #2 to retrieve delMua1 and delMua2
    [X_found, ~] = LMA_2_FittingDelMom(delOptProp_all, OptProp_base, n, rho, delMom, Mom_base, time_ns, cut_ind_mom, -1);

    Results.Mua1(j,:) = X_found(1:2);
end

% Repeating for Mua2:
OptProp_base = opt_prop_base(1:6); % [Mua1 Mua2 Mua3 Musp1 Musp2 Musp3]

for j = 1:length(ranges.mua1)
    OptProp_base(2:3) = ranges.mua2(j); % Vary the value of baseline Mua2
    
    % Use LMA #2 to retrieve delMua1 and delMua2
    [X_found, ~] = LMA_2_FittingDelMom(delOptProp_all, OptProp_base, n, rho, delMom, Mom_base, time_ns, cut_ind_mom, -1);

    Results.Mua2(j,:) = X_found(1:2);
end

disp([char(datetime('now','Format','HH:mm:ss')) '  Finished for Mua1 and Mua2'])

disp('DONE all')


%% Plot 

% set_dim = 1; % Setting this to 1 will make the figure sizes the same as we used for Publication 2023
set_dim = 0;

w_which_all = [1 2 3 4]; % Show all 4 panels in Fig. 3 of Publication 2023
% w_which_all = 1; % L
% w_which_all = 2; % Musp1 and Musp2
% w_which_all = 3; % n1 and n2
% w_which_all = 4; % Mua1 and Mua2

for w_which = w_which_all

    POS_3 = [36  50 770  600];
    hf = figure(5 + w_which); clf
    if set_dim == 1; hf.Position = POS_3; end

    if w_which == 4 % L
        X = ranges.L;
        X_true = L_known;
    
        Y1 = Results.L(:,1) * 10;
        Y2 = Results.L(:,2) * 10;
    
    elseif w_which == 2 % Musp1 and Musp2
        X = ranges.musp1 * 10;
        X_true = opt_prop_base(3) * 10;
        X_true_2 = opt_prop_base(4) * 10;
    
        Y1 = Results.Musp1(:,1) * 10;
        Y2 = Results.Musp1(:,2) * 10;
        Y3 = Results.Musp2(:,1) * 10;
        Y4 = Results.Musp2(:,2) * 10;
    
    elseif w_which == 3 % n1 and n2
        X = ranges.n1;
        X_true = n;
        X_true_2 = n;
    
        Y1 = Results.n1(:,1) * 10;
        Y2 = Results.n1(:,2) * 10;
        Y3 = Results.n2(:,1) * 10;
        Y4 = Results.n2(:,2) * 10;
    
    elseif w_which == 1 % Mua1 and Mua2
        X = ranges.mua1 * 10;
        X_true = opt_prop_base(1) * 10;
        X_true_2 = opt_prop_base(2) * 10;
        
        Y1 = Results.Mua1(:,1) * 10;
        Y2 = Results.Mua1(:,2) * 10;
        Y3 = Results.Mua2(:,1) * 10;
        Y4 = Results.Mua2(:,2) * 10;
    
    end
    
    lw = 2.5; lw2 = 2.5;
    ms = 15; ms2 = 12;
    
    if w_which == 4
        plot([X(1) X(end)],[1 1] * (opt_prop_after(1) - opt_prop_base(1)) * 10,'--','Color','green','LineWidth',2.5,'Handlevisibility','off')
        hold on
        plot([X(1) X(end)],[1 1] * (opt_prop_after(2) - opt_prop_base(2)) * 10,'--','Color','green','LineWidth',1.5,'Handlevisibility','off')

        plot([1 1] * X_true,[-0.2 0.2],'-','Color','green','LineWidth',1.5,'handlevisibility','off');
        
        H1 = plot(X, Y1, '-','Color','red','LineWidth',lw, 'MarkerSize',ms);
        H2 = plot(X, Y2, '-','Color','blue','LineWidth',lw, 'MarkerSize',ms);

        ylim([-0.0012 0.015] * 10)

%         xlim([X(1) X(end)])
%         xlim([10 20])
%         set(gca,'XTick',10:2:20)
        xlim([11 19])
        set(gca,'XTick',11:2:20)
        xlabel('Assumed  L / mm ')
        
        leg1 = legend([H1 H2],'\Delta\mu_{a,Sup.}','\Delta\mu_{a,Deep}','Location','NorthWest');
    
    else
        COL = [255 128 0] / 255;
        COL2 = [255 0 255] / 255;
        
        plot([X(1) X(end)],[1 1] * (opt_prop_after(1) - opt_prop_base(1)) * 10,'--','Color','green','LineWidth',1.5,'handlevisibility','off')
        hold on
        plot([X(1) X(end)],[1 1] * (opt_prop_after(2) - opt_prop_base(2)) * 10,'--','Color','green','LineWidth',1.5,'handlevisibility','off')

        plot([1 1] * X_true,[-0.2 0.2],'-','Color','green','LineWidth',1.5,'handlevisibility','off');

        plot([1 1] * X_true_2,[-0.2 0.2],'-','Color','green','LineWidth',1.5,'handlevisibility','off');
    
        H1 = plot(X, Y1, '-','Color','red','LineWidth',lw, 'MarkerSize',ms);
        H1.MarkerFaceColor = 'White';

        H2 = plot(X, Y2, '-','Color','blue','LineWidth',lw, 'MarkerSize',ms);
        H2.MarkerFaceColor = 'White';
    
        H3 = plot(X, Y3, '--','Color',COL2,'LineWidth',lw, 'MarkerSize',ms2);

        H4 = plot(X, Y4, '--','Color',COL,'LineWidth',lw, 'MarkerSize',ms2);
        
        ylim([-0.0012 0.015] * 10)

        set(gca,'YColor','black')
        
        xlim([X(1) X(end)])
    %     xlim([0.06 0.14])
        
        if w_which == 2
            xlabel('Assumed  \mu''_{s, Sup}  or  \mu''_{s, Deep} / cm^{-1}')
            xlim([8 12])
            leg1 = legend([H1 H2 H3 H4],'\Delta\mu_{a,Sup}   (for \mu''_{s,Sup})','\Delta\mu_{a,Deep} (for \mu''_{s,Sup})',...
                '\Delta\mu_{a,Sup}   (for \mu''_{s,Deep})','\Delta\mu_{a,Deep} (for \mu''_{s,Deep})','Location','NorthWest');

        elseif w_which == 3
            xlabel('Assumed  n_{Sup}  or  n_{Deep}')
            set(gca,'Xtick',[1.13 1.23 1.33 1.43 1.53])
            leg1 = legend([H1 H2 H3 H4],'\Delta\mu_{a,Sup}   (for n_{Sup})','\Delta\mu_{a,Deep} (for n_{Sup})',...
                '\Delta\mu_{a,Sup}   (for n_{Deep})','\Delta\mu_{a,Deep} (for n_{Deep})','Location','North');
            
        elseif w_which == 1
            xlabel('Assumed  \mu_{a, Sup}  or  \mu_{a, Deep} / cm^{-1}')
            xlim([0.06 0.14])
            leg1 = legend([H1 H2 H3 H4],'\Delta\mu_{a,Sup}   (for \mu_{a,Sup})','\Delta\mu_{a,Deep} (for \mu_{a,Sup})',...
                '\Delta\mu_{a,Sup}   (for \mu_{a,Deep})','\Delta\mu_{a,Deep} (for \mu_{a,Deep})','Location','North');
        end
    end

    grid on
    set(gca,'FontSize',18)

    ylabel('\Delta\mu_a / cm')
end