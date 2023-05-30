% Written by Aleh Sudakou (https://www.researchgate.net/profile/Aleh-Sudakou)
% Nałęcz Institute of Biocybernetics and Biomedical Engineering, Polish Academy of Sciences 
% This is one of the codes shared on  https://github.com/asudakou/Analyzing_TD-NIRS
% Last updated: 20 May 2023


%% Description:
% This script plots results of "Blood_main.m"
%
% Note, there are two ways to plot using this file: 
%     1. Call from "Blood_main.m"
%     2. Run the corresponding section in this script using "Run current section" (Ctrl + Enter)
%     This works because the names of variables inputted into this script are the same as those in "Blood_main.m"
%     And this script doesn't modify any variable, it only plots, so it can be called multiple times.
%
% print(gcf,'Name1.png','-dpng','-r600');  % A line of code for saving the current figure with high resolution in the curret directory


function Blood_plotting(Exp1_Time_min,Exp2_Time_min,Exp3_Time_min,    Exp1_DTOF,Exp2_DTOF,Exp3_DTOF,   Time_ns,   Exp1_IRF,Exp2_IRF,Exp3_IRF,   Base_index,...
            Exp1_Mom,Exp2_Mom,Exp3_Mom,  Exp1_LMA1_MuaMusp,Exp2_LMA1_MuaMusp,Exp3_LMA1_MuaMusp,...
            Exp1_LMA2_Mua1Mua2,Exp2_LMA2_Mua1Mua2,Exp3_LMA2_Mua1Mua2,...
            Exp1_Conc,Exp2_Conc,Exp3_Conc,  Exp1_StO2,Exp2_StO2,Exp3_StO2, set_dim, varargin)

% varargin  :  specifies which figure to plot. If it is empty, all figures will be plotted
% 
% set_dim   :  if it is 1, then the sizes of figure will be the same as we used for Publication 2023

% which_Fig must equal one of the following (or empty, for plotting all):
which_Fig = {'7a','7b','7c','8a','8b','9a','9b','9c','10','11',...
            '12a','12b','12c','13a','13b','13c','Optional_1', 'Optional_2', 'Optional_3'};

if ~isempty(varargin)
    temp = varargin{1};
    if sum(ismember(temp, which_Fig)) == 0
        error('Wrong input, check which figures can be plotted with this script: Fig. 7a to 13c')
    else
        which_Fig = temp;
    end
end

wavelengths = 680:12.5:868;

bg_TwoL = [0.5 1]; % Region for subtracting background noise
% bg_Deep = [0.5 1];
bg_irf = [7.3 8.3];


%% Fig 7 (a) - DTOFs
fig_num = 100;
if exist('which_Fig','var') &&  ~ismember('7a', which_Fig)

else

% for j_chan = 1:14 % To run in a loop and see results for each wavelength
for j_chan = 3 % For 2023 publication

    which_exp = 1;

    POS_1 = [70, 100, 800, 600];
    fig = figure(fig_num); clf
    if set_dim == 1; fig.Position = POS_1; end

    lw = 1.25;

    st = '-';

    plot_minutes = [0 12 47 44]; % StO2:  Both max, Deep min, Sup min, Both min.  Min StO2 corresponds to max Mua

    plot_index = zeros(length(plot_minutes),1);

    for j_1 = 1:size(plot_index,1)
        if     which_exp == 1; [~, plot_index(j_1)] = min(abs(plot_minutes(j_1) - Exp1_Time_min)); % Exp 1
        elseif which_exp == 2; [~, plot_index(j_1)] = min(abs(plot_minutes(j_1) - Exp2_Time_min)); % Exp 1
        elseif which_exp == 3; [~, plot_index(j_1)] = min(abs(plot_minutes(j_1) - Exp3_Time_min)); % Exp 3
        end
    end

    colcol = [0 00 204; 204 102 0; 0 180 0; 204 0 50]/255; % Colors
    
    for j = 1:4
        if     which_exp == 1; [ temp, ~ ] = DTOF_filter( squeeze(mean(Exp1_DTOF.TwoL(plot_index(j):plot_index(j)+10,j_chan,:),1)), Time_ns, bg_TwoL, 2, [0.25 0.03], 0);
        elseif which_exp == 2; [ temp, ~ ] = DTOF_filter( squeeze(mean(Exp2_DTOF.TwoL(plot_index(j):plot_index(j)+10,j_chan,:),1)), Time_ns, bg_TwoL, 2, [0.25 0.03], 0);
        elseif which_exp == 3; [ temp, ~ ] = DTOF_filter( squeeze(mean(Exp3_DTOF.TwoL(plot_index(j):plot_index(j)+10,j_chan,:),1)), Time_ns, bg_TwoL, 2, [0.25 0.03], 0);
        end
        semilogy(Time_ns, temp / max(temp),st,'LineWidth',lw,'Color',colcol(j,:),'Handlevisibility','off')
        hold on
        plot(-5,-5,st,'Color',colcol(j,:),'LineWidth',2)
    end

    xlim([1.5 10])
    set(gca,'Xtick',2:1:10)
    ylim([5e-5 1])
    set(gca,'Ytick',[1e-5 1e-4 1e-3 1e-2 1e-1 1])

    lw2 = 1.75;
    [ ~, temp] = DTOF_filter( squeeze(mean(Exp1_DTOF.TwoL(plot_index(1):plot_index(1)+10,j_chan,:),1)), Time_ns, bg_TwoL, 2, [0.25 0.03], 1);
    plot([1 1]*Time_ns(temp(1)),ylim,'--','Color',[0 0 204]/255,'LineWidth',lw2)
    plot([1 1]*Time_ns(temp(2)),ylim,'--','Color',[0 0 204]/255,'LineWidth',lw2,'Handlevisibility','off')

    grid on
    ylabel('Normalized photon count')
    xlabel('Time / ns')
    set(gca,'FontSize',24)

    % Plot also IRF (didn't plot in publication 2023)
    if     which_exp == 1; [ temp, ~ ] = DTOF_filter( squeeze(Exp1_IRF.TwoL_start(j_chan,:)), Time_ns, bg_irf, 2, [0 0], 0);
    elseif which_exp == 2; [ temp, ~ ] = DTOF_filter( squeeze(Exp2_IRF.TwoL_start(j_chan,:)), Time_ns, bg_irf, 2, [0 0], 0);
    elseif which_exp == 3; [ temp, ~ ] = DTOF_filter( squeeze(Exp3_IRF.TwoL_start(j_chan,:)), Time_ns, bg_irf, 2, [0 0], 0);
    end
    semilogy(Time_ns, temp / max(temp), '.','Color','green')

    legend('Both min','Deep max','Sup max','Both max','25% and 3%','IRF')
end
end


%% Fig 7 (b, c) - Moments
fig_num = 200;
if exist('which_Fig','var') &&  ~(ismember('7b', which_Fig) || ismember('7c', which_Fig))

else
    if ~exist('which_Fig','var')
        plot_here = [1 2]; % [1 2] to plot '7b' and '7c'
    else
        plot_here = [0 0];
        if ismember('7b', which_Fig); plot_here(1) = 1; end
        if ismember('7c', which_Fig); plot_here(2) = 2; end
    end
%     plot_here = 1;

% for j_chan = 1:14 % To run in a loop and see results for each wavelength
for j_chan = 3 % For 2023 publication

    for J = plot_here(plot_here~=0)
        if plot_here(J) > 0
            POS_1 = [70, 100, 1050, 600];
            fig = figure( fig_num + J ); clf
            if set_dim == 1; fig.Position = POS_1; end

            lw = 1.5;
%             lw = 3;
        
            if plot_here(J) == 1
                wm = 2; % Mean time of flight
            else
                wm = 3; % Variance
            end
            [~, temp] = min(abs(Exp1_Time_min));
            plot(Exp1_Time_min(1:end-11), Exp1_Mom.TwoL(1:end-11,j_chan,wm) - Exp1_Mom.TwoL(temp,j_chan,wm), '-', 'Color', 'red', 'LineWidth',lw)
            hold on
            plot(Exp1_Time_min(1:end-11), Exp1_Mom.Deep(1:end-11,j_chan,wm) - Exp1_Mom.Deep(temp,j_chan,wm), '-', 'Color', 'magenta', 'LineWidth',lw)

            [~, temp] = min(abs(Exp2_Time_min));
            plot(Exp2_Time_min, Exp2_Mom.TwoL(:,j_chan,wm) - Exp2_Mom.TwoL(temp,j_chan,wm), '-', 'Color', 'blue', 'LineWidth',lw)
            hold on
            plot(Exp2_Time_min, Exp2_Mom.Deep(:,j_chan,wm) - Exp2_Mom.Deep(temp,j_chan,wm), '-', 'Color', 'cyan', 'LineWidth',lw)

            [~, temp] = min(abs(Exp3_Time_min));
            plot(Exp3_Time_min, Exp3_Mom.TwoL(:,j_chan,wm) - Exp3_Mom.TwoL(temp,j_chan,wm), '-', 'Color', 'green', 'LineWidth',lw)
            hold on
            plot(Exp3_Time_min, Exp3_Mom.Deep(:,j_chan,wm) - Exp3_Mom.Deep(temp,j_chan,wm), '-', 'Color', 'black', 'LineWidth',lw)


            grid on
            set(gca,'FontSize',24)
            if J == 1
                ylim([-0.67 0.03])
                set(gca,'YTick',-0.6:0.2:0)
                ylabel('\Deltam_1 / ns')
            else
                ylim([-0.4 0.02])
                set(gca,'YTick',-0.4:0.1:0)
                ylabel('\DeltaV / ns^2')
            end
            xlim([-2 92])
            set(gca,'XTick',0:10:90)
            xlabel('Time / min')
        end
%         legend('L = 12 mm (Exp. #1)','L = 15 mm (Exp. #2)','L = 17 mm (Exp. #3)','Orientation','Horizontal')
    end
end
end


%% Fig 8 (a, b) - Spectra of Mua or Musp (curve-fitting)
fig_num = 300;
if exist('which_Fig','var') &&  ~(ismember('8a', which_Fig) || ismember('8b', which_Fig))

else
    if ~exist('which_Fig','var')
        plot_here = [1 2]; % [1 2] to plot '8a' and '8b'
    else
        plot_here = [0 0];
        if ismember('8a', which_Fig); plot_here(1) = 1; end
        if ismember('8b', which_Fig); plot_here(2) = 2; end
    end
%     plot_here = 2;


    % Water's Mua from literature  -  see "Ink_Nominal_mua.m"
    mua_Water_all = load('matcher94_nir_water_37.txt');
    mua_Water_all(:,2) = mua_Water_all(:,2) * log(10); % Absorption spectrum of Water from Matcher [OD cm-1]
    mua_Water = interp1(mua_Water_all(:,1),mua_Water_all(:,2),wavelengths); % At 16 wavelengths

%     % Another Mua spectra for water from literature:
%     mua_Water = GetExtinctions(wavelengths, 2); % [HbO Hb H2O lipid aa3]    [cm-1/(moles/liter)]   and   cm-1 for Water
%     mua_Water = mua_Water(:,3);

    Ext_coeff = [];
%     Ext_coeff(:,1) = 680:1:880;
    Ext_coeff(:,1) = wavelengths;
    wh_source = 1; % Gratzer
%     wh_source = 2; % Moaveni
%     wh_source = 3; % Takatani
    Ext_coeff(:,2:end+5) = GetExtinctions(Ext_coeff(:,1), wh_source); % [HbO Hb H2O lipid aa3]
    Ext_coeff(:,2:3) = Ext_coeff(:,2:3) .* log(10); %  [OD M-1 cm-1]    See "Ink_Nominal_mua.m"


    do_subt_water = 0; % Didn't subtract for the figure in 2023 Publication
%     do_subt_water = 1; % To subtract water mua

    if do_subt_water == 1
        mua_to_subtract = mua_Water(1:14) * (3366 + 180) / 3600; % Assuming SMOFlipid has the same Mua as Water
    else
        mua_to_subtract = 0;
    end

    for Mua_or_Musp = plot_here(plot_here~=0) % For plotting  Mua (set to 1), Musp (set to 2), or both (set to 1:2)
%     for Mua_or_Musp = 2
        
        fig = figure( fig_num + Mua_or_Musp ); clf
        if Mua_or_Musp == 1; POS_1 = [70, 96.5, 1060, 600]; else; POS_1 = [99   96.5  937  600]; end
        if set_dim == 1; fig.Position = POS_1; end
        
        lw = 2; lw_conn = 1.25;
        ms = 15;
        
        do_conn = 1; % Draw the connecting line?

        % Define colors:
        col_1 = {'red',[130 0 255]/255,'magenta';   'blue',[0 128 255]/255,'cyan'};
        st = {'+', 'x'};

        if Mua_or_Musp == 1
            yyaxis right
            plot(Ext_coeff(:,1), Ext_coeff(:,2) ./ 10^3, '.','Color','red','MarkerSize',25)
            hold on
            plot(Ext_coeff(:,1), Ext_coeff(:,2) ./ 10^3, '-','Color','red','LineWidth',1.5)
            plot(Ext_coeff(:,1), Ext_coeff(:,3) ./ 10^3, '.','Color','blue','MarkerSize',25)
            plot(Ext_coeff(:,1), Ext_coeff(:,3) ./ 10^3, '-','Color','blue','LineWidth',1.5)

            set(gca,'yColor','black')
            ylim([0 8])
            yyaxis left
            set(gca,'yColor','black')
        end
        
        j_chan_all = 1:14;
        for w_b = 1:2 % Loop for two baselines
            exp = 1; % Exp. 1
            plot(wavelengths(j_chan_all), squeeze(mean(Exp1_LMA1_MuaMusp.Deep(Base_index(w_b,1,exp):Base_index(w_b,2,exp),1:14, Mua_or_Musp),1)) * 10 - mua_to_subtract, st{w_b},'Color',col_1{w_b, exp},'LineWidth',lw,'MarkerSize',ms)
            hold on
            if do_conn
            plot(wavelengths(j_chan_all), squeeze(mean(Exp1_LMA1_MuaMusp.Deep(Base_index(w_b,1,exp):Base_index(w_b,2,exp),1:14, Mua_or_Musp), 1)) * 10 - mua_to_subtract, '--','Color',col_1{w_b, exp},'LineWidth',lw_conn);
            end
            
            
            exp = 2; % Exp. 2 (copy-pasted the above code)
            plot(wavelengths(j_chan_all), squeeze(mean(Exp2_LMA1_MuaMusp.Deep(Base_index(w_b,1,exp):Base_index(w_b,2,exp),1:14, Mua_or_Musp),1)) * 10 - mua_to_subtract, st{w_b},'Color',col_1{w_b, exp},'LineWidth',lw,'MarkerSize',ms)
            hold on
            if do_conn
            plot(wavelengths(j_chan_all), squeeze(mean(Exp2_LMA1_MuaMusp.Deep(Base_index(w_b,1,exp):Base_index(w_b,2,exp),1:14, Mua_or_Musp), 1)) * 10 - mua_to_subtract, '--','Color',col_1{w_b, exp},'LineWidth',lw_conn);
            end
            
            exp = 3; % Exp. 3 (copy-pasted the above code)
            plot(wavelengths(j_chan_all), squeeze(mean(Exp3_LMA1_MuaMusp.Deep(Base_index(w_b,1,exp):Base_index(w_b,2,exp),1:14, Mua_or_Musp),1)) * 10 - mua_to_subtract, st{w_b},'Color',col_1{w_b, exp},'LineWidth',lw,'MarkerSize',ms)
            hold on
            if do_conn
            plot(wavelengths(j_chan_all), squeeze(mean(Exp3_LMA1_MuaMusp.Deep(Base_index(w_b,1,exp):Base_index(w_b,2,exp),1:14, Mua_or_Musp), 1)) * 10 - mua_to_subtract, '--','Color',col_1{w_b, exp},'LineWidth',lw_conn);
            end
        end

        set(gca,'LooseInset',get(gca,'TightInset'));

%         legend('Exp. #1','Exp. #2','Exp. #3','Location','North')
%         legend('E of HbO2','E of Hb','Location','North')
        grid on
        box on
        set(gca,'FontSize',24)

        xlim([wavelengths(j_chan_all(1)) wavelengths(j_chan_all(end))])
        set(gca,'XTick',680:12.5*2:860)
        xlabel('Wavelength / nm')

        if Mua_or_Musp == 2
            ylabel('\mu''_s / cm^{-1}')
            ylim([9.0 12.0])
            set(gca,'YTick',8:1:12)
        else
            plot(mua_Water_all(:,1), mua_Water_all(:,2), '-','Color','green','LineWidth',1.5)

            ylim([0 0.24])
            set(gca,'YTick',0:0.06:0.3)
            ylabel('\mu_a / cm^{-1}')

            yyaxis right
            ylabel(['Molar Abs. Coef. (' char(949) ') / OD M^{-1} cm^{-1} 10^3'],'Color','black','FontSize',22)
        end
    end
end


%% Fig 10 - Musp (curve-fitting)
fig_num = 400;
if exist('which_Fig','var') &&  ~ismember('10', which_Fig)

else

%     for j_chan = 1:14
    for j_chan = [3 11] % Publication 2023
    
        fig = figure( fig_num + j_chan ); clf
        POS_1 = [70, 96.5, 1060, 600];
        if set_dim == 1; fig.Position = POS_1; end
    
        COL1 = {'red','blue','green'};
        COL2 = {'magenta','cyan','black'};
    
        lw = 1;
        
        plot(Exp1_Time_min, Exp1_LMA1_MuaMusp.TwoL(:,j_chan,2) * 10, 'Color',COL1{1},'LineWidth',lw)
        hold on
        
        plot(Exp2_Time_min, Exp2_LMA1_MuaMusp.TwoL(:,j_chan,2) * 10, 'Color',COL1{2},'LineWidth',lw)
    
        plot(Exp3_Time_min, Exp3_LMA1_MuaMusp.TwoL(:,j_chan,2) * 10, 'Color',COL1{3},'LineWidth',lw)
    
        plot(Exp1_Time_min, Exp1_LMA1_MuaMusp.Deep(:,j_chan,2) * 10, 'Color',COL2{1},'LineWidth',lw)
        
        plot(Exp2_Time_min, Exp2_LMA1_MuaMusp.Deep(:,j_chan,2) * 10, 'Color',COL2{2},'LineWidth',lw)
    
        plot(Exp3_Time_min, Exp3_LMA1_MuaMusp.Deep(:,j_chan,2) * 10, 'Color',COL2{3},'LineWidth',lw)
    
        grid on
    
        xlim([-2 92])
        set(gca,'XTick',0:10:90)
    
        set(gca,'FontSize',22)
    
        xlabel('Time / min')
        ylabel('\mu''_s / cm^{-1}')
    
    %     legend('Exp. #1 (L = 12 mm)','Exp. #2 (L = 15 mm)','Exp. #3 (L = 17 mm)','Exp. #1','Exp. #2','Exp. #3'); xlim([-50 -49]); grid off
    end
end


%% Fig 9 (a, b, c) - Determined Mua
fig_num = 500;
if exist('which_Fig','var') &&  ~(ismember('9a', which_Fig) || ismember('9b', which_Fig) || ismember('9c', which_Fig))

else
    if ~exist('which_Fig','var')
        plot_here = [1 2 3]; % [1 2 3] to plot '9a' and '9b' and '9c'
    else
        plot_here = [0 0 0];
        if ismember('9a', which_Fig); plot_here(1) = 1; end
        if ismember('9b', which_Fig); plot_here(2) = 2; end
        if ismember('9c', which_Fig); plot_here(3) = 3; end
    end

%     for j_chan = 3:11 % Use this and pause on each itteration to check for all wavelengths
    for j_chan = 3 % 3rd was plotted in Publication 2023
%     for j_chan = 1

        med_wind = 10;
%         med_wind = 1;

        lw = 1.5;
    
        for J = plot_here(plot_here~=0) % Exp1 Exp2 Exp3
%         for J = 1 % Only for Exp1

            fig = figure( fig_num + J ); clf
            POS_1 = [70, 100, 940, 600];
            if set_dim == 1; fig.Position = POS_1; end
            
            if J == 1
                temp = smooth(Exp1_LMA1_MuaMusp.Deep(:,j_chan,1) * 10, med_wind); temp(end-14:end) = str2double('NaN');
                plot(Exp1_Time_min, temp,'Color','black','LineWidth',lw)
                hold on
    
                temp = smooth(Exp1_LMA1_MuaMusp.TwoL(:,j_chan,1) * 10, med_wind); temp(end-14:end) = str2double('NaN');
                plot(Exp1_Time_min, temp,'Color','green','LineWidth',lw)
    
                temp = smooth(Exp1_LMA2_Mua1Mua2.TwoL(:,j_chan,1) * 10, med_wind); temp(end-14:end) = str2double('NaN');
                plot(Exp1_Time_min, temp,'Color','magenta','LineWidth',lw)
    
                temp = smooth(Exp1_LMA2_Mua1Mua2.TwoL(:,j_chan,2) * 10, med_wind); temp(end-14:end) = str2double('NaN');
                plot(Exp1_Time_min, temp,'Color','cyan','LineWidth',lw)

            elseif J == 2
                temp = smooth(Exp2_LMA1_MuaMusp.Deep(:,j_chan,1) * 10, med_wind); temp(end-med_wind:end) = str2double('NaN');
                plot(Exp2_Time_min, temp,'Color','black','LineWidth',lw)
                hold on
    
                temp = smooth(Exp2_LMA1_MuaMusp.TwoL(:,j_chan,1) * 10, med_wind); temp(end-med_wind:end) = str2double('NaN');
                plot(Exp2_Time_min, temp,'Color','green','LineWidth',lw)
    
                temp = smooth(Exp2_LMA2_Mua1Mua2.TwoL(:,j_chan,1) * 10, med_wind); temp(end-med_wind:end) = str2double('NaN');
                plot(Exp2_Time_min, temp,'Color','magenta','LineWidth',lw)
    
                temp = smooth(Exp2_LMA2_Mua1Mua2.TwoL(:,j_chan,2) * 10, med_wind); temp(end-med_wind:end) = str2double('NaN');
                plot(Exp2_Time_min, temp,'Color','cyan','LineWidth',lw)

            elseif J == 3
                temp = smooth(Exp3_LMA1_MuaMusp.Deep(:,j_chan,1) * 10, med_wind); temp(end-med_wind:end) = str2double('NaN');
                plot(Exp3_Time_min, temp,'Color','black','LineWidth',lw)
                hold on

                temp = smooth(Exp3_LMA1_MuaMusp.TwoL(:,j_chan,1) * 10, med_wind); temp(end-med_wind:end) = str2double('NaN');
                plot(Exp3_Time_min, temp,'Color','green','LineWidth',lw)
    
                temp = smooth(Exp3_LMA2_Mua1Mua2.TwoL(:,j_chan,1) * 10, med_wind); temp(end-med_wind:end) = str2double('NaN');
                plot(Exp3_Time_min, temp,'Color','magenta','LineWidth',lw)
    
                temp = smooth(Exp3_LMA2_Mua1Mua2.TwoL(:,j_chan,2) * 10, med_wind); temp(end-med_wind:end) = str2double('NaN');
                plot(Exp3_Time_min, temp,'Color','cyan','LineWidth',lw)
            end
    
            grid on

            ylim([0.035 0.190001])
            set(gca,'YTick',0.04:0.03:0.2)

            xlim([-2 92])
            set(gca,'XTick',0:10:90)

            set(gca,'FontSize',24)

            POSPOS = [0.15   0.1981    0.7758    0.7596];
            if set_dim == 1; set(gca,'Position',POSPOS); end

            xlabel('Time / min')
            ylabel('\mu_a / cm^{-1}')

%             print(gcf,['Name' num2str(J) '.png'],'-dpng','-r600')
        end
    end
end


%% Fig 11 (a to f) - Determined Conc
fig_num = 600;
if exist('which_Fig','var') &&  ~ismember('11', which_Fig)

else
    
    lw = 1.5;
    med_wind = 10;

    for J = 1:6 % for 6 panels in figure 11
%     for J = 6

%         set_dim = 0;
%         set_dim = 1;

        POS_1 = [90, 96.5, 950, 600];

        if J == 1 % panel (a)
            fig = figure( fig_num + J ); clf
            if set_dim == 1; fig.Position = POS_1; end

            temp = smooth(Exp1_Conc.Deep_LMA1(:,1), med_wind); temp(end-14:end) = str2double('NaN');
            plot(Exp1_Time_min, temp,'-','Color','red','LineWidth',lw) % Oxy
            hold on
        
            temp = smooth(Exp1_Conc.Deep_LMA1(:,2), med_wind); temp(end-14:end) = str2double('NaN');
            plot(Exp1_Time_min, temp,'-','Color','blue','LineWidth',lw) % Deoxy
        
            temp = smooth(Exp1_Conc.TwoL_LMA2(:,1,2), med_wind); temp(end-14:end) = str2double('NaN');
            plot(Exp1_Time_min, temp,'-','Color','magenta','LineWidth',lw) % Oxy
        
            temp = smooth(Exp1_Conc.TwoL_LMA2(:,2,2), med_wind); temp(end-14:end) = str2double('NaN');
            plot(Exp1_Time_min, temp,'-','Color','cyan','LineWidth',lw) % Deoxy

        elseif J == 2 % panel (b)
            fig = figure( fig_num + J ); clf
            if set_dim == 1; fig.Position = POS_1; end

            temp = smooth(Exp1_Conc.TwoL_LMA1(:,1), med_wind); temp(end-14:end) = str2double('NaN');
            plot(Exp1_Time_min, temp,'-','Color','red','LineWidth',lw) % Oxy
            hold on
        
            temp = smooth(Exp1_Conc.TwoL_LMA1(:,2), med_wind); temp(end-14:end) = str2double('NaN');
            plot(Exp1_Time_min, temp,'-','Color','blue','LineWidth',lw) % Deoxy
        
            temp = smooth(Exp1_Conc.TwoL_LMA2(:,1,1), med_wind); temp(end-14:end) = str2double('NaN');
            plot(Exp1_Time_min, temp,'-','Color','magenta','LineWidth',lw) % Oxy
        
            temp = smooth(Exp1_Conc.TwoL_LMA2(:,2,1), med_wind); temp(end-14:end) = str2double('NaN');
            plot(Exp1_Time_min, temp,'-','Color','cyan','LineWidth',lw) % Deoxy


            % Repeat the above for all other experiments:
        elseif J == 3 % panel (c)
            fig = figure( fig_num + J ); clf
            if set_dim == 1; fig.Position = POS_1; end

            temp = smooth(Exp2_Conc.Deep_LMA1(:,1), med_wind); temp(end-med_wind:end) = str2double('NaN');
            plot(Exp2_Time_min, temp,'-','Color','red','LineWidth',lw) % Oxy
            hold on
        
            temp = smooth(Exp2_Conc.Deep_LMA1(:,2), med_wind); temp(end-med_wind:end) = str2double('NaN');
            plot(Exp2_Time_min, temp,'-','Color','blue','LineWidth',lw) % Deoxy
        
            temp = smooth(Exp2_Conc.TwoL_LMA2(:,1,2), med_wind); temp(end-med_wind:end) = str2double('NaN');
            plot(Exp2_Time_min, temp,'-','Color','magenta','LineWidth',lw) % Oxy
        
            temp = smooth(Exp2_Conc.TwoL_LMA2(:,2,2), med_wind); temp(end-med_wind:end) = str2double('NaN');
            plot(Exp2_Time_min, temp,'-','Color','cyan','LineWidth',lw) % Deoxy

        elseif J == 4 % panel (d)
            fig = figure( fig_num + J ); clf
            if set_dim == 1; fig.Position = POS_1; end

            temp = smooth(Exp2_Conc.TwoL_LMA1(:,1), med_wind); temp(end-med_wind:end) = str2double('NaN');
            plot(Exp2_Time_min, temp,'-','Color','red','LineWidth',lw) % Oxy
            hold on
        
            temp = smooth(Exp2_Conc.TwoL_LMA1(:,2), med_wind); temp(end-med_wind:end) = str2double('NaN');
            plot(Exp2_Time_min, temp,'-','Color','blue','LineWidth',lw) % Deoxy
        
            temp = smooth(Exp2_Conc.TwoL_LMA2(:,1,1), med_wind); temp(end-med_wind:end) = str2double('NaN');
            plot(Exp2_Time_min, temp,'-','Color','magenta','LineWidth',lw) % Oxy
        
            temp = smooth(Exp2_Conc.TwoL_LMA2(:,2,1), med_wind); temp(end-med_wind:end) = str2double('NaN');
            plot(Exp2_Time_min, temp,'-','Color','cyan','LineWidth',lw) % Deoxy

        elseif J == 5 % panel (e)
            fig = figure( fig_num + J ); clf
            if set_dim == 1; fig.Position = POS_1; end

            temp = smooth(Exp3_Conc.Deep_LMA1(:,1), med_wind); temp(end-med_wind:end) = str2double('NaN');
            plot(Exp3_Time_min, temp,'-','Color','red','LineWidth',lw) % Oxy
            hold on
        
            temp = smooth(Exp3_Conc.Deep_LMA1(:,2), med_wind); temp(end-med_wind:end) = str2double('NaN');
            plot(Exp3_Time_min, temp,'-','Color','blue','LineWidth',lw) % Deoxy
        
            temp = smooth(Exp3_Conc.TwoL_LMA2(:,1,2), med_wind); temp(end-med_wind:end) = str2double('NaN');
            plot(Exp3_Time_min, temp,'-','Color','magenta','LineWidth',lw) % Oxy
        
            temp = smooth(Exp3_Conc.TwoL_LMA2(:,2,2), med_wind); temp(end-med_wind:end) = str2double('NaN');
            plot(Exp3_Time_min, temp,'-','Color','cyan','LineWidth',lw) % Deoxy

        elseif J == 6 % panel (f)
            fig = figure( fig_num + J ); clf
            if set_dim == 1; fig.Position = POS_1; end

            temp = smooth(Exp3_Conc.TwoL_LMA1(:,1), med_wind); temp(end-med_wind:end) = str2double('NaN');
            plot(Exp3_Time_min, temp,'-','Color','red','LineWidth',lw) % Oxy
            hold on
        
            temp = smooth(Exp3_Conc.TwoL_LMA1(:,2), med_wind); temp(end-med_wind:end) = str2double('NaN');
            plot(Exp3_Time_min, temp,'-','Color','blue','LineWidth',lw) % Deoxy
        
            temp = smooth(Exp3_Conc.TwoL_LMA2(:,1,1), med_wind); temp(end-med_wind:end) = str2double('NaN');
            plot(Exp3_Time_min, temp,'-','Color','magenta','LineWidth',lw) % Oxy
        
            temp = smooth(Exp3_Conc.TwoL_LMA2(:,2,1), med_wind); temp(end-med_wind:end) = str2double('NaN');
            plot(Exp3_Time_min, temp,'-','Color','cyan','LineWidth',lw) % Deoxy

        end

        xlim([-2 92])
        set(gca,'XTick',0:10:90)

        ylim([-5 45])

        grid on
        box on

        set(gca,'FontSize',24)

        POSPOS = [0.15   0.1981    0.7758    0.7596];
        if set_dim == 1; set(gca,'Position',POSPOS); end

        xlabel('Time / min')
        ylabel(['Concentration / ' char(181) 'M'])

%         print(gcf,['Name' num2str(J) '.png'],'-dpng','-r600')
    end
end


%% Fig 12 (a, b, c) - Determined StO2
fig_num = 700;
if exist('which_Fig','var') &&  ~(ismember('12a', which_Fig) || ismember('12b', which_Fig) || ismember('12c', which_Fig))

else
    if ~exist('which_Fig','var')
        plot_here = [1 2 3]; % [1 2 3] to plot '9a' and '9b' and '9c'
    else
        plot_here = [0 0 0];
        if ismember('12a', which_Fig); plot_here(1) = 1; end
        if ismember('12b', which_Fig); plot_here(2) = 2; end
        if ismember('12c', which_Fig); plot_here(3) = 3; end
    end
    
    med_wind = 10;
%     med_wind = 1;

    lw = 1.5;
%     st = 'x'; st2 = '+';
    st = '-'; st2 = '-';

    for J = plot_here(plot_here~=0) % Exp1 Exp2 Exp3
%     for J = 2

        fig = figure( fig_num + J ); clf
        POS_1 = [70, 100, 940, 600];
        if set_dim == 1; fig.Position = POS_1; end
        
        if J == 1
            j_cycle_all = 1:size(Exp1_StO2.Deep_LMA1,1);

            temp = smooth(Exp1_StO2.Deep_LMA1, med_wind); temp(end-14:end) = str2double('NaN');
            plot(Exp1_Time_min(j_cycle_all), temp(j_cycle_all),st2,'Color','black','LineWidth',lw)
            hold on

            temp = smooth(Exp1_StO2.TwoL_LMA1, med_wind); temp(end-14:end) = str2double('NaN');
            plot(Exp1_Time_min(j_cycle_all), temp(j_cycle_all),st2,'Color','green','LineWidth',lw)

            temp = smooth(Exp1_StO2.TwoL_LMA2(:,1), med_wind); temp(end-14:end) = str2double('NaN');
            plot(Exp1_Time_min(j_cycle_all), temp(j_cycle_all),st,'Color','magenta','LineWidth',lw)

            temp = smooth(Exp1_StO2.TwoL_LMA2(:,2), med_wind); temp(end-14:end) = str2double('NaN');
            plot(Exp1_Time_min(j_cycle_all), temp(j_cycle_all),st,'Color','cyan','LineWidth',lw)

        elseif J == 2
            j_cycle_all = 1:size(Exp2_StO2.Deep_LMA1,1);

            temp = smooth(Exp2_StO2.Deep_LMA1, med_wind); temp(end-med_wind:end) = str2double('NaN');
            plot(Exp2_Time_min(j_cycle_all), temp(j_cycle_all),st2,'Color','black','LineWidth',lw)
            hold on

            temp = smooth(Exp2_StO2.TwoL_LMA1, med_wind); temp(end-med_wind:end) = str2double('NaN');
            plot(Exp2_Time_min(j_cycle_all), temp(j_cycle_all),st2,'Color','green','LineWidth',lw)

            temp = smooth(Exp2_StO2.TwoL_LMA2(:,1), med_wind); temp(end-med_wind:end) = str2double('NaN');
            plot(Exp2_Time_min(j_cycle_all), temp(j_cycle_all),st,'Color','magenta','LineWidth',lw)

            temp = smooth(Exp2_StO2.TwoL_LMA2(:,2), med_wind); temp(end-med_wind:end) = str2double('NaN');
            plot(Exp2_Time_min(j_cycle_all), temp(j_cycle_all),st,'Color','cyan','LineWidth',lw)

        elseif J == 3
            j_cycle_all = 1:size(Exp3_StO2.Deep_LMA1,1);

            temp = smooth(Exp3_StO2.Deep_LMA1, med_wind); temp(end-med_wind:end) = str2double('NaN');
            plot(Exp3_Time_min(j_cycle_all), temp(j_cycle_all),st2,'Color','black','LineWidth',lw)
            hold on

            temp = smooth(Exp3_StO2.TwoL_LMA1, med_wind); temp(end-med_wind:end) = str2double('NaN');
            plot(Exp3_Time_min(j_cycle_all), temp(j_cycle_all),st2,'Color','green','LineWidth',lw)

            temp = smooth(Exp3_StO2.TwoL_LMA2(:,1), med_wind); temp(end-med_wind:end) = str2double('NaN');
            plot(Exp3_Time_min(j_cycle_all), temp(j_cycle_all),st,'Color','magenta','LineWidth',lw)

            temp = smooth(Exp3_StO2.TwoL_LMA2(:,2), med_wind); temp(end-med_wind:end) = str2double('NaN');
            plot(Exp3_Time_min(j_cycle_all), temp(j_cycle_all),st,'Color','cyan','LineWidth',lw)
        end

        grid on

        xlim([-2 92])
        set(gca,'XTick',0:10:90)

        ylim([-10 105])

        grid on
        box on

        set(gca,'FontSize',24)

        POSPOS = [0.15   0.1981    0.7758    0.7596];
        if set_dim == 1; set(gca,'Position',POSPOS); end

        xlabel('Time / min')
        if J == 1
            ylabel('StO_2 / %')
        end

%         print(gcf,['Name' num2str(J) '.png'],'-dpng','-r600')
    end
end



%% Fig 13 (a, b, c) - Determined StO2  vs.  StO2
fig_num = 800;
if exist('which_Fig','var') &&  ~(ismember('13a', which_Fig) || ismember('13b', which_Fig) || ismember('13c', which_Fig))

else
    if ~exist('which_Fig','var')
        plot_here = [1 2 3]; % [1 2 3] to plot '9a' and '9b' and '9c'
    else
        plot_here = [0 0 0];
        if ismember('13a', which_Fig); plot_here(1) = 1; end
        if ismember('13b', which_Fig); plot_here(2) = 2; end
        if ismember('13c', which_Fig); plot_here(3) = 3; end
    end
    
%     plot_minutes = [0 10; 15 29.5; 47 59.5; 62 76.5; 79 90];
    plot_minutes = [0 12.8;   15.2 29.8;   47.2 59.8;   62.2 76.8;   79.2 90]; % Publication 2023

    mean_wind = 10;

    colcol = {'red','blue','green','black','magenta'};

    lw = 1.5;

    for J = plot_here(plot_here~=0) % Exp1 Exp2 Exp3
%     for J = 1
        
        fig = figure( fig_num + J ); clf
        POS_1 = [70, 96.5, 1060, 600];
        if set_dim == 1; fig.Position = POS_1; end

        plot([-30 130], [-30 130],':','Color',[1 1 1]*150/255,'LineWidth',2,'handlevisibility','off')
        hold on

        if J == 1
            Y_1 = smooth(Exp1_StO2.Deep_LMA1, mean_wind);
            Y_2 = smooth(Exp1_StO2.TwoL_LMA2(:,2), mean_wind);

            for j_1 = 1:size(plot_minutes,1)
                [~, plot_ind_1] = min(abs(plot_minutes(j_1,1) - Exp1_Time_min)); % Exp 1
                [~, plot_ind_2] = min(abs(plot_minutes(j_1,2) - Exp1_Time_min));

                plot( Y_1(plot_ind_1:plot_ind_2), Y_2(plot_ind_1:plot_ind_2), '-', 'Color', colcol{j_1},'MarkerSize',7, 'LineWidth',lw)
                hold on
            end

        elseif J == 2
            Y_1 = smooth(Exp2_StO2.Deep_LMA1, mean_wind);
            Y_2 = smooth(Exp2_StO2.TwoL_LMA2(:,2), mean_wind);

            for j_1 = 1:size(plot_minutes,1)
                [~, plot_ind_1] = min(abs(plot_minutes(j_1,1) - Exp2_Time_min)); % Exp 2
                [~, plot_ind_2] = min(abs(plot_minutes(j_1,2) - Exp2_Time_min));

                plot( Y_1(plot_ind_1:plot_ind_2), Y_2(plot_ind_1:plot_ind_2), '-', 'Color', colcol{j_1},'MarkerSize',7, 'LineWidth',lw)
                hold on
            end

        elseif J == 3
            Y_1 = smooth(Exp3_StO2.Deep_LMA1, mean_wind);
            Y_2 = smooth(Exp3_StO2.TwoL_LMA2(:,2), mean_wind);

            for j_1 = 1:size(plot_minutes,1)
                [~, plot_ind_1] = min(abs(plot_minutes(j_1,1) - Exp3_Time_min)); % Exp 3
                [~, plot_ind_2] = min(abs(plot_minutes(j_1,2) - Exp3_Time_min));

                plot( Y_1(plot_ind_1:plot_ind_2), Y_2(plot_ind_1:plot_ind_2), '-', 'Color', colcol{j_1},'MarkerSize',7, 'LineWidth',lw)
                hold on
            end
        end

%         legend('Cycle 1','Cycle 2','Cycle 4','Cycle 5','Cycle 6','Identity line','Location','NorthWest')

        ylabel('StO_2 for Pos. #2 (2-L. Deep) / %')
        xlabel('StO_2 for Pos. #1 (Homog.) / %')

        grid on

        set(gca,'FontSize',22)

        axis equal

        ylim([-5 105])
        xlim([-5 105])
        set(gca,'YTick',-20:20:105)
        set(gca,'XTick',-20:20:105)

%         print(gcf,['Name' num2str(J) '.png'],'-dpng','-r600')
    end
end
end