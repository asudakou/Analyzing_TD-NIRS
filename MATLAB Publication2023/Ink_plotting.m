% Written by Aleh Sudakou (https://www.researchgate.net/profile/Aleh-Sudakou)
% Nałęcz Institute of Biocybernetics and Biomedical Engineering, Polish Academy of Sciences 
% This is one of the codes shared on  https://github.com/asudakou/Analyzing_TD-NIRS
% Last updated: 20 May 2023


%% Description:
% This script plots results of "Ink_main.m"
%
% Note, there are two ways to plot using this file: 
%     1. Call from "Ink_main.m"
%     2. Run the corresponding section in this script using "Run current section" (Ctrl + Enter)
%     This works because the names of variables inputted into this script are the same as those in "Ink_main.m"
%     And this script doesn't modify any variable, it only plots, so it can be called multiple times.
%
% print(gcf,'Name1.png','-dpng','-r600');  % A line of code for saving the current figure with high resolution in the curret directory


function Ink_plotting( DTOF, Time_ns, IRF, delMom, LMA1_MuaMusp, LMA2_Mua1Mua2, LMA1_Optional, set_dim, which_chan, varargin )


% varargin  :  specifies which figure to plot. If it is empty, all figures will be plotted
% 
% set_dim   :  if it is 1, then the sizes of figure will be the same as we used for Publication 2023

% which_Fig must equal one of the following (or empty, for plotting all):
which_Fig = {'4a','4b','5a','5b','6a','6b', 'Optional_1', 'Optional_2', 'Optional_3'};

if ~isempty(varargin)
    temp = varargin{1};
    if sum(ismember(temp, which_Fig)) == 0
        error('Wrong input, check which figures can be plotted with this script: Fig. 4a to 6b or Optional 1 to 3')
    else
        which_Fig = temp;
    end
end

wavelengths = 680:12.5:868;

bg_Red = [0.49 0.92]; % Region for subtracting background noise
bg_Blu = [0.49 0.92]; % For the second detector
bg_irf = [7.3 8.3]; % For IRF

mua_nom = Ink_Nominal_mua(wavelengths); % Fo x-axis in some figures


%% Fig 5 (a) - DTOFs and IRFs
fig_num = 50;
if exist('which_Fig','var') &&  ~ismember('5a', which_Fig)

else
    fig = figure(fig_num); clf
    POS_1 = [69, 96.5, 850, 630];
    if set_dim == 1; fig.Position = POS_1; end
    
%     j_chan = 8; % 8 for 2023 publication
    j_chan = which_chan; 
    
    j_MuaStep = 6; % 6, which is the 5th Mua step, for 2023 publication
    
    lw_dtof = 1.25;
    lw_ind = 2.2;
    st_dtof = '-';

    offset = 0.06/physconst('Lightspeed')*1e9; % ns. Time offset of IRF due the 6 cm distance between fibers when measuring IRF
    
    temp_col = {'red','blue','green','black','blue'};
    
    [ dtof_1, ~ ] = DTOF_filter( squeeze(DTOF.Exp1_Deep_Red(j_MuaStep,j_chan,:)), Time_ns, bg_Red, 2, 0, 0);
    [ dtof_2, ~ ] = DTOF_filter( squeeze(DTOF.Exp1_TwoL_Blu(j_MuaStep,j_chan,:)), Time_ns, bg_Blu, 2, 0, 0);
    [ ~, ind_mom] = DTOF_filter( squeeze(DTOF.Exp1_TwoL_Blu(11,j_chan,:)), Time_ns, bg_Blu, 2, [0.25 0.03], 1);
    [ irf_1, ~ ] =  DTOF_filter( squeeze(IRF.Exp1_Red(:,j_chan)), Time_ns, bg_irf, 2, 0, 0);
    [ irf_2, ~ ] =  DTOF_filter( squeeze(IRF.Exp1_Blu(:,j_chan)), Time_ns, bg_irf, 2, 0, 0);

    % Note, only for this figure, we shift DTOFs and IRFs for 1 detector by 0.61 ns so that it aligns with the other detector. 
    temp = 0.61; 

    semilogy(Time_ns + temp, dtof_1 / max(dtof_1),st_dtof,'LineWidth',lw_dtof,'Color',temp_col{1},'Handlevisibility','off')
    hold on

    semilogy(Time_ns + temp - offset, irf_1 / max(irf_1),st_dtof,'LineWidth',lw_dtof,'Color',temp_col{3},'Handlevisibility','off')

    semilogy(Time_ns, dtof_2 / max(dtof_2),st_dtof,'LineWidth',lw_dtof,'Color',temp_col{2},'Handlevisibility','off')

    semilogy(Time_ns - offset, irf_2 / max(irf_2),st_dtof,'LineWidth',lw_dtof,'Color',temp_col{4},'Handlevisibility','off')

    plot([1 1]*Time_ns(ind_mom(1)), ylim, '--','Color',temp_col{5},'LineWidth',lw_ind,'Handlevisibility','off')
    plot([1 1]*Time_ns(ind_mom(2)), ylim, '--','Color',temp_col{5},'LineWidth',lw_ind,'Handlevisibility','off')

    % For the legend:
    plot(-5,0,st_dtof,'Color',temp_col{1},'LineWidth',2,'MarkerSize',30)
    plot(-5,0,st_dtof,'Color',temp_col{2},'LineWidth',2,'MarkerSize',30)
    plot(-5,0,st_dtof,'Color',temp_col{3},'LineWidth',2,'MarkerSize',30)
    plot(-5,0,st_dtof,'Color',temp_col{4},'LineWidth',2,'MarkerSize',30)
    plot(-5,0,'--','Color',temp_col{5},'LineWidth',1.85)

    legend('DTOF (on Pos. #1)','DTOF (on Pos. #2)','IRF (Detector #1)','IRF (Detector #2)','25% and 3%','FontSize',19)

    xlim([1 8])
    ylim([1e-5 1])
    set(gca,'Ytick',[1e-5 1e-4 1e-3 1e-2 1e-1 1])

    grid on
    ylabel('Normalized photon count')
    xlabel('Time / ns')

    POSPOS = [0.1632    0.15    0.7758    0.7596];
    if set_dim == 1; set(gca,'Position',POSPOS); end
    set(gca,'FontSize',22)

    title(['Spectral channel : ' num2str(round(wavelengths(j_chan))) ' nm'])
end


%% Fig 4 (a) - Spectra of Mua or Musp (curve-fitting)
fig_num = 150;
if exist('which_Fig','var') &&  ~ismember('4a', which_Fig)

else
    for Mua_or_Musp = 1 % For plotting  Mua (set to 1), Musp (set to 2), or both (set to 1:2)
%     for Mua_or_Musp = 2
    
        fig = figure( fig_num + Mua_or_Musp ); clf
        POS_1 = [69 + (Mua_or_Musp-1)*200, 96.5, 1050, 600];
        if set_dim == 1; fig.Position = POS_1; end
        
        j_chan_all = 1:14; % Which channels to show
        
        j_steps_all = [1 11 21]; % For 2023 Publication
%         j_steps_all = [1 11 16 21]; % Which Mua steps to show
        
        lw_dtof = 2; lw_conn = 1;
        ms = 15; ms2 = 15;
        
        do_conn = 1; % Draw the connecting line?

        % Water's Mua from literature  -  see "Ink_Nominal_mua.m"
        mua_Water_all = load('matcher94_nir_water_37.txt');
        mua_Water_all(:,2) = mua_Water_all(:,2) * log(10); % [cm-1]   Absorption spectrum of Water from Matcher in OD/cm
        mua_Water = interp1(mua_Water_all(:,1),mua_Water_all(:,2),wavelengths); % At 16 wavelengths

        % Another Mua spectra for water from literature:
%         mua_Water = GetExtinctions(wavelengths, 2); % [HbO Hb H2O lipid aa3]    [cm-1/(moles/liter)]   and   cm-1 for Water
%         mua_Water = mua_Water(:,3);

        do_subt_water = 0; % Didn't subtract for the figure in 2023 Publication 
%         do_subt_water = 1; % To subtract water mua

        if do_subt_water == 1
            mua_to_subtract = mua_Water(j_chan_all)';
        else
            mua_to_subtract = 0;
        end
        
        % List of colors:
        temp_col = {'black','red',[130 0 255]/255,'magenta',[255 128 0]/255,'blue',[0 128 255]/255,'cyan','green','cyan','magenta'};
        if length(temp_col) < 2 + length(j_steps_all) * 2
            for j = length(temp_col)+1 : 2 + length(j_steps_all) * 2
                temp_col{j} = 'red';
            end
        end
    
        % 1. Baseline Superficial
        temp_mua = squeeze(LMA1_MuaMusp.ConstLayer_Exp1_SupL_Blu(1,j_chan_all,Mua_or_Musp)) * 10 - mua_to_subtract;
        plot(wavelengths(j_chan_all), temp_mua,'+','Color',temp_col{1},'LineWidth',lw_dtof,'MarkerSize',ms)
        hold on
        if do_conn; plot(wavelengths(j_chan_all), temp_mua,'--','Color',temp_col{1},'LineWidth',lw_conn,'Handlevisibility','off'); end
    
        % 2. Deep for Mua steps
        for j_step = 1:length(j_steps_all)
            temp_mua = squeeze(LMA1_MuaMusp.Exp1_Deep_Red(j_steps_all(j_step),j_chan_all,Mua_or_Musp)) * 10 - mua_to_subtract;
    
            plot(wavelengths(j_chan_all), temp_mua,'+','Color',temp_col{1 + j_step},'LineWidth',lw_dtof,'MarkerSize',ms2)
            if do_conn; plot(wavelengths(j_chan_all), temp_mua,'--','Color',temp_col{1 + j_step},'LineWidth',lw_conn,'Handlevisibility','off'); end
        end
        j_step = j_step + 2;
    
        % 3. Baseline Deep
        temp_mua = squeeze(LMA1_MuaMusp.ConstLayer_Exp2_Deep_Blu(1,j_chan_all,Mua_or_Musp)) * 10 - mua_to_subtract;
        plot(wavelengths(j_chan_all), temp_mua,'x','Color',temp_col{j_step},'LineWidth',lw_dtof,'MarkerSize',ms)
        if do_conn; plot(wavelengths(j_chan_all), temp_mua,'--','Color',temp_col{j_step},'LineWidth',lw_conn,'Handlevisibility','off'); end
    
        % 4. Superficial for Mua steps
        for jj2 = 1:length(j_steps_all)
            plot(mua_Water_all(:,1),mua_Water_all(:,2),'-','Color','green','LineWidth',1.5)

            temp_mua = squeeze(LMA1_MuaMusp.Exp2_SupL_Red(j_steps_all(jj2),j_chan_all,Mua_or_Musp)) * 10 - mua_to_subtract;
    
            plot(wavelengths(j_chan_all), temp_mua,'x','Color',temp_col{j_step+jj2},'LineWidth',lw_dtof,'MarkerSize',ms2)
            if do_conn; plot(wavelengths(j_chan_all), temp_mua,'--','Color',temp_col{j_step+jj2},'LineWidth',lw_conn,'Handlevisibility','off'); end
        end
    
        if Mua_or_Musp == 1
            ylabel('\mu_a / cm^{-1}')
            ylim([0 0.25001])
            set(gca,'YTick',0.00:0.05:0.25)
        else
            ylabel('\mu''_s / cm^{-1}')
            ylim([9 14])
        end
        
        grid on
        xlim([wavelengths(min(j_chan_all)) wavelengths(max(j_chan_all))])
        set(gca,'XTick',round(wavelengths(1:2:end)))
        xlabel('Wavelength / nm')
        set(gca,'FontSize',24)
        box on
    end
%     legend('Sup. (Constant)','Deep (0^{th} step)','Deep (10^{th} step)','Deep (20^{th} step)') % First legend
%     legend('Deep (Constant)','Sup. (0^{th} step)','Sup. (10^{th} step)','Sup. (20^{th} step)') % Second legend
end


%% Fig 4 (b) - Mua or Musp versus nominal Mua step (curve-fitting)
fig_num = 250;
if exist('which_Fig','var') &&  ~ismember('4b', which_Fig)

else
    for Mua_or_Musp = 2 % For plotting  Mua (set to 1), Musp (set to 2), or both (set to 1:2)
%     for Mua_or_Musp = 1:2
    
        which_Exp = 1; % 1st, 2nd, or both experiments. For 2023 plotted 1st
%         which_Exp = 2;
%         which_Exp = [1 2]; % For both on one figure

        fig = figure(fig_num + 3*sum(which_Exp) + Mua_or_Musp); clf
        POS_1 = [70, 96.5, 850, 600];
        if set_dim == 1; fig.Position = POS_1; end
        
        j_chan_all = [3 5 7 9 11 13]; % For 2023 Publication
%         j_chan_all = [3 8]; % Which channels to show
%         j_chan_all = 7;
        
        j_steps_all = 1:21; % Which Mua steps to show (1 to 21 for 2023 Publication)
        
        lw_dtof = 2;
        ms2 = 15;
        
        % Define colors:
        temp_col = {'black','red',[130 0 255]/255,'magenta',[255 128 0]/255,'blue',[0 128 255]/255,'cyan','green','cyan','magenta'};
        temp_col_2 = {'black','blue','black','blue','black','blue','black','blue','black','blue','black','blue','black'};
        if length(temp_col) < 2 + length(j_steps_all) * 2
            for j = length(temp_col)+1 : 2 + length(j_steps_all) * 2
                temp_col{j} = 'red';
            end
        end
        
        % Plot:
        for j_chan = 1:length(j_chan_all)
            temp_X = mua_nom(:,j_chan_all(j_chan), 1);
            
            if ismember(1, which_Exp) == 1 % Deep, Exp. 1
                temp_musp = squeeze(LMA1_MuaMusp.Exp1_Deep_Red(j_steps_all,j_chan_all(j_chan),Mua_or_Musp)) * 10; % Exp. 1
                plot(temp_X(j_steps_all), temp_musp,'x--','Color',temp_col{1 + j_chan},'LineWidth',lw_dtof,'MarkerSize',ms2)
                hold on
            end
            if ismember(2, which_Exp) == 1 % Superficial, Exp. 2
                temp_musp = squeeze(LMA1_MuaMusp.Exp2_SupL_Red(j_steps_all,j_chan_all(j_chan),Mua_or_Musp)) * 10; % For Exp. 2
                plot(temp_X(j_steps_all), temp_musp,'o--','Color',temp_col_2{1 + j_chan},'LineWidth',lw_dtof,'MarkerSize',ms2)
                hold on
            end
        end
        
        temp = ''; % Legend
        for j = 1:length(j_chan_all)
            temp{j} = [num2str(wavelengths(j_chan_all(j))) ' nm'];
        end
        if Mua_or_Musp == 2
            if sum(which_Exp) < 3; legend(temp); end
            if ~ismember(2, which_Exp) == 1
                ylim([9.0 12])
            end
            set(gca,'YTick',9.0:1:13)
            ylabel('\mu''_s / cm^{-1}')
        else
            if sum(which_Exp) < 3; legend(temp,'location','southeast'); end
            ylim([0 0.25001])
            set(gca,'YTick',0.00:0.05:0.25)
            ylabel('\mu_a / cm^{-1}')
        end
        
        grid on
        xlim([0.025 0.25])
        set(gca,'XTick',[0.02 0.05:0.05:0.25])
        xlabel('Nominal \mu_a / cm^{-1}')
        set(gca,'FontSize',23)
        box on
    end
end


%% Fig 5 (b) - changes in moments
fig_num = 350; 
if exist('which_Fig','var') &&  ~ismember('5b', which_Fig)

else
    figure(fig_num); close(fig_num); % If using plotyyy, better to close the whole figure to clear it
    fig = figure(fig_num);
    POS_1 = [69, 96.5, 1250, 500];
    if set_dim == 1; fig.Position = POS_1; end
    
%     j_chan = 8; % For 2023 publication
    j_chan = which_chan;
    
    temp_X   = mua_nom(:,j_chan,1); % For Exp. #1
    temp_X_2 = mua_nom(:,j_chan,2); % For Exp. #2
    
    xlimits_s = [0.0497 0.24];
    
%     ylabels = {'\DeltaA','\Deltam_1','\DeltaV'};
    ylabels = {'','',''};
    
    % This wont plot anything, we use it for creating 3 y-axes and for legends
    [ax, ~] = plotyyy(-10, -10, -10, -10, -10, -10, ylabels, xlimits_s, fig_num);
    
    ylim_var(1,:) = [-1.5 1.5];
    ylim_var(2,:) = [-205 206]*0.90;
    ylim_var(3,:) = [-45 45];
    ylim_var(:,1) = ylim_var(:,1) * 1.4;
    
    lw = 2;
    col = {'red','blue','green'};
    
    for j_mom = 1:3
%         axes(ax(j_mom)); cla; hold on
        
        if j_mom == 1; ttt = 1; else; ttt = -1*10^3; end % Rescale m1 and V from +pico seconds to -nano seconds (See Publication 2023 Fig. 5b)
        
        % Exp. #1
        plot(ax(j_mom), temp_X, delMom.Exp1_TwoL_Blu(:,j_chan,j_mom)*ttt, 'o', 'Color', col{j_mom}, 'LineWidth',lw, 'MarkerSize', 10);
        plot(ax(j_mom), temp_X, delMom.Exp1_TwoL_Blu(:,j_chan,j_mom)*ttt, '--', 'Color', col{j_mom}, 'LineWidth',1, 'MarkerSize', 10);
    
        % Exp. #2
        if j_mom == 2; st = '+'; else; st = 'x'; end
        plot(ax(j_mom), temp_X_2, delMom.Exp2_TwoL_Blu(:,j_chan,j_mom)*ttt, st, 'Color', col{j_mom}, 'LineWidth',lw, 'MarkerSize', 14);
        plot(ax(j_mom), temp_X_2, delMom.Exp2_TwoL_Blu(:,j_chan,j_mom)*ttt, '--', 'Color', col{j_mom}, 'LineWidth',1, 'MarkerSize', 14);
        
        set(ax(j_mom),'YLim',ylim_var(j_mom,:))
        
        if j_mom==1
            val_YTick = -0.5*5:0.5:1.5;
            set(ax(j_mom),'YTick',val_YTick)
            set(ax(j_mom),'YColor', [235 0 0]/255)
            grid(ax(j_mom),'on')
            set(ax(j_mom),'GridAlpha',0.25)
            set(ax(j_mom),'XTick',[0.05:0.05:0.24 0.24])
        elseif j_mom==2
            set(ax(j_mom),'YTick',round( val_YTick / ylim_var(1,1) * ylim_var(2,1) ))
            set(ax(j_mom),'YColor', [0 0 204]/255)
            set(ax(j_mom),'XTick',[])
        elseif j_mom==3
            set(ax(j_mom),'YTick',round( val_YTick / ylim_var(1,1) * ylim_var(3,1) ))
            set(ax(j_mom),'YColor', [0 190 0]/255)
            set(ax(j_mom),'XTick',[])
        end

        set(ax(j_mom),'GridColor',[0 0 0])

        set(ax(j_mom),'FontSize',18)
%         set(gca,'YColor', 'Black')
        set(ax(j_mom),'Box','off')
    end
    title(['Spectral channel : ' num2str(round(wavelengths(j_chan))) ' nm'])
end


%% Fig 6 (a) and (b) - Determined Mua
fig_num = 450;
if exist('which_Fig','var') && ~(ismember('6a', which_Fig) || ismember('6b', which_Fig))
    
else
    if ~exist('which_Fig','var')
        plot_here = [1 2]; % [1 2] to plot '6a' and '6b'
    else
        plot_here = [0 0];
        if ismember('6a', which_Fig); plot_here(1) = 1; end
        if ismember('6b', which_Fig); plot_here(2) = 2; end
    end
%     plot_here = 1; % 1 is for plotting only '6a'.   2 is for plotting only '6b'.  [1 2] is for both

% for j_chan = 1:14 % Use this and pause on each itteration to check for all wavelengths
% for j_chan = 8 % 8th channel was plotted in 2023 publication
for j_chan = which_chan

    for J = 1:length(plot_here)
        if plot_here(J) > 0
            if plot_here(J) == 1 % 6a
                POS_1 = [50, 96.5, 1050, 600];
            else % 6b
                POS_1 = [250, 96.5, 1050, 600];
            end
            fig = figure( fig_num + J ); clf
            if set_dim == 1; fig.Position = POS_1; end

        
            ms = 15; ms2 = 11;
            lw = 2;
            
            c_gray = [155 155 155]/255; lw_gray = 1.3; % Unity line (y = x)
            plot([0 1], [0 1],'-','Color',c_gray,'LineWidth',lw_gray,'handlevisibility','off')
            hold on
        
            if plot_here(J) == 1 % 6a   Exp. 1
                temp_X   = mua_nom(:,j_chan,1); 
        
                plot(temp_X, LMA1_MuaMusp.Exp1_Deep_Red(:,j_chan,1) * 10,'+','MarkerSize',ms,'Color','black','LineWidth',lw) % 1. Hom.
        
                plot(temp_X, LMA1_MuaMusp.Exp1_TwoL_Blu(:,j_chan,1) * 10,'o','Color','green','MarkerSize',ms2,'LineWidth',lw) % 2. Hom.
        
                plot(temp_X, LMA2_Mua1Mua2.Exp1_TwoL_Blu(:,j_chan,1) * 10,'x','MarkerSize',ms,'Color','red','LineWidth',lw) % 3. Sup.
        
                plot(temp_X, LMA2_Mua1Mua2.Exp1_TwoL_Blu(:,j_chan,2) * 10,'x','MarkerSize',ms,'Color','blue','LineWidth',lw) % 4. Deep
        
                legend('Pos. #1  Homog.','Pos. #2  Homog.','Pos. #2  2-L. Sup.','Pos. #2  2-L. Deep','Location','NorthWest')
        
                xlabel('Nominal \mu_{a,Deep} / cm^{-1}')
        
            else % 6b   Exp. 2
                temp_X   = mua_nom(:,j_chan,2);
        
                plot(temp_X, LMA1_MuaMusp.Exp2_SupL_Red(:,j_chan,1) * 10,'+','MarkerSize',ms,'Color','black','LineWidth',lw) % 1. Hom.
        
                plot(temp_X, LMA1_MuaMusp.Exp2_TwoL_Blu(:,j_chan,1) * 10,'o','Color','green','MarkerSize',ms2,'LineWidth',lw) % 2. Hom.
        
                plot(temp_X, LMA2_Mua1Mua2.Exp2_TwoL_Blu(:,j_chan,1) * 10,'x','MarkerSize',ms,'Color','red','LineWidth',lw) % 3. Sup.
        
                plot(temp_X, LMA2_Mua1Mua2.Exp2_TwoL_Blu(:,j_chan,2) * 10,'x','MarkerSize',ms,'Color','blue','LineWidth',lw) % 4. Deep
        
                legend('Pos. #3  Homog.','Pos. #2  Homog.','Pos. #2  2-L. Sup.','Pos. #2  2-L. Deep','Location','NorthWest')
        
                xlabel('Nominal \mu_{a,Sup.} / cm^{-1}')
            end
            ylabel('\mu_{a} / cm^{-1}')
        
        %     axis equal
            box on
            grid on
            set(gca,'FontSize',24)
        
            xlim([0.0497 0.24])
            xlim([0.0488 0.24])
            set(gca,'XTick',[0.05:0.05:0.24 0.24])
        
            ylim([0.0497 0.24])
            ylim([0.0488 0.24])
            set(gca,'YTick',[0.05:0.05:0.2, 0.24])
        
            title(['Spectral channel : ' num2str(round(wavelengths(j_chan))) ' nm'])
        end
    end
end
end



%% Fig. Optional_1 - Check how all DTOFs look (allows to check the chosen region for background noise)
fig_num = 550;
if exist('which_Fig','var') &&  ~ismember('Optional_1', which_Fig)

else
    % Plot all DTOFs for one of 4 sets (2 experiments x 2 source-detector pairs)
    wh_show = 4; % Set value from 1 to 4
    
    figure(fig_num); clf
    
%     j_chan = 8;
    j_chan = which_chan;
    
    st_dtof = '-'; 
%     st_dtof = '.'; 
%     st_dtof = '--';

    do_norm = 'norm';
%     do_norm = 0;
    
    for j_MuaStep = 1:21
        if j_MuaStep == 1; do_raw = 'raw'; else; do_raw = 0; end % < Display unfiltered DTOF for first cycle
    
        if wh_show == 1     % DTOFs for Red detector, for Exp. #1  (Measuring on Deep, increasing Deep)
            [ ~, ~ ] = DTOF_filter( squeeze(mean(DTOF.Exp1_Deep_Red(j_MuaStep,j_chan,:),1)), Time_ns, bg_Red, 2, 0, 0, {fig_num,do_norm,do_raw,st_dtof});
    
        elseif wh_show == 2 % DTOFs for Red detector, for Exp. #2  (Measuring on Sup, increasing Sup)
            [ ~, ~ ] = DTOF_filter( squeeze(mean(DTOF.Exp2_SupL_Red(j_MuaStep,j_chan,:),1)), Time_ns, bg_Red, 2, 0, 0, {fig_num,do_norm,do_raw,st_dtof});
    
        elseif wh_show == 3 % DTOFs for Blu detector, for Exp. #1  (Measuring on Two-layered, increasing Deep)
            [ ~, ~ ] = DTOF_filter( squeeze(mean(DTOF.Exp1_TwoL_Blu(j_MuaStep,j_chan,:),1)), Time_ns, bg_Blu, 2, 0, 0, {fig_num,do_norm,do_raw,st_dtof});
    
        elseif wh_show == 4 % DTOFs for Blu detector, for Exp. #2  (Measuring on Two-layered, increasing Sup)
            [ ~, ~ ] = DTOF_filter( squeeze(mean(DTOF.Exp2_TwoL_Blu(j_MuaStep,j_chan,:),1)), Time_ns, bg_Blu, 2, 0, 0, {fig_num,do_norm,do_raw,st_dtof});
        end
        
        if isequal(do_norm,'norm')
            ylim([1e-5 1])
        end
        xlim([0 11.5])
    end
    legend('Step 1 Unfiltered','Step 1 filtered','Step 2 filtered','Step 3 filtered','...','FontSize',20)
    
    title(['Spectral channel : ' num2str(round(wavelengths(j_chan))) ' nm'])
end


%% Fig. Optional_2  Check how IRFs compare for the two experiments (performed one after the other on the same day)
fig_num = 560;
if exist('which_Fig','var') &&  ~ismember('Optional_2', which_Fig)

else
    figure(fig_num); clf
    
%     j_chan = 8;
    j_chan = which_chan;
    
    do_shift = 1; % To shift IRF by 0.4 ONLY FOR DISPLAYING PURPOSES HERE

    if do_shift == 0
        % IRFs for RED detector
        [ ~, ~ ] = DTOF_filter( squeeze(IRF.Exp1_Red(:,j_chan)), Time_ns, bg_Red, 2, [0.01 0.001], 0, {fig_num, 'blue','norm'}); % At start of Exp. 1
        [ ~, ~ ] = DTOF_filter( squeeze(IRF.Exp2_Red(:,j_chan)), Time_ns, bg_Blu, 2, [0.01 0.001], 0, {fig_num, 'red','norm'});  % At start of Exp. 2

        % IRFs for BLUE detector
        [ ~, ~ ] = DTOF_filter( squeeze(IRF.Exp1_Blu(:,j_chan)), Time_ns, bg_Blu, 2, [0.01 0.001], 0, {fig_num, 'green','norm'}); % At start of Exp. 1
        [ ~, ~ ] = DTOF_filter( squeeze(IRF.Exp2_Blu(:,j_chan)), Time_ns, bg_Blu, 2, [0.01 0.001], 0, {fig_num, 'black','norm'}); % At start of Exp. 2
    else
        % IRFs for RED detector
        [ ~, ~ ] = DTOF_filter( squeeze(IRF.Exp1_Red(:,j_chan)), Time_ns, bg_Red, 2, [0.01 0.001], 0, {fig_num, 'blue','norm'}); % At start of Exp. 1
        [ ~, ~ ] = DTOF_filter( squeeze(IRF.Exp2_Red(:,j_chan)), Time_ns + 0.4*(Time_ns(2)-Time_ns(1)), bg_Blu, 2, [0.01 0.001], 0, {fig_num, 'red','norm'}); % At start of Exp. 2
        
        % IRFs for BLUE detector
        [ ~, ~ ] = DTOF_filter( squeeze(IRF.Exp1_Blu(:,j_chan)), Time_ns, bg_Blu, 2, [0.01 0.001], 0, {fig_num, 'green','norm'}); % At start of Exp. 1
        [ ~, ~ ] = DTOF_filter( squeeze(IRF.Exp2_Blu(:,j_chan)), Time_ns + 0.2*(Time_ns(2)-Time_ns(1)), bg_Blu, 2, [0.01 0.001], 0, {fig_num, 'black','norm'}); % At start of Exp. 2
    end
    ylim([1e-2 1])
    xlim([0.95 2.1])
    
    if do_shift == 1
        title(['MANUALLY SHIFTED IRFs BY 0.4 time channels for displaying purposes here.  Spectral channel : ' num2str(round(wavelengths(j_chan))) ' nm'])
    else
        title(['Raw IRFs.   Spectral channel : ' num2str(round(wavelengths(j_chan))) ' nm'])
    end
end


%% Fig. Optional_3 - Determined Mua1 and Mua2 using Curve-fitting (single DTOF)
fig_num = 570;
if exist('which_Fig','var') &&  ~ismember('Optional_3', which_Fig)

else
    POS_1 = [50, 96.5, 1050, 600];
    fig = figure(fig_num); clf
    if set_dim == 1; fig.Position = POS_1; end

    ms = 15;
    lw = 2;
    
%     j_chan = 8;
    j_chan = which_chan;
    
    c_gray = [155 155 155]/255; lw_gray = 1.3;
    plot([0 1], [0 1],'-','Color',c_gray,'LineWidth',lw_gray,'handlevisibility','off')
    hold on

    temp_X   = mua_nom(:,j_chan,1); % Exp. 1
    plot(temp_X, LMA1_Optional.Exp1_TwoL_Blu(:,j_chan,1) * 10,'x','MarkerSize',ms,'Color','red','LineWidth',lw) % Exp. 1. Sup

    plot(temp_X, LMA1_Optional.Exp1_TwoL_Blu(:,j_chan,2) * 10,'x','MarkerSize',ms,'Color','blue','LineWidth',lw) % Exp. 1. Deep

    plot(temp_X, LMA1_MuaMusp.Exp1_Deep_Red(:,j_chan,1) * 10,'+','MarkerSize',ms,'Color','black','LineWidth',lw) % Exp. 1. On a single layer


    temp_X   = mua_nom(:,j_chan,2); % Exp. 2
    plot(temp_X, LMA1_Optional.Exp2_TwoL_Blu(:,j_chan,1) * 10,'+','MarkerSize',ms,'Color','magenta','LineWidth',lw) % Exp. 2. Sup

    plot(temp_X, LMA1_Optional.Exp2_TwoL_Blu(:,j_chan,2) * 10,'+','MarkerSize',ms,'Color','cyan','LineWidth',lw) % Exp. 2. Deep

    plot(temp_X, LMA1_MuaMusp.Exp2_SupL_Red(:,j_chan,1) * 10,'x','MarkerSize',ms,'Color','green','LineWidth',lw) % Exp. 2. On a single layer

    legend('Exp.1 Sup','Exp.1 Deep','Exp.1 Single Layer','Exp.2 Sup.','Exp.2 Deep','Exp.2 Single Layer','Location','NorthWest')

    xlabel('Nominal \mu_a / cm^{-1}')
    ylabel('\mu_{a} / cm^{-1}')

%     axis equal
    box on
    grid on
    set(gca,'FontSize',24)

    xlim([0.0499 0.24])
    set(gca,'XTick',[0.05:0.05:0.24 0.24])

    ylim([0.044 0.24])
    set(gca,'YTick',[0.05:0.05:0.2, 0.24])

    title(['Spectral channel : ' num2str(round(wavelengths(j_chan))) ' nm'])
end
end