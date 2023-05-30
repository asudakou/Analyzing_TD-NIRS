% Written by Aleh Sudakou (https://www.researchgate.net/profile/Aleh-Sudakou)
% Nałęcz Institute of Biocybernetics and Biomedical Engineering, Polish Academy of Sciences 
% This is one of the codes shared on  https://github.com/asudakou/Analyzing_TD-NIRS
% Last updated: 20 May 2023


%% Description:
% Testing the uniquness of delN, delm1, and delV for delMua1 and delMua2.
% This script calculates changes in moments for all combinations of delMua1 and delMua2
% Then plots a figure, showing that for a combination of delMua1 and
% delMua2, there exists a unique set of changes in any 2 moments
%
% The results are presented in Fig. 2 (a) in the 2023 Publication.
 

%% Main

clear; % clc

% Folder where to save the results, since this can be a lengthy script to run:
folder_calc = 'C:\Users\asudakou\Documents\01 Work mine\MatLab codes\Wrote by me\TR analyses\Two-layered analysis\Aleh LMA fitting\Open-access\DATA calculated';

file_name = 'Other_uniqueness_vr1'; % used L_known = 12  and  [0.013 0.013 0.013 1 1 1 L_known 40]

do_load = 1;

if do_load == 1
    load(file_name)
    
else
    do_save = 1;

%     cut_lim_mom = [0.25 0.03];
    cut_lim_mom = [0.01 0.01];

    rho = 30;

    n = 1.33;

    L_known = 12;

    opt_prop_base = [0.013 0.013 0.013 1 1 1 L_known 40]; % [Mua1 Mua2 Mua3 Musp1 Musp2 Musp3 L1 L2]

    [R_base, time_ns] = DTOF_generate_Liemert(opt_prop_base,  n, rho, -1);

    [ ~, cut_ind_mom ] = DTOF_filter( R_base, time_ns, 0, 0, cut_lim_mom, 1);

    Mom_base = DTOF_CentralMom(time_ns(cut_ind_mom(1):cut_ind_mom(2)), R_base(cut_ind_mom(1):cut_ind_mom(2)));

    % Range of changes in Mua1 and Mua2:
    delMua1 = -0.01:0.001:0.01;
    delMua2 = -0.01:0.001:0.01;

    % Results:
    delN = zeros(length(delMua1), length(delMua2));
    delm1 = zeros(size(delN));
    delV = zeros(size(delN));
    delm3c = zeros(size(delN));

    % The long loop, where we calculate changes in moments for changes in Mua
    opt_prop_base_temp = opt_prop_base;
    for j_1 = 1:length(delMua1)
        opt_prop_base_temp(1) = opt_prop_base(1) + delMua1(j_1); % [Mua1 Mua2 Mua3 Musp1 Musp2 Musp3 L1 L2]

        for j_2 = 1:length(delMua2)
            opt_prop_base_temp(2:3) = opt_prop_base(2) + delMua2(j_2); % [Mua1 Mua2 Mua3 Musp1 Musp2 Musp3 L1 L2]

            % New DTOF
            [R_after, ~] = DTOF_generate_Liemert(opt_prop_base_temp,  n, rho, -1);

            Mom_after = DTOF_CentralMom(time_ns(cut_ind_mom(1):cut_ind_mom(2)), R_after(cut_ind_mom(1):cut_ind_mom(2)));

            temp = DTOF_DelMom(Mom_base, Mom_after);
            delN(j_1, j_2) = temp(1);
            delm1(j_1, j_2) = temp(2);
            delV(j_1, j_2) = temp(3);
            delm3c(j_1, j_2) = temp(4);
        end
    end

    if do_save == 1
        save([folder_calc '\' file_name], 'delN','delm1','delV','delm3c','delMua1','delMua2')

        disp([char(datetime('now','Format','HH:mm:ss')) '  Saved into variable ''' file_name ''''])
    end
end


%% Plot    Note, we can also plot for the 3rd central moment

figure(5); clf

lw = 2;

% do_show_text = 'off';
do_show_text = 'on'; % Uncomment to see texts

if isequal(do_show_text,'on')
    text_font = 20;
    text_spacing = 144*3;
end

lines_spacing = 0.8; % Total number of photons (Ntot)
[C, h] = contour(delMua1*10, delMua2*10, delN', -lines_spacing*4:lines_spacing:lines_spacing*4, ':', 'LineWidth',lw,'ShowText',do_show_text,'Color','red');
if isequal(do_show_text,'on')
    clabel(C,h,'FontSize',text_font,'Color','black')
    h.LabelSpacing = text_spacing;
end
hold on

lines_spacing = 0.10 * 10^3; % Mean time of flight (m1)
[C, h] = contour(delMua1*10, delMua2*10, delm1' * 10^3, -lines_spacing*4:lines_spacing:lines_spacing*4, '--', 'LineWidth',lw,'ShowText',do_show_text,'Color','blue');
if isequal(do_show_text,'on')
    clabel(C,h,'FontSize',text_font,'Color','black')
    h.LabelSpacing = text_spacing;
end

lines_spacing = 0.03 * 10^3; % Variance (V)
[C, h] = contour(delMua1*10, delMua2*10, delV' * 10^3, -lines_spacing*4:lines_spacing:lines_spacing*4, '-', 'LineWidth',lw,'ShowText',do_show_text,'Color','green');
if isequal(do_show_text,'on')
    clabel(C,h,'FontSize',text_font,'Color','black')
    h.LabelSpacing = text_spacing;
end

% lines_spacing = 0.005 * 10^3; % Third central moment (m3c)
% [C, h] = contour(delMua1*10, delMua2*10, delm3c' * 10^3, -lines_spacing*4:lines_spacing:lines_spacing*4, '--', 'LineWidth',lw,'ShowText',do_show_text,'Color','magenta');
% if isequal(do_show_text,'on')
%     clabel(C,h,'FontSize',text_font,'Color','black')
%     h.LabelSpacing = text_spacing;
% end


xlabel('\Delta\mu_{a,Sup} / cm^{-1}')
ylabel('\Delta\mu_{a,Deep} / cm^{-1}')

% xlim([min(delMua1) max(delMua1)]*10)
% ylim([min(delMua2) max(delMua2)]*10)
xlim([-0.01 0.0100001]*10)
ylim([-0.01 0.0100001]*10)

set(gca,'XTick',(-0.01:0.005:0.015)*10)
set(gca,'YTick',(-0.01:0.005:0.015)*10)

axis equal
grid on
set(gca,'FontSize',24)