% Written by Aleh Sudakou (https://www.researchgate.net/profile/Aleh-Sudakou)
% Nałęcz Institute of Biocybernetics and Biomedical Engineering, Polish Academy of Sciences 
% This is one of the codes shared on  https://github.com/asudakou/Analyzing_TD-NIRS
% Last updated: 20 May 2023


%% Description:
% Calculate Chi^2 (similar to a residual) for all combinations of delMua1 and delMua2
% 
% The results are presented in Fig. 2 (b, c) in the 2023 Publication.


%% Main

clear; % clc

% Folder where to save the results, since this can be a lengthy script to run:
folder_calc = 'C:\Users\asudakou\Documents\01 Work mine\MatLab codes\Wrote by me\TR analyses\Two-layered analysis\Aleh LMA fitting\Open-access\DATA calculated';

file_name = 'Other_Chi2_vr1'; % used L_known = 12  and  [0.013 0.013 0.013 1 1 1 L_known 40]

do_load = 1;

if do_load == 1
    load(file_name)

else
    do_save = 1;

    cut_lim_mom = [0.25 0.03];
%     cut_lim_mom = [0.01 0.01];

    rho = 30;

    n = 1.33;

    L_known = 12; 

    opt_prop_base = [0.013 0.013 0.013 1 1 1 L_known 40]; % [Mua1 Mua2 Mua3 Musp1 Musp2 Musp3 L1 L2]

    [R_base, time_ns] = DTOF_generate_Liemert(opt_prop_base,  n, rho, -1);

    [ ~, cut_ind_mom ] = DTOF_filter( R_base, time_ns, 0, 0, cut_lim_mom, 1);

    Mom_base = DTOF_CentralMom(time_ns(cut_ind_mom(1):cut_ind_mom(2)), R_base(cut_ind_mom(1):cut_ind_mom(2)));

    std_base = 1; % Will use weights (uncertainties of moments) to make Chi unitless
    std_base(2) = ( Mom_base(3)  );
    std_base(3) = ( (Mom_base(5) - Mom_base(3)^2)  );

    % Range of changes in Mua1 and Mua2:
    delMua1 = -0.01:0.0001:0.01;
    delMua2 = -0.01:0.0001:0.01;
    
    Chi2_withN = zeros(size(delMua1,1), size(delMua2,2));
    Chi2_withoutN = zeros(size(delMua1,1), size(delMua2,2));

    opt_prop_base_temp = opt_prop_base;
    for j_1 = 1:length(delMua1)
        opt_prop_base_temp(1) = opt_prop_base(1) + delMua1(j_1); % [Mua1 Mua2 Mua3 Musp1 Musp2 Musp3 L1 L2]

        for j_2 = 1:length(delMua2)
            opt_prop_base_temp(2:3) = opt_prop_base(2) + delMua2(j_2); % [Mua1 Mua2 Mua3 Musp1 Musp2 Musp3 L1 L2]

            [R_after, ~] = DTOF_generate_Liemert(opt_prop_base_temp,  n, rho, -1);

            Mom_after = DTOF_CentralMom(time_ns(cut_ind_mom(1):cut_ind_mom(2)), R_after(cut_ind_mom(1):cut_ind_mom(2)));

            delMom = DTOF_DelMom(Mom_base, Mom_after);

            delMom(1:3) = delMom(1:3) ./ sqrt(std_base');

            % Calculate Chi^2:
            Chi2_withN(j_1, j_2) = sum(  delMom(1).^2 + delMom(2).^2 + delMom(3).^2  );

            Chi2_withoutN(j_1, j_2) = sum(  delMom(2).^2 + delMom(3).^2  );

        end
    end

    if do_save == 1
        save([folder_calc '\' file_name], 'Chi2_withN','Chi2_withoutN','delMua1','delMua2')

        disp([char(datetime('now','Format','HH:mm:ss')) '  Saved into variable ''' file_name ''''])
    end
end


%% Plot the Chi^2 distribution for LMA 2

figure(5); clf

h = pcolor(delMua1 * 10, delMua2 * 10, Chi2_withN');
h.EdgeColor = 'none';
hold on

h = colorbar;
clim([0 0.5])
colormap('jet');
% colormap('hsv');

ylabel(h, '2 / N_{tot}   \chi^2')
xlabel('\Delta\mu_{a,Sup} / cm^{-1}')
ylabel('\Delta\mu_{a,Deep} / cm^{-1}')

plot([1 1] * 0 * 10, [1 1] * 0 * 10, '+','MarkerSize',15,'LineWidth',4,'Color','White')
% set(gca, 'XTick', 0.1:0.02:0.16)

set(gca,'FontSize',24)

axis square


% Repeat the above, but for Chi2_withoutN
figure(6); clf

h = pcolor(delMua1 * 10, delMua2 * 10, Chi2_withoutN');
h.EdgeColor = 'none';
hold on

h = colorbar;
clim([0 0.05])
colormap('jet');
% colormap('hsv');

ylabel(h, '2 / N_{tot}   \chi^2')
xlabel('\Delta\mu_{a,Sup} / cm^{-1}')
ylabel('\Delta\mu_{a,Deep} / cm^{-1}')

plot([1 1]*0*10, [1 1]*0*10, '+','MarkerSize',15,'LineWidth',4,'Color','White')
% set(gca, 'XTick', 0.1:0.02:0.16)

set(gca,'FontSize',24)

axis square