% Written by Aleh Sudakou (https://www.researchgate.net/profile/Aleh-Sudakou)
% Nałęcz Institute of Biocybernetics and Biomedical Engineering, Polish Academy of Sciences 
% This is one of the codes shared on  https://github.com/asudakou/Analyzing_TD-NIRS
% Last updated: 20 May 2023


%% Description:
% This script does data analysis for two experiments with ink.
% For full details, please refer to the publication.
%
% "Ink_plotting.m" plots all figures presented in Publication 2023
% 
% In Exp1 we increased Mua in the Deep layer
% In Exp2 we increased Mua in the Superficial layer


%% Step 1. Load data and define empty variables
% These are all of the variables that will contain results or any relevant data

clear; % clc;

% Folder where to save/load the calculated results, as some calculations can take some time:
folder_calc = 'C:\Users\asudakou\Documents\01 Work mine\MatLab codes\Wrote by me\TR analyses\Two-layered analysis\Aleh LMA fitting\Open-access\DATA calculated';

% For Step 1 - Load the measured DTOF for 2 experiments with ink:
           load('data_Ink_Pub2023')
           DTOF;
           DTOF_SameND;
           DTOF_ConstLayer; % At start and end of measurement
           IRF;
           Time_ns;
           % Additional variable: filtered and shifted IRF:
           IRF_shifted.Exp1_Blu = zeros(1024,16);
           IRF_shifted.Exp1_Red = zeros(size(IRF_shifted.Exp1_Blu)); % Blu and Red correspond to two different source-detector pairs
           IRF_shifted.Exp2_Blu = zeros(size(IRF_shifted.Exp1_Blu));
           IRF_shifted.Exp2_Red = zeros(size(IRF_shifted.Exp1_Blu));

% For Step 2 - Perform Homogeneous fitting (LMA 1) to obtain Mua and Musp
           LMA1_MuaMusp = [];
           LMA1_MuaMusp.Exp1_TwoL_Blu = zeros(21,14,2);
           LMA1_MuaMusp.Exp1_Deep_Red = zeros(21,14,2);
           LMA1_MuaMusp.Exp2_TwoL_Blu = zeros(21,14,2);
           LMA1_MuaMusp.Exp2_SupL_Red = zeros(21,14,2);
           LMA1_MuaMusp.ConstLayer_Exp1_SupL_Blu = zeros(2,14,2);
           LMA1_MuaMusp.ConstLayer_Exp2_Deep_Blu = zeros(2,14,2);

% For Step 3 - Calculate moments and perform 2-layered fitting (LMA 2) to determine delMua1 and delMua2
           Mom = [];
           Mom.Exp1_TwoL_Blu = zeros(21,14,5);
           Mom.Exp1_Deep_Red = zeros(21,14,5);
           Mom.Exp2_TwoL_Blu = zeros(21,14,5);
           Mom.Exp2_SupL_Red = zeros(21,14,5);
           delMom = [];
           delMom.Exp1_TwoL_Blu = zeros(21,14,4);
           delMom.Exp1_Deep_Red = zeros(21,14,4);
           delMom.Exp2_TwoL_Blu = zeros(21,14,4);
           delMom.Exp2_SupL_Red = zeros(21,14,4);
           LMA2_Mua1Mua2 = [];
           LMA2_Mua1Mua2.Exp1_TwoL_Blu = zeros(21,14,2);
           LMA2_Mua1Mua2.Exp2_TwoL_Blu = zeros(21,14,2);

% For Step 4 (Optional, not in publication) - Use LMA 1 to determine Mua1 and Mua2, instead of Mua and Musp
           LMA1_Optional = [];


% Define all constants and filter IRFs:

% set_dim = 1; % Setting this to 1 will make the figure sizes the same as we used for Publication 2023
set_dim = 0;

wavelengths = 680:12.5:868; % 16 wavelengths, unit nm

mua_nom = Ink_Nominal_mua(wavelengths); % Nominal Mua values based on concentration of ink - used only by "Ink_plotting.m" (x-axis)

rho = 30; % Units: mm. Source-detector distance, it was always set to 30 mm

offset = 0.06/physconst('Lightspeed')*1e9; % ns. Time offset of IRF due the 6 cm distance between fibers when measuring IRF

n = 1.33; % Refractive index. Assuming it is the same as water because we measured on liquid phantoms where main component was water

bg_Red = [0.49 0.92]; % Region for subtracting background signal
bg_Blu = [0.49 0.92]; % For the second detector
bg_irf = [7.3 8.3]; % For IRF

% Filter IRFs:
for j_chan = 1:16
    % Remove background and cut at 0% on both sides of the maximum, although this filtering of IRF makes little difference to the results of data analysis
    [IRF_shifted.Exp1_Blu(:,j_chan), ~ ] = DTOF_filter( squeeze(IRF.Exp1_Blu(:,j_chan)), Time_ns, bg_irf, 2, [0 0], 1);
    [IRF_shifted.Exp1_Red(:,j_chan), ~ ] = DTOF_filter( squeeze(IRF.Exp1_Red(:,j_chan)), Time_ns, bg_irf, 2, [0 0], 1);
    [IRF_shifted.Exp2_Blu(:,j_chan), ~ ] = DTOF_filter( squeeze(IRF.Exp2_Blu(:,j_chan)), Time_ns, bg_irf, 2, [0 0], 1);
    [IRF_shifted.Exp2_Red(:,j_chan), ~ ] = DTOF_filter( squeeze(IRF.Exp2_Red(:,j_chan)), Time_ns, bg_irf, 2, [0 0], 1);

    % Shift by offset
    IRF_shifted.Exp1_Blu(:,j_chan) = interp1(Time_ns - offset, IRF_shifted.Exp1_Blu(:,j_chan), Time_ns);
    IRF_shifted.Exp1_Red(:,j_chan) = interp1(Time_ns - offset, IRF_shifted.Exp1_Red(:,j_chan), Time_ns);
    IRF_shifted.Exp2_Blu(:,j_chan) = interp1(Time_ns - offset, IRF_shifted.Exp2_Blu(:,j_chan), Time_ns);
    IRF_shifted.Exp2_Red(:,j_chan) = interp1(Time_ns - offset, IRF_shifted.Exp2_Red(:,j_chan), Time_ns);
end

disp([char(datetime('now','Format','HH:mm:ss')) '  Step 1. Loaded measured data. And defined empty variables for all results'])


% PLOT:
which_chan = 8; % For "Ink_plotting.m"
show_chan = 8; % 8th is the channel shown in Publication 2023
% show_chan = 5; % To show for any other channel
show_fig = {'5a','Optional_1','Optional_2'};
Ink_plotting(DTOF,Time_ns,IRF,delMom,LMA1_MuaMusp,LMA2_Mua1Mua2,LMA1_Optional, set_dim, show_chan, show_fig);


%% Step 2.  Perform Homogeneous fitting for all DTOF

% For saving/loading the calculated data:
file_name = 'Ink_LMA1_MuaMusp_vr1';  % cut_lim_fit = [0.85 0.01]
% file_name = 'Ink_LMA1_MuaMusp_vr2';

do_load = 1; % 1 == load previously calculated results. 0 == Do calculation
% do_load = 0; 

do_save = 1; % 1 == Save result in 'folder_calc'. Will use 'file_name' for the saved file
% do_save = 0;

if do_load == 1
%     load([folder_calc '\' file_name])
    load(file_name)

else
    cut_lim_fit = [0.85 0.01]; % Region for fitting

    OptProp_all = [-7 -8 -8 -7 -8 -8 40 40]; % [Mua1 Mua2 Mua3 Musp1 Musp2 Musp3 L1 L2]

    j_chan_all = 1:14; % For which spectral channels to do the calculation


    % PART A:  Perform Homogeneous fitting for the constant layer
    for j_chan = j_chan_all
        % Exp. 1,  On Superficial at start and at end of experiment
        for j_cycle = 1:2
            [ dtof, cut_ind_fit] = DTOF_filter( squeeze(DTOF_ConstLayer.Exp1_SupL_Blu(j_cycle,j_chan,:)), Time_ns, bg_Blu, 2, cut_lim_fit, 1);
            
            [X_found, ~] = LMA_1_FittingDTOF(OptProp_all, n, rho, dtof, Time_ns, cut_ind_fit, IRF_shifted.Exp1_Blu(:,j_chan));
            LMA1_MuaMusp.ConstLayer_Exp1_SupL_Blu(j_cycle,j_chan,1:2) = X_found(1:2);
        end
        % Exp. 2,  On Deep at start and at end of experiment
        counter = 1;
        for j_cycle = 1:2
            [ dtof, cut_ind_fit] = DTOF_filter( squeeze(DTOF_ConstLayer.Exp2_Deep_Blu(j_cycle,j_chan,:)), Time_ns, bg_Blu, 2, cut_lim_fit, 1);

            [X_found, ~] = LMA_1_FittingDTOF(OptProp_all, n, rho, dtof, Time_ns, cut_ind_fit, IRF_shifted.Exp2_Blu(:,j_chan));
            LMA1_MuaMusp.ConstLayer_Exp2_Deep_Blu(counter,j_chan,1:2) = X_found(1:2); counter = counter + 1;
        end
    end


    % PART B:  Perform Homogeneous fitting (LMA 1) for all 21 steps of Mua
    j_cycle_all = 1:21; 

    for j_chan = j_chan_all
        % Exp. 1, on Two-Layered, Blu detector
        for j_cycle = j_cycle_all
            [ dtof, cut_ind_fit] = DTOF_filter( squeeze(DTOF.Exp1_TwoL_Blu(j_cycle,j_chan,:)), Time_ns, bg_Blu, 2, cut_lim_fit, 1);

            [X_found, ~] = LMA_1_FittingDTOF(OptProp_all, n, rho, dtof, Time_ns, cut_ind_fit, IRF_shifted.Exp1_Blu(:,j_chan));
            LMA1_MuaMusp.Exp1_TwoL_Blu(j_cycle,j_chan,1:2) = X_found(1:2);
        end

        % The following 3 loops are repetitions of the above loop but for the other sets of DTOF
        
        for j_cycle = j_cycle_all % Exp. 1, on Deep, Red detector
            [ dtof, cut_ind_fit] = DTOF_filter( squeeze(DTOF.Exp1_Deep_Red(j_cycle,j_chan,:)), Time_ns, bg_Red, 2, cut_lim_fit, 1);
            
            [X_found, ~] = LMA_1_FittingDTOF(OptProp_all, n, rho, dtof, Time_ns, cut_ind_fit, IRF_shifted.Exp1_Red(:,j_chan));
            LMA1_MuaMusp.Exp1_Deep_Red(j_cycle,j_chan,1:2) = X_found(1:2);
        end

        for j_cycle = j_cycle_all % Exp. 2, on Two-Layered, Blu detector
            [ dtof, cut_ind_fit] = DTOF_filter( squeeze(DTOF.Exp2_TwoL_Blu(j_cycle,j_chan,:)), Time_ns, bg_Blu, 2, cut_lim_fit, 1);

            [X_found, ~] = LMA_1_FittingDTOF(OptProp_all, n, rho, dtof, Time_ns, cut_ind_fit, IRF_shifted.Exp2_Blu(:,j_chan));
            LMA1_MuaMusp.Exp2_TwoL_Blu(j_cycle,j_chan,1:2) = X_found(1:2);
        end

        for j_cycle = j_cycle_all % Exp. 2, on Superficial, Red detector
            [ dtof, cut_ind_fit] = DTOF_filter( squeeze(DTOF.Exp2_SupL_Red(j_cycle,j_chan,:)), Time_ns, bg_Red, 2, cut_lim_fit, 1);
            [X_found, ~] = LMA_1_FittingDTOF(OptProp_all, n, rho, dtof, Time_ns, cut_ind_fit, IRF_shifted.Exp2_Red(:,j_chan));
            LMA1_MuaMusp.Exp2_SupL_Red(j_cycle,j_chan,1:2) = X_found(1:2);
        end

        disp([char(datetime('now','Format','HH:mm:ss')) '  Finished for channel ' num2str(j_chan)])
    end

    % Saving:
    if do_save == 1
        save([folder_calc '\' file_name], 'LMA1_MuaMusp','cut_lim_fit')

        disp([char(datetime('now','Format','HH:mm:ss')) '  Saved ''LMA1_MuaMusp'' into variable ''' file_name ''''])
    end
end
if do_load == 0
    disp([char(datetime('now','Format','HH:mm:ss')) '  Step 2. Calculated Mua and Musp (homogeneous model) for all DTOF'])
else
    disp([char(datetime('now','Format','HH:mm:ss')) '  Step 2. LOADED previously determined Mua and Musp (homogeneous model) for all DTOF'])
end


% PLOT:
show_chan = 8; % 8th is the channel shown in Publication 2023
show_fig = {'4a','4b'};
Ink_plotting(DTOF,Time_ns,IRF,delMom,LMA1_MuaMusp,LMA2_Mua1Mua2,LMA1_Optional, set_dim, show_chan, show_fig);


%% Step 3.  Calculate moments, changes in moments, and perform 2-layered fitting for delMoments

% For saving/loading the calculated data:
file_name = 'Ink_LMA2_Mua1Mua2_vr1'; % cut_lim_mom = [0.25 0.03];
% file_name = 'Ink_LMA2_Mua1Mua2_vr2';

do_load = 1; % 1 == load previously calculated results. 0 == Do calculation
% do_load = 0; 

do_save = 1; % 1 == Save result in 'folder_calc'. Will use 'file_name' for the saved file
% do_save = 0;

if do_load == 1
%     load([folder_calc '\' file_name])
    load(file_name)

else
    cut_lim_mom = [0.25 0.03]; % Used this for Publication 2023

    c_base = [11 11]; % Which Mua step to use as Baseline for 1st and 2nd experiments

    L_known = 14.5;

    delOptProp_all = [-7 -7 -8 0 0 0 L_known 40]; % [delMua1 delMua2 delMua3 delMusp1 delMusp2 delMusp3 L1 L2]

    Options_LMA = [];

    OptProp_base = zeros(6,1);

    j_cycle_all = 1:21;

    j_chan_all = 1:14;
%     j_chan_all = 8;

    for j_chan = j_chan_all
        
        % ----------- For Exp. #1 (start) -----------
        % Step 1.  Set baseline optical properties of both layers
        OptProp_base(1) = LMA1_MuaMusp.ConstLayer_Exp1_SupL_Blu(1,j_chan,1); % Mua  of Superficial layer
        OptProp_base(4) = LMA1_MuaMusp.ConstLayer_Exp1_SupL_Blu(1,j_chan,2); % Musp of Superficial layer
        OptProp_base(2:3) = LMA1_MuaMusp.Exp1_Deep_Red(c_base(1),j_chan,1); % Mua  of Deep layer
        OptProp_base(5:6) = LMA1_MuaMusp.Exp1_Deep_Red(c_base(1),j_chan,2); % Musp of Deep layer

        [ ~, cut_ind_mom]     = DTOF_filter( squeeze(DTOF.Exp1_TwoL_Blu(c_base(1),j_chan,:)), Time_ns, bg_Blu, 2, cut_lim_mom, 1); % Region of time-channels for calculating moments
        [ ~, cut_ind_mom_red] = DTOF_filter( squeeze(DTOF.Exp1_Deep_Red(c_base(1),j_chan,:)), Time_ns, bg_Red, 2, cut_lim_mom, 1); % For the other source-detector (red)
        
        % Step 2.  Calculate moments 
        correct_N = 1;
        correct_N_red = 1;
        for j_cycle = 1:21
            % For Blu detector on two-layered window:
            [dtof,~] = DTOF_filter( squeeze(DTOF.Exp1_TwoL_Blu(j_cycle,j_chan,:)), Time_ns, bg_Blu, 2, cut_ind_mom, 3);
            Mom.Exp1_TwoL_Blu(j_cycle,j_chan,:) = DTOF_CentralMom(Time_ns, dtof);
            
            % Correct Ntot if changed laser intensity - this doesn't affect any calculations since we chose not to use Ntot
            if sum(DTOF_SameND.Exp1_TwoL_Blu(j_cycle,1,:)) ~= 0
                [dtof, ~] = DTOF_filter( squeeze(DTOF_SameND.Exp1_TwoL_Blu(j_cycle,j_chan,:)), Time_ns, bg_Blu, 2, cut_ind_mom, 3);
                mom_temp = DTOF_CentralMom(Time_ns, dtof);
                
                correct_N = correct_N * (mom_temp(1) / Mom.Exp1_TwoL_Blu(j_cycle,j_chan,1));
            end
            Mom.Exp1_TwoL_Blu(j_cycle,j_chan,1) = Mom.Exp1_TwoL_Blu(j_cycle,j_chan,1) * correct_N;

            % Repeat for Red detector on deep layer:
            [dtof,~] = DTOF_filter( squeeze(DTOF.Exp1_Deep_Red(j_cycle,j_chan,:)), Time_ns, bg_Red, 2, cut_ind_mom_red, 3);
            Mom.Exp1_Deep_Red(j_cycle,j_chan,:) = DTOF_CentralMom(Time_ns, dtof);

            if sum(DTOF_SameND.Exp1_Deep_Red(j_cycle,1,:)) ~= 0
                [dtof,~] = DTOF_filter( squeeze(DTOF_SameND.Exp1_Deep_Red(j_cycle,j_chan,:)), Time_ns, bg_Red, 2, cut_ind_mom_red, 3);
                mom_temp = DTOF_CentralMom(Time_ns, dtof);

                correct_N_red = correct_N_red * (mom_temp(1) / Mom.Exp1_Deep_Red(j_cycle,j_chan,1));
            end
            Mom.Exp1_Deep_Red(j_cycle,j_chan,1) = Mom.Exp1_Deep_Red(j_cycle,j_chan,1) * correct_N_red;
        end
        
        % Step 3.  Calculate changes in moments relative to chosen baseline
        for j_cycle = 1:21
            delMom.Exp1_TwoL_Blu(j_cycle,j_chan,:) = DTOF_DelMom(Mom.Exp1_TwoL_Blu(c_base(1),j_chan,:), Mom.Exp1_TwoL_Blu(j_cycle,j_chan,:));
            delMom.Exp1_Deep_Red(j_cycle,j_chan,:) = DTOF_DelMom(Mom.Exp1_Deep_Red(c_base(1),j_chan,:), Mom.Exp1_Deep_Red(j_cycle,j_chan,:));
        end

        % Step 4.  Calculate photon noise for each moment, i.e. weight (although makes little difference to the results of LMA)
        temp1 = sqrt( 1  );
        temp1(2) = sqrt( Mom.Exp1_TwoL_Blu(c_base(1),j_chan,3)  );
        temp1(3) = sqrt( (Mom.Exp1_TwoL_Blu(c_base(1),j_chan,5) - Mom.Exp1_TwoL_Blu(c_base(1),j_chan,3)^2)  );
        Options_LMA.weight_noise = temp1;
        
        % Step 5.  Call LMA function to determine changes in optical properties
        for j_cycle = j_cycle_all
            [X_found, ~] = LMA_2_FittingDelMom(delOptProp_all, OptProp_base, n, rho, squeeze(delMom.Exp1_TwoL_Blu(j_cycle,j_chan,:)),...
                squeeze(Mom.Exp1_TwoL_Blu(j_cycle,j_chan,:)), Time_ns, cut_ind_mom, IRF_shifted.Exp1_Blu(:,j_chan), Options_LMA);

            X_found(1) = X_found(1) + OptProp_base(1);
            X_found(2) = X_found(2) + OptProp_base(2);
            LMA2_Mua1Mua2.Exp1_TwoL_Blu(j_cycle,j_chan,1:2) = X_found;
        end
        % ----------- For Exp. #1 (end) -----------
        
        

        % The code for Exp. 2 is the same as for Exp. 1 (copy-pasted the above)
        % ----------- For Exp. #2 (start) -----------
        % Step 1.  Set baseline optical properties of both layers
        OptProp_base(1) = LMA1_MuaMusp.Exp2_SupL_Red(c_base(2),j_chan,1); % Mua  of Superficial layer
        OptProp_base(4) = LMA1_MuaMusp.Exp2_SupL_Red(c_base(2),j_chan,2); % Musp of Superficial layer
        OptProp_base(2:3) = LMA1_MuaMusp.ConstLayer_Exp2_Deep_Blu(1,j_chan,1); % Mua  of Deep layer
        OptProp_base(5:6) = LMA1_MuaMusp.ConstLayer_Exp2_Deep_Blu(1,j_chan,2); % Musp of Deep layer

        [ ~, cut_ind_mom]     = DTOF_filter( squeeze(DTOF.Exp2_TwoL_Blu(c_base(2),j_chan,:)), Time_ns, bg_Blu, 2, cut_lim_mom, 1); % Region of time-channels for calculating moments
        [ ~, cut_ind_mom_red] = DTOF_filter( squeeze(DTOF.Exp2_SupL_Red(c_base(2),j_chan,:)), Time_ns, bg_Red, 2, cut_lim_mom, 1);
        
        % Step 2.  Calculate moments 
        correct_N = 1;
        correct_N_red = 1;
        for j_cycle = 1:21
            [dtof,~] = DTOF_filter( squeeze(DTOF.Exp2_TwoL_Blu(j_cycle,j_chan,:)), Time_ns, bg_Blu, 2, cut_ind_mom, 3);
            Mom.Exp2_TwoL_Blu(j_cycle,j_chan,:) = DTOF_CentralMom(Time_ns, dtof);
            
            % Correct Ntot if changed laser intensity - this doesnt affect any calculations since we dont use Ntot in any calculations
            if sum(DTOF_SameND.Exp2_TwoL_Blu(j_cycle,1,:)) ~= 0 
                [dtof, ~] = DTOF_filter( squeeze(DTOF_SameND.Exp2_TwoL_Blu(j_cycle,j_chan,:)), Time_ns, bg_Blu, 2, cut_ind_mom, 3);
                mom_temp = DTOF_CentralMom(Time_ns, dtof);
                
                correct_N = correct_N * (mom_temp(1) / Mom.Exp2_TwoL_Blu(j_cycle,j_chan,1));
            end
            Mom.Exp2_TwoL_Blu(j_cycle,j_chan,1) = Mom.Exp2_TwoL_Blu(j_cycle,j_chan,1) * correct_N;

            % Repeat for Red detector:
            [dtof,~] = DTOF_filter( squeeze(DTOF.Exp2_SupL_Red(j_cycle,j_chan,:)), Time_ns, bg_Red, 2, cut_ind_mom_red, 3);
            Mom.Exp2_SupL_Red(j_cycle,j_chan,:) = DTOF_CentralMom(Time_ns, dtof);

            if sum(DTOF_SameND.Exp2_SupL_Red(j_cycle,1,:)) ~= 0
                [dtof,~] = DTOF_filter( squeeze(DTOF_SameND.Exp2_SupL_Red(j_cycle,j_chan,:)), Time_ns, bg_Red, 2, cut_ind_mom_red, 3);
                mom_temp = DTOF_CentralMom(Time_ns, dtof);

                correct_N_red = correct_N_red * (mom_temp(1) / Mom.Exp2_SupL_Red(j_cycle,j_chan,1));
            end
            Mom.Exp2_SupL_Red(j_cycle,j_chan,1) = Mom.Exp2_SupL_Red(j_cycle,j_chan,1) * correct_N_red;
        end
        
        % Step 3.  Calculate changes in moments relative to chosen baseline
        for j_cycle = j_cycle_all
            delMom.Exp2_TwoL_Blu(j_cycle,j_chan,:) = DTOF_DelMom(Mom.Exp2_TwoL_Blu(c_base(1),j_chan,:), Mom.Exp2_TwoL_Blu(j_cycle,j_chan,:));
            delMom.Exp2_SupL_Red(j_cycle,j_chan,:) = DTOF_DelMom(Mom.Exp2_SupL_Red(c_base(1),j_chan,:), Mom.Exp2_SupL_Red(j_cycle,j_chan,:));
        end

        % Step 4.  Calculate photon noise for each moment, i.e. weight (although makes little difference to the results of LMA)
        temp1 = sqrt( 1  );
        temp1(2) = sqrt( Mom.Exp2_TwoL_Blu(c_base(2),j_chan,3)  );
        temp1(3) = sqrt( (Mom.Exp2_TwoL_Blu(c_base(2),j_chan,5) - Mom.Exp2_TwoL_Blu(c_base(2),j_chan,3)^2)  );
        Options_LMA.weight_noise = temp1;
        
        % Step 5.  Call LMA function to determine changes in optical properties
        for j_cycle = j_cycle_all
            [X_found, ~] = LMA_2_FittingDelMom(delOptProp_all, OptProp_base, n, rho, squeeze(delMom.Exp2_TwoL_Blu(j_cycle,j_chan,:)),...
                            squeeze(Mom.Exp2_TwoL_Blu(j_cycle,j_chan,:)), Time_ns, cut_ind_mom, IRF_shifted.Exp2_Blu(:,j_chan), Options_LMA);

            X_found(1) = X_found(1) + OptProp_base(1);
            X_found(2) = X_found(2) + OptProp_base(2);
            LMA2_Mua1Mua2.Exp2_TwoL_Blu(j_cycle,j_chan,1:2) = X_found;
        end
        % ----------- For Exp. #2 (end) -----------

        disp([char(datetime('now','Format','HH:mm:ss')) '  Finished channel ' num2str(j_chan)])
    end

    % Saving:
    if do_save == 1
        save([folder_calc '\' file_name], 'Mom', 'delMom', 'LMA2_Mua1Mua2', 'cut_lim_mom', 'c_base', 'L_known')
        
        disp([char(datetime('now','Format','HH:mm:ss')) '  SAVED ''LMA2_Mua1Mua2'' INTO VARIABLE ' file_name])
    end
end
if do_load == 0
    disp([char(datetime('now','Format','HH:mm:ss')) '  Step 3. Determined two-layered absorption'])
else
    disp([char(datetime('now','Format','HH:mm:ss')) '  Step 3. LOADED previously Determined two-layered absorption'])
end


% PLOT:
show_chan = 8; % 8th is the channel shown in Publication 2023
show_fig = {'5b','6a','6b'};
Ink_plotting(DTOF,Time_ns,IRF,delMom,LMA1_MuaMusp,LMA2_Mua1Mua2,LMA1_Optional, set_dim, show_chan, show_fig);


return
%% Step 4. (not in Publication)    Curve-fitting for two-layered model

% For saving/loading the calculated data:
file_name = 'Ink_LMA1_Optional_vr1'; % cut_lim_fit = [0.85 0.01]
% file_name = 'Ink_LMA1_Optional_vr2';

do_load = 1; % 1 == load previously calculated results. 0 == Do calculation
% do_load = 0; 

do_save = 1; % 1 == Save result in 'folder_calc'. Will use 'file_name' for the saved file
% do_save = 0;

if do_load == 1
%     load([folder_calc '\' file_name])
    load(file_name)

else
    L_known = 14.5;
    
    which_case = 1;

    % Choose which optical properties to determine:
    % Note, OptProp_base === [Mua1 Mua2 Mua3 Musp1 Musp2 Musp3 L1 L2]
    if which_case == 1
        OptProp_base = [-7 -7 -8 -7 -8 -8 L_known 40]; % Determine Mua1 Mua2
    elseif which_case == 2
        OptProp_base = [-7 -7 -8 -7 -8 -8 L_known 40]; % Determine Mua1 Mua2 and Musp
    end
    
    x_unknowns = sum(OptProp_base==-7);
    LMA1_Optional = [];
    LMA1_Optional.Exp1_TwoL_Blu = zeros(21,14,x_unknowns);
    LMA1_Optional.Exp2_TwoL_Blu = zeros(21,14,x_unknowns);
    
    cut_lim_fit = [0.85 0.01]; % Region for fitting
    
    j_chan_all = 1:14;
    j_chan_all = 8;
    
    j_cycle_all = 1:21; 
%     j_cycle_all = [1 11 21]; 
    
    for j_chan = j_chan_all
        if which_case == 1 % Exp. 1, on Two-Layered, Blu detector
            OptProp_base(4) = LMA1_MuaMusp.ConstLayer_Exp1_SupL_Blu(1,j_chan,2); % Musp of Superficial
            OptProp_base(5) = LMA1_MuaMusp.Exp1_Deep_Red(c_base(1),j_chan,2); % Musp of Deep
        end
        for j_cycle = j_cycle_all
            [ dtof, cut_ind_fit] = DTOF_filter( squeeze(DTOF.Exp1_TwoL_Blu(j_cycle,j_chan,:)), Time_ns, bg_Blu, 2, cut_lim_fit, 1);
    
            [X_found, ~] = LMA_1_FittingDTOF(OptProp_base, n, rho, dtof, Time_ns, cut_ind_fit, IRF_shifted.Exp1_Blu(:,j_chan));
            LMA1_Optional.Exp1_TwoL_Blu(j_cycle,j_chan,1:2) = X_found(1:2);
        end
    
        if which_case == 1 % Exp. 2, on Two-Layered, Blu detector
            OptProp_base(4) = LMA1_MuaMusp.Exp2_SupL_Red(c_base(2),j_chan,2); % Musp of Superficial
            OptProp_base(5) = LMA1_MuaMusp.ConstLayer_Exp2_Deep_Blu(1,j_chan,2); % Musp of Deep
        end
        for j_cycle = j_cycle_all
            [ dtof, cut_ind_fit] = DTOF_filter( squeeze(DTOF.Exp2_TwoL_Blu(j_cycle,j_chan,:)), Time_ns, bg_Blu, 2, cut_lim_fit, 1);
    
            [X_found, ~] = LMA_1_FittingDTOF(OptProp_base, n, rho, dtof, Time_ns, cut_ind_fit, IRF_shifted.Exp2_Blu(:,j_chan));
            LMA1_Optional.Exp2_TwoL_Blu(j_cycle,j_chan,1:2) = X_found(1:2);
        end
    
        disp([char(datetime('now','Format','HH:mm:ss')) '  Finished channel ' num2str(j_chan)])
    end
    
    % Saving:
    if do_save == 1
        save([folder_calc '\' file_name], 'LMA1_Optional')

        disp([char(datetime('now','Format','HH:mm:ss')) '  Saved ''LMA1_Optional'' into variable ' file_name])
    end
    
    disp([char(datetime('now','Format','HH:mm:ss')) '  Finished Curve-fitting for two-layered model'])
end
if do_load == 0
    disp([char(datetime('now','Format','HH:mm:ss')) '  Step 4 (Optional). Determined two-layered absorption using LMA 1 (Curve-fitting)'])
else
    disp([char(datetime('now','Format','HH:mm:ss')) '  Step 4 (Optional). LOADED previously Determined two-layered absorption using LMA 1 (Curve-fitting)'])
end


% PLOT:
show_chan = 8;
show_fig = {'Optional_3'};
Ink_plotting(DTOF,Time_ns,IRF,delMom,LMA1_MuaMusp,LMA2_Mua1Mua2,LMA1_Optional, set_dim, show_chan, show_fig);