% Written by Aleh Sudakou (https://www.researchgate.net/profile/Aleh-Sudakou)
% Nałęcz Institute of Biocybernetics and Biomedical Engineering, Polish Academy of Sciences 
% This is one of the codes shared on  https://github.com/asudakou/Analyzing_TD-NIRS
% Last updated: 20 May 2023


%% Description:
% This script does data analysis for three experiments with blood.
% For full details, please refer to the publication.
%
% "Blood_plotting.m" plots all figures presented in Publication 2023
% 
% Note, for each step in this script, the codes for Exp. 2 and Exp. 3 are just like for Exp. 1 (copy-pasted two times)


%% Step 1. Load data and define empty variables
% These are all of the variables that will contain results or any relevant data

clear; % clc;

% Folder where to save the results, since this can be a very lengthy script to run:
folder_calc = 'C:\Users\asudakou\Documents\01 Work mine\MatLab codes\Wrote by me\TR analyses\Two-layered analysis\Aleh LMA fitting\Open-access\DATA calculated';

which_Exp = [1 2 3]; % For which experiments to do all calculations

% For Step 1 - Load the measured DTOF for 3 experiments with blood:
            if ismember(1, which_Exp)
                load('data_Blood_Exp1_Pub2023')
                Exp1_DTOF;
                Exp1_IRF; Exp1_IRF_shifted = [];
                Time_ns;
                CollectionTime;

                % Set the time 0 to the time when yeast was added to the deep compartment:
                Exp1_Time_min = (1:size(Exp1_DTOF.TwoL,1)) * CollectionTime / 60;
                Exp1_Time_min = Exp1_Time_min - 16.65; % 16.65 was measured by a stopwatch
            end
            if ismember(2, which_Exp) % For Exp. #2, copy-pasted the above 
                load('data_Blood_Exp2_Pub2023')
                Exp2_DTOF;
                Exp2_IRF; Exp2_IRF_shifted = [];

                Exp2_Time_min = (1:size(Exp2_DTOF.TwoL,1)) * CollectionTime / 60;
                Exp2_Time_min = Exp2_Time_min - 24.8;
            end
            if ismember(3, which_Exp) % For Exp. #3, copy-pasted the above 
                load('data_Blood_Exp3_Pub2023')
                Exp3_DTOF;
                Exp3_IRF; Exp3_IRF_shifted = [];

                Exp3_Time_min = (1:size(Exp3_DTOF.TwoL,1)) * CollectionTime / 60;
                Exp3_Time_min = Exp3_Time_min - 20.85;
            end

            % Moments
            if ismember(1, which_Exp)
                Exp1_Mom = [];
                Exp1_Mom.TwoL = zeros(size(Exp1_DTOF.TwoL,1),14,5);
                Exp1_Mom.Deep = zeros(size(Exp1_DTOF.TwoL,1),14,5);
            end
            if ismember(2, which_Exp)
                Exp2_Mom = [];
                Exp2_Mom.TwoL = zeros(size(Exp2_DTOF.TwoL,1),14,5);
                Exp2_Mom.Deep = zeros(size(Exp2_DTOF.TwoL,1),14,5);
            end
            if ismember(3, which_Exp)
                Exp3_Mom = [];
                Exp3_Mom.TwoL = zeros(size(Exp3_DTOF.TwoL,1),14,5);
                Exp3_Mom.Deep = zeros(size(Exp3_DTOF.TwoL,1),14,5);
            end

% For Step 2 - Perform Homogeneous fitting (LMA 1) to obtain Mua and Musp
            if ismember(1, which_Exp)
                Exp1_LMA1_MuaMusp = [];
                Exp1_LMA1_MuaMusp.TwoL = zeros(size(Exp1_DTOF.TwoL,1),14,2);
                Exp1_LMA1_MuaMusp.Deep = zeros(size(Exp1_DTOF.TwoL,1),14,2);
            end
            if ismember(2, which_Exp)
                Exp2_LMA1_MuaMusp = [];
                Exp2_LMA1_MuaMusp.TwoL = zeros(size(Exp2_DTOF.TwoL,1),14,2);
                Exp2_LMA1_MuaMusp.Deep = zeros(size(Exp2_DTOF.TwoL,1),14,2);
            end
            if ismember(3, which_Exp)
                Exp3_LMA1_MuaMusp = [];
                Exp3_LMA1_MuaMusp.TwoL = zeros(size(Exp3_DTOF.TwoL,1),14,2);
                Exp3_LMA1_MuaMusp.Deep = zeros(size(Exp3_DTOF.TwoL,1),14,2);
            end

% For Step 3 - Perform 2-layered fitting (LMA 2) to determine delMua1 and delMua2
            if ismember(1, which_Exp)
                Exp1_LMA2_Mua1Mua2 = [];
                Exp1_LMA2_Mua1Mua2.TwoL = zeros(size(Exp1_DTOF.TwoL,1),14,2);
            end
            if ismember(2, which_Exp)
                Exp2_LMA2_Mua1Mua2 = [];
                Exp2_LMA2_Mua1Mua2.TwoL = zeros(size(Exp2_DTOF.TwoL,1),14,2);
            end
            if ismember(3, which_Exp)
                Exp3_LMA2_Mua1Mua2 = [];
                Exp3_LMA2_Mua1Mua2.TwoL = zeros(size(Exp3_DTOF.TwoL,1),14,2);
            end
            
% For Step 4 - Calculate concentrations of Oxy and Deoxy
            if ismember(1, which_Exp)
                Exp1_Conc = [];
                Exp1_Conc.TwoL_LMA1 = zeros(size(Exp1_DTOF.TwoL,1),2);
                Exp1_Conc.TwoL_LMA2 = zeros(size(Exp1_DTOF.TwoL,1),2,2);
                Exp1_Conc.Deep_LMA1 = zeros(size(Exp1_DTOF.TwoL,1),2);
            end
            if ismember(2, which_Exp)
                Exp2_Conc = [];
                Exp2_Conc.TwoL_LMA1 = zeros(size(Exp2_DTOF.TwoL,1),2);
                Exp2_Conc.TwoL_LMA2 = zeros(size(Exp2_DTOF.TwoL,1),2,2);
                Exp2_Conc.Deep_LMA1 = zeros(size(Exp2_DTOF.TwoL,1),2);
            end
            if ismember(3, which_Exp)
                Exp3_Conc = [];
                Exp3_Conc.TwoL_LMA1 = zeros(size(Exp3_DTOF.TwoL,1),2);
                Exp3_Conc.TwoL_LMA2 = zeros(size(Exp3_DTOF.TwoL,1),2,2);
                Exp3_Conc.Deep_LMA1 = zeros(size(Exp3_DTOF.TwoL,1),2);
            end

% For Step 5 - Calculate StO2
            if ismember(1, which_Exp)
                Exp1_StO2 = [];
                Exp1_StO2.TwoL_LMA1 = zeros(size(Exp1_DTOF.TwoL,1),1);
                Exp1_StO2.TwoL_LMA2 = zeros(size(Exp1_DTOF.TwoL,1),2);
                Exp1_StO2.Deep_LMA1 = zeros(size(Exp1_DTOF.TwoL,1),1);
            end
            if ismember(2, which_Exp)
                Exp2_StO2 = [];
                Exp2_StO2.TwoL_LMA1 = zeros(size(Exp2_DTOF.TwoL,1),1);
                Exp2_StO2.TwoL_LMA2 = zeros(size(Exp2_DTOF.TwoL,1),2);
                Exp2_StO2.Deep_LMA1 = zeros(size(Exp2_DTOF.TwoL,1),1);
            end
            if ismember(3, which_Exp)
                Exp3_StO2 = [];
                Exp3_StO2.TwoL_LMA1 = zeros(size(Exp3_DTOF.TwoL,1),1);
                Exp3_StO2.TwoL_LMA2 = zeros(size(Exp3_DTOF.TwoL,1),2);
                Exp3_StO2.Deep_LMA1 = zeros(size(Exp3_DTOF.TwoL,1),1);
            end

% Define all constants and filter IRFs:

% set_dim = 1; % Setting this to 1 will make the figure sizes the same as we used for Publication 2023
set_dim = 0;

wavelengths = 680:12.5:868; % 16 wavelengths, unit nm

rho = 30; % Units: mm. Source-detector distance was always set to 30 mm

offset = 0.06/physconst('Lightspeed')*1e9; % Units: ns. Time offset of IRF due the 6 cm distance between fibers

n = 1.33; % Refractive index. Assuming it is the same as water because we measured on liquid phantoms where main component was water

bg_TwoL = [0.5 1]; % Region for subtracting background noise
bg_Deep = [0.5 1];
bg_irf = [7.3 8.3];

% Filter IRFs:  and choosing to use the IRF that was measured at the start
if ismember(1, which_Exp)
    Exp1_IRF_shifted.TwoL = zeros(1024, 16);
    Exp1_IRF_shifted.Deep = zeros(size(Exp1_IRF_shifted.TwoL));
    for j_chan = 1:16
        % Remove background and cut at 0% on both sides of the maximum, although this filtering of IRF makes little difference to the results of data analysis
        [Exp1_IRF_shifted.TwoL(:,j_chan), ~ ] = DTOF_filter( squeeze(Exp1_IRF.TwoL_start(j_chan,:)), Time_ns, bg_irf, 2, [0 0], 1);
        [Exp1_IRF_shifted.Deep(:,j_chan), ~ ] = DTOF_filter( squeeze(Exp1_IRF.Deep_start(j_chan,:)), Time_ns, bg_irf, 2, [0 0], 1);

        % Shift by offset
        Exp1_IRF_shifted.TwoL(:,j_chan) = interp1(Time_ns - offset, Exp1_IRF_shifted.TwoL(:,j_chan), Time_ns);
        Exp1_IRF_shifted.Deep(:,j_chan) = interp1(Time_ns - offset, Exp1_IRF_shifted.Deep(:,j_chan), Time_ns);
    end
end
if ismember(2, which_Exp) % For Exp. #2, copy-pasted the above 
    Exp2_IRF_shifted.TwoL = zeros(1024, 16);
    Exp2_IRF_shifted.Deep = zeros(size(Exp2_IRF_shifted.TwoL));
    for j_chan = 1:16
        % Remove background and cut at 0% on both sides of the maximum, although this filtering of IRF makes little difference to the results of data analysis
        [Exp2_IRF_shifted.TwoL(:,j_chan), ~ ] = DTOF_filter( squeeze(Exp2_IRF.TwoL_start(j_chan,:)), Time_ns, bg_irf, 2, [0 0], 1);
        [Exp2_IRF_shifted.Deep(:,j_chan), ~ ] = DTOF_filter( squeeze(Exp2_IRF.Deep_start(j_chan,:)), Time_ns, bg_irf, 2, [0 0], 1);

        % Shift by offset
        Exp2_IRF_shifted.TwoL(:,j_chan) = interp1(Time_ns - offset, Exp2_IRF_shifted.TwoL(:,j_chan), Time_ns);
        Exp2_IRF_shifted.Deep(:,j_chan) = interp1(Time_ns - offset, Exp2_IRF_shifted.Deep(:,j_chan), Time_ns);
    end
end
if ismember(3, which_Exp) % For Exp. #3, copy-pasted the above 
    Exp3_IRF_shifted.TwoL = zeros(1024, 16);
    Exp3_IRF_shifted.Deep = zeros(size(Exp3_IRF_shifted.TwoL));
    for j_chan = 1:16
        % Remove background and cut at 0% on both sides of the maximum, although this filtering of IRF makes little difference to the results of data analysis
        [Exp3_IRF_shifted.TwoL(:,j_chan), ~ ] = DTOF_filter( squeeze(Exp3_IRF.TwoL_start(j_chan,:)), Time_ns, bg_irf, 2, [0 0], 1);
        [Exp3_IRF_shifted.Deep(:,j_chan), ~ ] = DTOF_filter( squeeze(Exp3_IRF.Deep_start(j_chan,:)), Time_ns, bg_irf, 2, [0 0], 1);

        % Shift by offset
        Exp3_IRF_shifted.TwoL(:,j_chan) = interp1(Time_ns - offset, Exp3_IRF_shifted.TwoL(:,j_chan), Time_ns);
        Exp3_IRF_shifted.Deep(:,j_chan) = interp1(Time_ns - offset, Exp3_IRF_shifted.Deep(:,j_chan), Time_ns);
    end
end


% Optional: How much moments of IRF change between start and finish:
do_this = 0; 
if do_this == 1
    j_chan = 7; % For which channel
    Mom_IRF = zeros(5,2);

    if ismember(1, which_Exp)
        [dtof, ~ ] = DTOF_filter( squeeze(Exp1_IRF.TwoL_start(j_chan,:)), Time_ns, bg_irf, 2, [0 0], 1);
        Mom_IRF(:,1) = DTOF_CentralMom(Time_ns, dtof);
        [dtof, ~ ] = DTOF_filter( squeeze(Exp1_IRF.TwoL_end(j_chan,:)), Time_ns, bg_irf, 2, [0 0], 1);
        Mom_IRF(:,2) = DTOF_CentralMom(Time_ns, dtof);
        delMom_IRF = DTOF_DelMom(Mom_IRF(:,1),Mom_IRF(:,2));
        disp(['For Exp.1 the moments of IRF chagned by : ' num2str(delMom_IRF')])

        [dtof, ~ ] = DTOF_filter( squeeze(Exp1_IRF.Deep_start(j_chan,:)), Time_ns, bg_irf, 2, [0 0], 1);
        Mom_IRF(:,1) = DTOF_CentralMom(Time_ns, dtof);
        [dtof, ~ ] = DTOF_filter( squeeze(Exp1_IRF.Deep_end(j_chan,:)), Time_ns, bg_irf, 2, [0 0], 1);
        Mom_IRF(:,2) = DTOF_CentralMom(Time_ns, dtof);
        delMom_IRF = DTOF_DelMom(Mom_IRF(:,1),Mom_IRF(:,2));
        disp(['And for the other laser-detector pair by : ' num2str(delMom_IRF')])
    end
    
    if ismember(2, which_Exp) % For Exp. #2, copy-pasted the above 
        [dtof, ~ ] = DTOF_filter( squeeze(Exp2_IRF.TwoL_start(j_chan,:)), Time_ns, bg_irf, 2, [0 0], 1);
        Mom_IRF(:,1) = DTOF_CentralMom(Time_ns, dtof);
        [dtof, ~ ] = DTOF_filter( squeeze(Exp2_IRF.TwoL_end(j_chan,:)), Time_ns, bg_irf, 2, [0 0], 1);
        Mom_IRF(:,2) = DTOF_CentralMom(Time_ns, dtof);
        delMom_IRF = DTOF_DelMom(Mom_IRF(:,1),Mom_IRF(:,2));
        disp(['For Exp.2 the moments of IRF chagned by : ' num2str(delMom_IRF')])

        [dtof, ~ ] = DTOF_filter( squeeze(Exp2_IRF.Deep_start(j_chan,:)), Time_ns, bg_irf, 2, [0 0], 1);
        Mom_IRF(:,1) = DTOF_CentralMom(Time_ns, dtof);
        [dtof, ~ ] = DTOF_filter( squeeze(Exp2_IRF.Deep_end(j_chan,:)), Time_ns, bg_irf, 2, [0 0], 1);
        Mom_IRF(:,2) = DTOF_CentralMom(Time_ns, dtof);
        delMom_IRF = DTOF_DelMom(Mom_IRF(:,1),Mom_IRF(:,2));
        disp(['And for the other laser-detector pair by : ' num2str(delMom_IRF')])
    end
    
    if ismember(3, which_Exp) % For Exp. #3, copy-pasted the above 
        [dtof, ~ ] = DTOF_filter( squeeze(Exp3_IRF.TwoL_start(j_chan,:)), Time_ns, bg_irf, 2, [0 0], 1);
        Mom_IRF(:,1) = DTOF_CentralMom(Time_ns, dtof);
        [dtof, ~ ] = DTOF_filter( squeeze(Exp3_IRF.TwoL_end(j_chan,:)), Time_ns, bg_irf, 2, [0 0], 1);
        Mom_IRF(:,2) = DTOF_CentralMom(Time_ns, dtof);
        delMom_IRF = DTOF_DelMom(Mom_IRF(:,1),Mom_IRF(:,2));
        disp(['For Exp.3 the moments of IRF chagned by : ' num2str(delMom_IRF')])

        [dtof, ~ ] = DTOF_filter( squeeze(Exp3_IRF.Deep_start(j_chan,:)), Time_ns, bg_irf, 2, [0 0], 1);
        Mom_IRF(:,1) = DTOF_CentralMom(Time_ns, dtof);
        [dtof, ~ ] = DTOF_filter( squeeze(Exp3_IRF.Deep_end(j_chan,:)), Time_ns, bg_irf, 2, [0 0], 1);
        Mom_IRF(:,2) = DTOF_CentralMom(Time_ns, dtof);
        delMom_IRF = DTOF_DelMom(Mom_IRF(:,1),Mom_IRF(:,2));
        disp(['And for the other laser-detector pair by : ' num2str(delMom_IRF')])
    end
end


%% Select region for baselines:
Base_minutes =   [14.6 15.6 31;  58.8 59.8 99]; % In minutes (start, end, until what time to use this baseline)
Base_minutes_3 = [14.6 15.6 30;  58.8 59.8 99]; 

% Base_minutes =   [-2.2 -0.2 31;  43.8 44.8 99]; % Another possible region for baselines
% Base_minutes_2 = [-2.2 -0.2 30;  43.8 44.8 99]; 

Base_index = zeros(size(Base_minutes,1),3,3);

for j_1 = 1:size(Base_minutes,1)
    for j_2 = 1:size(Base_minutes,2)
        [~, Base_index(j_1,j_2,1)] = min(abs(Base_minutes(j_1,j_2) - Exp1_Time_min)); % Exp 1

        [~, Base_index(j_1,j_2,2)] = min(abs(Base_minutes(j_1,j_2) - Exp2_Time_min)); % Exp 2

        [~, Base_index(j_1,j_2,3)] = min(abs(Base_minutes_3(j_1,j_2) - Exp3_Time_min)); % Exp 3
    end
end


%% Check how DTOFs look like during two baselines (max and min StO2, which corresponds to max and min MUA)
% This allows checking what region is chosen for calculating background noise, which is later subtracted

do_this = 0; 
if do_this == 1
    cut_bg = 0; % 2 = subtract,   0 = not subtract
    j_chan_all = 3:11;
    
    fig = 1:3;
    for j = fig; figure(j); clf; end
    for j_chan = j_chan_all
        temp = Base_index(1,1,1) : Base_index(1,2,1); % Exp.1  baseline 1
        [ ~, ~] = DTOF_filter( squeeze(mean(Exp1_DTOF.TwoL(temp,j_chan,:),1)), Time_ns, bg_TwoL, cut_bg, 0, 0, {fig(1), 'red'});
        temp = Base_index(2,1,1) : Base_index(2,2,1); % Exp.1  baseline 2
        [ ~, ~] = DTOF_filter( squeeze(mean(Exp1_DTOF.TwoL(temp,j_chan,:),1)), Time_ns, bg_TwoL, cut_bg, 0, 0, {fig(1), 'blue'});
    end
    plot([1 1]*bg_TwoL(1), ylim*0.8,'Color','red','LineWidth',2)
    plot([1 1]*bg_TwoL(2), ylim*0.8,'Color','red','LineWidth',2); title('Exp. 1   Two-layered medium')
    for j_chan = j_chan_all
        temp = Base_index(1,1,2) : Base_index(1,2,2); % Exp.2  baseline 1
        [ ~, ~] = DTOF_filter( squeeze(mean(Exp2_DTOF.TwoL(temp,j_chan,:),1)), Time_ns, bg_TwoL, cut_bg, 0, 0, {fig(2),'norm'});
        temp = Base_index(2,1,2) : Base_index(2,2,2); % Exp.2  baseline 2
        [ ~, ~] = DTOF_filter( squeeze(mean(Exp2_DTOF.TwoL(temp,j_chan,:),1)), Time_ns, bg_TwoL, cut_bg, 0, 0, {fig(2),'norm'});
    end
    plot([1 1]*bg_TwoL(1), ylim*0.8,'Color','red','LineWidth',2)
    plot([1 1]*bg_TwoL(2), ylim*0.8,'Color','red','LineWidth',2); title('Exp. 2   Two-layered medium')
    for j_chan = j_chan_all
        temp = Base_index(1,1,3) : Base_index(1,2,3); % Exp.3  baseline 1
        [ ~, ~] = DTOF_filter( squeeze(mean(Exp3_DTOF.TwoL(temp,j_chan,:),1)), Time_ns, bg_TwoL, cut_bg, 0, 0, {fig(3)});
        temp = Base_index(2,1,3) : Base_index(2,2,3); % Exp.3  baseline 2
        [ ~, ~] = DTOF_filter( squeeze(mean(Exp3_DTOF.TwoL(temp,j_chan,:),1)), Time_ns, bg_TwoL, cut_bg, 0, 0, {fig(3)});
    end
    plot([1 1]*bg_TwoL(1), ylim*0.8,'Color','red','LineWidth',2)
    plot([1 1]*bg_TwoL(2), ylim*0.8,'Color','red','LineWidth',2); title('Exp. 3   Two-layered medium')
    

    % Repeat for the second source-detector pair (on deep layer):
    fig = 4:6;
    for j = fig; figure(j); clf; end
    for j_chan = j_chan_all
        temp = Base_index(1,1,1) : Base_index(1,2,1); % Exp.1  baseline 1
        [ ~, ~] = DTOF_filter( squeeze(mean(Exp1_DTOF.Deep(temp,j_chan,:),1)), Time_ns, bg_Deep, cut_bg, 0, 0, {fig(1), 'red'});
        temp = Base_index(2,1,1) : Base_index(2,2,1); % Exp.1  baseline 2
        [ ~, ~] = DTOF_filter( squeeze(mean(Exp1_DTOF.Deep(temp,j_chan,:),1)), Time_ns, bg_Deep, cut_bg, 0, 0, {fig(1), 'blue'});
    end
    plot([1 1]*bg_Deep(1), ylim*0.8,'Color','red','LineWidth',2)
    plot([1 1]*bg_Deep(2), ylim*0.8,'Color','red','LineWidth',2); title('Exp. 3   Deep layer')
    for j_chan = j_chan_all
        temp = Base_index(1,1,2) : Base_index(1,2,2); % Exp.2  baseline 1
        [ ~, ~] = DTOF_filter( squeeze(mean(Exp2_DTOF.Deep(temp,j_chan,:),1)), Time_ns, bg_Deep, cut_bg, 0, 0, {fig(2),'norm'});
        temp = Base_index(2,1,2) : Base_index(2,2,2); % Exp.2  baseline 2
        [ ~, ~] = DTOF_filter( squeeze(mean(Exp2_DTOF.Deep(temp,j_chan,:),1)), Time_ns, bg_Deep, cut_bg, 0, 0, {fig(2),'norm'});
    end
    plot([1 1]*bg_Deep(1), ylim*0.8,'Color','red','LineWidth',2)
    plot([1 1]*bg_Deep(2), ylim*0.8,'Color','red','LineWidth',2); title('Exp. 3   Deep layer')
    for j_chan = j_chan_all
        temp = Base_index(1,1,3) : Base_index(1,2,3); % Exp.3  baseline 1
        [ ~, ~] = DTOF_filter( squeeze(mean(Exp3_DTOF.Deep(temp,j_chan,:),1)), Time_ns, bg_Deep, cut_bg, 0, 0, {fig(3)});
        temp = Base_index(2,1,3) : Base_index(2,2,3); % Exp.3  baseline 2
        [ ~, ~] = DTOF_filter( squeeze(mean(Exp3_DTOF.Deep(temp,j_chan,:),1)), Time_ns, bg_Deep, cut_bg, 0, 0, {fig(3)});
    end
    plot([1 1]*bg_Deep(1), ylim*0.8,'Color','red','LineWidth',2)
    plot([1 1]*bg_Deep(2), ylim*0.8,'Color','red','LineWidth',2); title('Exp. 3   Deep layer')
end


%% Calculate moments and plot

% For saving/loading the calculated data:
file_name_1 = 'Blood_Exp1_Mom_vr1';
file_name_2 = 'Blood_Exp2_Mom_vr1';
file_name_3 = 'Blood_Exp3_Mom_vr1';
% file_name_1 = 'Blood_Exp1_Mom_vr2';
% file_name_2 = 'Blood_Exp2_Mom_vr2';
% file_name_3 = 'Blood_Exp3_Mom_vr2';

% do_load = 1; % 1 == load previously calculated results. 0 == Do calculation
do_load = 0;

% do_save = 1; % 1 == Save result in 'folder_calc'. Will use 'file_name' for the saved file
do_save = 0;

j_chan_all = 1:14; 
% j_chan_all = 3;

cut_lim_mom = [0.25 0.03]; % Used for Publication 2023
% cut_lim_mom = [0.25 0.05];
% cut_lim_mom = [0.001 0.001];

to_skip = 1;

for j_Exp = which_Exp % For which experiment to perform the calculation
% for j_Exp = 1

    if isequal(1, j_Exp) % FIRST EXPERIMENT
        if do_load == 1
            load(file_name_1)
        else
            j_cycle_all = 1:to_skip:size(Exp1_DTOF.TwoL,1);
    
            for j_chan = j_chan_all
                temp = Base_index(1,1,1) : Base_index(1,2,1);
                [ ~, cut_ind_mom] = DTOF_filter( squeeze(mean(Exp1_DTOF.TwoL(temp,j_chan,:),1)), Time_ns, bg_TwoL, 2, cut_lim_mom, 1);
                [ ~, cut_ind_mom_2] = DTOF_filter( squeeze(mean(Exp1_DTOF.Deep(temp,j_chan,:),1)), Time_ns, bg_Deep, 2, cut_lim_mom, 1);
    
                for j_cycle = j_cycle_all
                    % First source-detector pair:
                    [ dtof, ~] = DTOF_filter( squeeze(Exp1_DTOF.TwoL(j_cycle,j_chan,:)), Time_ns, bg_TwoL, 2, cut_ind_mom, 0);
                    Exp1_Mom.TwoL(j_cycle,j_chan,:) = DTOF_CentralMom(Time_ns(cut_ind_mom(1):cut_ind_mom(2)), dtof(cut_ind_mom(1):cut_ind_mom(2)));
    
                    % Second source-detector pair:
                    [ dtof, ~] = DTOF_filter( squeeze(Exp1_DTOF.Deep(j_cycle,j_chan,:)), Time_ns, bg_Deep, 2, cut_ind_mom_2, 0);
                    Exp1_Mom.Deep(j_cycle,j_chan,:) = DTOF_CentralMom(Time_ns(cut_ind_mom_2(1):cut_ind_mom_2(2)), dtof(cut_ind_mom_2(1):cut_ind_mom_2(2)));
                end
            end
            disp([char(datetime('now','Format','HH:mm:ss')) '  Exp. ' num2str(j_Exp) '.  Finished calculation moments'])
            if do_save == 1
                save([folder_calc '\' file_name_1], 'Exp1_Mom')
                disp([char(datetime('now','Format','HH:mm:ss')) '  Saved ''Exp1_Mom'' into variable ''' file_name_1 ''''])
            end
        end

        
    elseif isequal(2, j_Exp) % SECOND EXPERIMENT (copy-pasted the above)
        if do_load == 1
            load(file_name_2)
        else
            j_cycle_all = 1:to_skip:size(Exp2_DTOF.TwoL,1);

            for j_chan = j_chan_all
                temp = Base_index(1,1,1) : Base_index(1,2,1); % Use the same indeces as for Exp. 1
                [ ~, cut_ind_mom] = DTOF_filter( squeeze(mean(Exp1_DTOF.TwoL(temp,j_chan,:),1)), Time_ns, bg_TwoL, 2, cut_lim_mom, 1);
                [ ~, cut_ind_mom_2] = DTOF_filter( squeeze(mean(Exp1_DTOF.Deep(temp,j_chan,:),1)), Time_ns, bg_Deep, 2, cut_lim_mom, 1);
%                 temp = Base_index(1,1,2) : Base_index(1,2,2);
%                 [ ~, cut_ind_mom] = DTOF_filter( squeeze(mean(Exp2_DTOF.TwoL(temp,j_chan,:),1)), Time_ns, bg_TwoL, 2, cut_lim_mom, 1);

                for j_cycle = j_cycle_all
                    % First source-detector pair:
                    [ dtof, ~] = DTOF_filter( squeeze(Exp2_DTOF.TwoL(j_cycle,j_chan,:)), Time_ns, bg_TwoL, 2, cut_ind_mom, 0);
                    Exp2_Mom.TwoL(j_cycle,j_chan,:) = DTOF_CentralMom(Time_ns(cut_ind_mom(1):cut_ind_mom(2)), dtof(cut_ind_mom(1):cut_ind_mom(2)));

                    % Second source-detector pair:
                    [ dtof, ~] = DTOF_filter( squeeze(Exp2_DTOF.Deep(j_cycle,j_chan,:)), Time_ns, bg_Deep, 2, cut_ind_mom_2, 0);
                    Exp2_Mom.Deep(j_cycle,j_chan,:) = DTOF_CentralMom(Time_ns(cut_ind_mom_2(1):cut_ind_mom_2(2)), dtof(cut_ind_mom_2(1):cut_ind_mom_2(2)));
                end
            end
            disp([char(datetime('now','Format','HH:mm:ss')) '  Exp. ' num2str(j_Exp) '.  Finished calculation moments'])
            if do_save == 1
                save([folder_calc '\' file_name_2], 'Exp2_Mom')
                disp([char(datetime('now','Format','HH:mm:ss')) '  Saved ''Exp2_Mom'' into variable ''' file_name_2 ''''])
            end
        end

    elseif isequal(3, j_Exp) % THIRD EXPERIMENT (copy-pasted the above)
        if do_load == 1
            load(file_name_3)
        else
            j_cycle_all = 1:to_skip:size(Exp3_DTOF.TwoL,1);

            for j_chan = j_chan_all
                temp = Base_index(1,1,1) : Base_index(1,2,1); % Use the same indeces as for Exp. 1
                [ ~, cut_ind_mom] = DTOF_filter( squeeze(mean(Exp1_DTOF.TwoL(temp,j_chan,:),1)), Time_ns, bg_TwoL, 2, cut_lim_mom, 1);
                [ ~, cut_ind_mom_2] = DTOF_filter( squeeze(mean(Exp1_DTOF.Deep(temp,j_chan,:),1)), Time_ns, bg_Deep, 2, cut_lim_mom, 1);
%                 temp = Base_index(1,1,3) : Base_index(1,2,3);
%                 [ ~, cut_ind_mom] = DTOF_filter( squeeze(mean(Exp3_DTOF.TwoL(temp,j_chan,:),1)), Time_ns, bg_TwoL, 2, cut_lim_mom, 1);

                for j_cycle = j_cycle_all
                    % First source-detector pair:
                    [ dtof, ~] = DTOF_filter( squeeze(Exp3_DTOF.TwoL(j_cycle,j_chan,:)), Time_ns, bg_TwoL, 2, cut_ind_mom, 0);
                    Exp3_Mom.TwoL(j_cycle,j_chan,:) = DTOF_CentralMom(Time_ns(cut_ind_mom(1):cut_ind_mom(2)), dtof(cut_ind_mom(1):cut_ind_mom(2)));

                    % Second source-detector pair:
                    [ dtof, ~] = DTOF_filter( squeeze(Exp3_DTOF.Deep(j_cycle,j_chan,:)), Time_ns, bg_Deep, 2, cut_ind_mom_2, 0);
                    Exp3_Mom.Deep(j_cycle,j_chan,:) = DTOF_CentralMom(Time_ns(cut_ind_mom_2(1):cut_ind_mom_2(2)), dtof(cut_ind_mom_2(1):cut_ind_mom_2(2)));
                end
            end
            disp([char(datetime('now','Format','HH:mm:ss')) '  Exp. ' num2str(j_Exp) '.  Finished calculation moments'])
            if do_save == 1
                save([folder_calc '\' file_name_3], 'Exp3_Mom')
                disp([char(datetime('now','Format','HH:mm:ss')) '  Saved ''Exp3_Mom'' into variable ''' file_name_3 ''''])
            end
        end
    end
end


show_fig = {'7a', '7b', '7c'};
Blood_plotting(Exp1_Time_min,Exp2_Time_min,Exp3_Time_min,    Exp1_DTOF,Exp2_DTOF,Exp3_DTOF,   Time_ns,   Exp1_IRF,Exp2_IRF,Exp3_IRF,   Base_index,...
            Exp1_Mom,Exp2_Mom,Exp3_Mom,  Exp1_LMA1_MuaMusp,Exp2_LMA1_MuaMusp,Exp3_LMA1_MuaMusp,...
            Exp1_LMA2_Mua1Mua2,Exp2_LMA2_Mua1Mua2,Exp3_LMA2_Mua1Mua2,...
            Exp1_Conc,Exp2_Conc,Exp3_Conc,  Exp1_StO2,Exp2_StO2,Exp3_StO2, set_dim, show_fig )


%% Step 2.  Perform Homogeneous fitting for all DTOF

% For saving/loading the calculated data:
file_name_1 = 'Blood_Exp1_LMA1_MuaMusp_vr1'; % Used: cut_lim_fit = [0.85 0.01]
file_name_2 = 'Blood_Exp2_LMA1_MuaMusp_vr1';
file_name_3 = 'Blood_Exp3_LMA1_MuaMusp_vr1';
% file_name_1 = 'Blood_Exp1_LMA1_MuaMusp_vr2';
% file_name_2 = 'Blood_Exp2_LMA1_MuaMusp_vr2';
% file_name_3 = 'Blood_Exp3_LMA1_MuaMusp_vr2';

do_load = 1; % 1 == load previously calculated results. 0 == Do calculation
% do_load = 0;

do_save = 1; % 1 == Save result in 'folder_calc'. Will use 'file_name' for the saved file
% do_save = 0;

j_chan_all = 1:14;
% j_chan_all = 3;

do_use_parpool = 1;
% do_use_parpool = 0;

if do_load == 0 && do_use_parpool == 1
    parpool
end

cut_lim_fit = [0.85 0.01];

OptProp_all = [-7 -8 -8 -7 -8 -8 40 40]; % [Mua1 Mua2 Mua3 Musp1 Musp2 Musp3 L1 L2]

to_skip = 1;
% to_skip = 500; % To calculate for not every data point


for j_Exp = which_Exp % For which experiment to perform the calculation
% for j_Exp = 1

    if isequal(1, j_Exp) % FIRST EXPERIMENT
        if do_load == 1
            load(file_name_1)
        else
            j_cycle_all = 1:to_skip:size(Exp1_DTOF.TwoL,1);
    %         j_cycle_all = [Base_index(1,1,1) : Base_index(1,2,1), Base_index(2,1,1) : Base_index(2,2,1)];
%             j_cycle_all = 2161; % For a specific data point
    
            for j_chan = j_chan_all
    
                if do_use_parpool == 0
                    for j_cycle = j_cycle_all
                        % First source-detector pair:
                        [ dtof, cut_ind_fit] = DTOF_filter( squeeze(Exp1_DTOF.TwoL(j_cycle,j_chan,:)), Time_ns, bg_TwoL, 2, cut_lim_fit, 1);
                        [X_found, ~] = LMA_1_FittingDTOF(OptProp_all, n, rho, dtof, Time_ns, cut_ind_fit, Exp1_IRF_shifted.TwoL(:,j_chan));
                        Exp1_LMA1_MuaMusp.TwoL(j_cycle,j_chan,:) = X_found(1:2);
    
                        % Second source-detector pair:
                        [ dtof, cut_ind_fit] = DTOF_filter( squeeze(Exp1_DTOF.Deep(j_cycle,j_chan,:)), Time_ns, bg_Deep, 2, cut_lim_fit, 1);
                        [X_found, ~] = LMA_1_FittingDTOF(OptProp_all, n, rho, dtof, Time_ns, cut_ind_fit, Exp1_IRF_shifted.Deep(:,j_chan));
                        Exp1_LMA1_MuaMusp.Deep(j_cycle,j_chan,:) = X_found(1:2);
                    end
                else
                    % The EXACT SAME CODE AS THE ABOVE, BUT WITH THE USE OF PARPOOL:
                    temp_res = zeros(length(j_cycle_all),2);
                    temp_res_2 = zeros(length(j_cycle_all),2);
                    temp_dtof = squeeze(Exp1_DTOF.TwoL(j_cycle_all,j_chan,:));
                    temp_dtof_2 = squeeze(Exp1_DTOF.Deep(j_cycle_all,j_chan,:));
                    temp_irf = Exp1_IRF_shifted.TwoL(:,j_chan);
                    temp_irf_2 = Exp1_IRF_shifted.Deep(:,j_chan);
                    parfor j = 1:length(j_cycle_all)
                        % First source-detector pair:
                        [ dtof, cut_ind_fit] = DTOF_filter( temp_dtof(j,:), Time_ns, bg_TwoL, 2, cut_lim_fit, 1);
                        [X_found, ~] = LMA_1_FittingDTOF(OptProp_all, n, rho, dtof, Time_ns, cut_ind_fit, temp_irf);
                        temp_res(j,:) = X_found(1:2);
    
                        % Second source-detector pair:
                        [ dtof, cut_ind_fit] = DTOF_filter( temp_dtof_2(j,:), Time_ns, bg_Deep, 2, cut_lim_fit, 1);
                        [X_found, ~] = LMA_1_FittingDTOF(OptProp_all, n, rho, dtof, Time_ns, cut_ind_fit, temp_irf_2);
                        temp_res_2(j,:) = X_found(1:2);
                    end
                    Exp1_LMA1_MuaMusp.TwoL(j_cycle_all,j_chan,1) = temp_res(:,1);
                    Exp1_LMA1_MuaMusp.TwoL(j_cycle_all,j_chan,2) = temp_res(:,2);
                    Exp1_LMA1_MuaMusp.Deep(j_cycle_all,j_chan,1) = temp_res_2(:,1);
                    Exp1_LMA1_MuaMusp.Deep(j_cycle_all,j_chan,2) = temp_res_2(:,2);
                end

                disp([char(datetime('now','Format','HH:mm:ss')) '  Exp. ' num2str(j_Exp) '. Finished LMA1 for channel ' num2str(j_chan)])
            end
            if do_save == 1
                save([folder_calc '\' file_name_1], 'Exp1_LMA1_MuaMusp')
                disp([char(datetime('now','Format','HH:mm:ss')) '  Saved ''Exp1_LMA1_MuaMusp'' into variable ''' file_name_1 ''''])
            end
        end
    

    elseif isequal(2, j_Exp) % SECOND EXPERIMENT (copy-pasted the above)
        if do_load == 1
            load(file_name_2)
        else
            j_cycle_all = 1:to_skip:size(Exp2_DTOF.TwoL,1);
    
            for j_chan = j_chan_all
    
                if do_use_parpool == 0
                    for j_cycle = j_cycle_all
                        % First source-detector pair:
                        [ dtof, cut_ind_fit] = DTOF_filter( squeeze(Exp2_DTOF.TwoL(j_cycle,j_chan,:)), Time_ns, bg_TwoL, 2, cut_lim_fit, 1);
                        [X_found, ~] = LMA_1_FittingDTOF(OptProp_all, n, rho, dtof, Time_ns, cut_ind_fit, Exp2_IRF_shifted.TwoL(:,j_chan));
                        Exp2_LMA1_MuaMusp.TwoL(j_cycle,j_chan,:) = X_found(1:2);
    
                        % Second source-detector pair:
                        [ dtof, cut_ind_fit] = DTOF_filter( squeeze(Exp2_DTOF.Deep(j_cycle,j_chan,:)), Time_ns, bg_Deep, 2, cut_lim_fit, 1);
                        [X_found, ~] = LMA_1_FittingDTOF(OptProp_all, n, rho, dtof, Time_ns, cut_ind_fit, Exp2_IRF_shifted.Deep(:,j_chan));
                        Exp2_LMA1_MuaMusp.Deep(j_cycle,j_chan,:) = X_found(1:2);
                    end
                else
                    % The EXACT SAME CODE AS THE ABOVE, BUT WITH THE USE OF PARPOOL:
                    temp_res = zeros(length(j_cycle_all),2);
                    temp_res_2 = zeros(length(j_cycle_all),2);
                    temp_dtof = squeeze(Exp2_DTOF.TwoL(j_cycle_all,j_chan,:));
                    temp_dtof_2 = squeeze(Exp2_DTOF.Deep(j_cycle_all,j_chan,:));
                    temp_irf = Exp2_IRF_shifted.TwoL(:,j_chan);
                    temp_irf_2 = Exp2_IRF_shifted.Deep(:,j_chan);
                    parfor j = 1:length(j_cycle_all)
                        % First source-detector pair:
                        [ dtof, cut_ind_fit] = DTOF_filter( temp_dtof(j,:), Time_ns, bg_TwoL, 2, cut_lim_fit, 1);
                        [X_found, ~] = LMA_1_FittingDTOF(OptProp_all, n, rho, dtof, Time_ns, cut_ind_fit, temp_irf);
                        temp_res(j,:) = X_found(1:2);
    
                        % Second source-detector pair:
                        [ dtof, cut_ind_fit] = DTOF_filter( temp_dtof_2(j,:), Time_ns, bg_Deep, 2, cut_lim_fit, 1);
                        [X_found, ~] = LMA_1_FittingDTOF(OptProp_all, n, rho, dtof, Time_ns, cut_ind_fit, temp_irf_2);
                        temp_res_2(j,:) = X_found(1:2);
                    end
                    Exp2_LMA1_MuaMusp.TwoL(j_cycle_all,j_chan,1) = temp_res(:,1);
                    Exp2_LMA1_MuaMusp.TwoL(j_cycle_all,j_chan,2) = temp_res(:,2);
                    Exp2_LMA1_MuaMusp.Deep(j_cycle_all,j_chan,1) = temp_res_2(:,1);
                    Exp2_LMA1_MuaMusp.Deep(j_cycle_all,j_chan,2) = temp_res_2(:,2);
                end
    
                disp([char(datetime('now','Format','HH:mm:ss')) '  Exp. ' num2str(j_Exp) '. Finished LMA1 for channel ' num2str(j_chan)])
            end
            if do_save == 1
                save([folder_calc '\' file_name_2], 'Exp2_LMA1_MuaMusp')
                disp([char(datetime('now','Format','HH:mm:ss')) '  Saved ''Exp2_LMA1_MuaMusp'' into variable ''' file_name_2 ''''])
            end
        end
    
    elseif isequal(3, j_Exp) % THIRD EXPERIMENT (copy-pasted the above)
        if do_load == 1
            load(file_name_3)
        else
            j_cycle_all = 1:to_skip:size(Exp3_DTOF.TwoL,1);
    %         j_cycle_all = [Base_index(1,1,3) : Base_index(1,2,3), Base_index(2,1,3) : Base_index(2,2,3)];
    
            for j_chan = j_chan_all
    
                if do_use_parpool == 0
                    for j_cycle = j_cycle_all
                        % First source-detector pair:
                        [ dtof, cut_ind_fit] = DTOF_filter( squeeze(Exp3_DTOF.TwoL(j_cycle,j_chan,:)), Time_ns, bg_TwoL, 2, cut_lim_fit, 1);
                        [X_found, ~] = LMA_1_FittingDTOF(OptProp_all, n, rho, dtof, Time_ns, cut_ind_fit, Exp3_IRF_shifted.TwoL(:,j_chan));
                        Exp3_LMA1_MuaMusp.TwoL(j_cycle,j_chan,:) = X_found(1:2);
    
                        % Second source-detector pair:
                        [ dtof, cut_ind_fit] = DTOF_filter( squeeze(Exp3_DTOF.Deep(j_cycle,j_chan,:)), Time_ns, bg_Deep, 2, cut_lim_fit, 1);
                        [X_found, ~] = LMA_1_FittingDTOF(OptProp_all, n, rho, dtof, Time_ns, cut_ind_fit, Exp3_IRF_shifted.Deep(:,j_chan));
                        Exp3_LMA1_MuaMusp.Deep(j_cycle,j_chan,:) = X_found(1:2);
                    end
                else
                    % The EXACT SAME CODE AS THE ABOVE, BUT WITH THE USE OF PARPOOL:
                    temp_res = zeros(length(j_cycle_all),2);
                    temp_res_2 = zeros(length(j_cycle_all),2);
                    temp_dtof = squeeze(Exp3_DTOF.TwoL(j_cycle_all,j_chan,:));
                    temp_dtof_2 = squeeze(Exp3_DTOF.Deep(j_cycle_all,j_chan,:));
                    temp_irf = Exp3_IRF_shifted.TwoL(:,j_chan);
                    temp_irf_2 = Exp3_IRF_shifted.Deep(:,j_chan);
                    parfor j = 1:length(j_cycle_all)
                        % First source-detector pair:
                        [ dtof, cut_ind_fit] = DTOF_filter( temp_dtof(j,:), Time_ns, bg_TwoL, 2, cut_lim_fit, 1);
                        [X_found, ~] = LMA_1_FittingDTOF(OptProp_all, n, rho, dtof, Time_ns, cut_ind_fit, temp_irf);
                        temp_res(j,:) = X_found(1:2);
    
                        % Second source-detector pair:
                        [ dtof, cut_ind_fit] = DTOF_filter( temp_dtof_2(j,:), Time_ns, bg_Deep, 2, cut_lim_fit, 1);
                        [X_found, ~] = LMA_1_FittingDTOF(OptProp_all, n, rho, dtof, Time_ns, cut_ind_fit, temp_irf_2);
                        temp_res_2(j,:) = X_found(1:2);
                    end
                    Exp3_LMA1_MuaMusp.TwoL(j_cycle_all,j_chan,1) = temp_res(:,1);
                    Exp3_LMA1_MuaMusp.TwoL(j_cycle_all,j_chan,2) = temp_res(:,2);
                    Exp3_LMA1_MuaMusp.Deep(j_cycle_all,j_chan,1) = temp_res_2(:,1);
                    Exp3_LMA1_MuaMusp.Deep(j_cycle_all,j_chan,2) = temp_res_2(:,2);
                end
    
                disp([char(datetime('now','Format','HH:mm:ss')) '  Exp. ' num2str(j_Exp) '. Finished LMA1 for channel ' num2str(j_chan)])
            end
            if do_save == 1
                save([folder_calc '\' file_name_3], 'Exp3_LMA1_MuaMusp')
                disp([char(datetime('now','Format','HH:mm:ss')) '  Saved ''Exp3_LMA1_MuaMusp'' into variable ''' file_name_3 ''''])
            end
        end
    end
end
if do_load == 0 && do_use_parpool == 1 % Terminate parpool
    delete(gcp)
end
if do_load == 0
    disp([char(datetime('now','Format','HH:mm:ss')) '  Finished calculating LMA 1 for all experiments'])
else
    disp([char(datetime('now','Format','HH:mm:ss')) '  LOADED LMA 1'])
end


show_fig = {'8a', '8b', '10'};
Blood_plotting(Exp1_Time_min,Exp2_Time_min,Exp3_Time_min,    Exp1_DTOF,Exp2_DTOF,Exp3_DTOF,   Time_ns,   Exp1_IRF,Exp2_IRF,Exp3_IRF,   Base_index,...
            Exp1_Mom,Exp2_Mom,Exp3_Mom,  Exp1_LMA1_MuaMusp,Exp2_LMA1_MuaMusp,Exp3_LMA1_MuaMusp,...
            Exp1_LMA2_Mua1Mua2,Exp2_LMA2_Mua1Mua2,Exp3_LMA2_Mua1Mua2,...
            Exp1_Conc,Exp2_Conc,Exp3_Conc,  Exp1_StO2,Exp2_StO2,Exp3_StO2, set_dim, show_fig )



%% Step 3.  Calculate moments, changes in moments, and perform 2-layered fitting for delMoments

% For saving/loading the calculated data:
file_name_1 = 'Blood_Exp1_LMA2_Mua1Mua2_vr1'; % Used: cut_lim_mom = [0.25 0.03];   Base_minutes = [14.6 15.6 31; 58.8 59.8 99];
file_name_2 = 'Blood_Exp2_LMA2_Mua1Mua2_vr1';
file_name_3 = 'Blood_Exp3_LMA2_Mua1Mua2_vr1';
% file_name_1 = 'Blood_Exp1_LMA2_Mua1Mua2_vr2';
% file_name_2 = 'Blood_Exp2_LMA2_Mua1Mua2_vr2';
% file_name_3 = 'Blood_Exp3_LMA2_Mua1Mua2_vr2';

do_load = 1; % 1 == load previously calculated results. 0 == Do calculation
% do_load = 0;

do_save = 1; % 1 == Save result in 'folder_calc'. Will use 'file_name' for the saved file
% do_save = 0;

% j_chan_all = 1:14;
j_chan_all = 3:11;
% j_chan_all = 3;

% do_use_parpool = 1;
do_use_parpool = 0;

if do_load == 0 && do_use_parpool == 1
    parpool
end

cut_lim_mom = [0.25 0.03]; % Used this for Publication 2023

L_known = [12 15 17];

delOptProp_all = [-7 -7 -8 0 0 0 0 40]; % [delMua1 delMua2 delMua3 delMusp1 delMusp2 delMusp3 L1 L2]

Options_LMA = [];

OptProp_base = zeros(6,1);

to_skip = 1;
% to_skip = 500; % To calculate for not every data point

for j_Exp = which_Exp % For which experiment to perform the calculation
% for j_Exp = 2

    if isequal(1, j_Exp) % FIRST EXPERIMENT
        delOptProp_all(7) = L_known(j_Exp);
        if do_load == 1
            load(file_name_1)
        else

            for j_chan = j_chan_all
                for j_base = 1:size(Base_index,1) % For each baseline
                    OptProp_base(1) = mean(Exp1_LMA1_MuaMusp.TwoL(Base_index(j_base,1,j_Exp):Base_index(j_base,2,j_Exp),j_chan, 1),1); % Mua  of Superficial layer (Measured on 2-layered)
                    OptProp_base(4) = mean(Exp1_LMA1_MuaMusp.TwoL(Base_index(j_base,1,j_Exp):Base_index(j_base,2,j_Exp),j_chan, 2),1); % Musp of Superficial layer (Measured on 2-layered)
                    OptProp_base(2:3) = mean(Exp1_LMA1_MuaMusp.Deep(Base_index(j_base,1,j_Exp):Base_index(j_base,2,j_Exp),j_chan, 1),1); % Mua  of Deep layer
                    OptProp_base(5:6) = mean(Exp1_LMA1_MuaMusp.Deep(Base_index(j_base,1,j_Exp):Base_index(j_base,2,j_Exp),j_chan, 2),1); % Musp of Deep layer
    
                    [ dtof, cut_ind_mom] = DTOF_filter( squeeze(mean(Exp1_DTOF.TwoL(Base_index(j_base,1,j_Exp):Base_index(j_base,2,j_Exp),j_chan,:),1)), Time_ns, bg_TwoL, 2, cut_lim_mom, 1); % Baseline
    
                    Mom_base = DTOF_CentralMom(Time_ns, dtof); % Baseline
    
                    if j_base == 1
                        j_cycle_all = 1:to_skip:Base_index(j_base,3,j_Exp);
                    else
                        j_cycle_all = 1+Base_index(j_base-1,3,j_Exp):to_skip:Base_index(j_base,3,j_Exp);
                    end
        %                     temp_min = 1:1:90; % For a quick test uncomment this and for figure(1000) below, will calculate for temp_min
        %                     j_cycle_all = zeros(1, length(temp_min));
        %                     for temp = 1:length(temp_min)
        %                         [~, j_cycle_all(temp)] = min(abs(temp_min(temp) - Exp1_Time_min));
        %                     end

                    if do_use_parpool == 0
                        for j_cycle = j_cycle_all
                            [ dtof, ~] = DTOF_filter( squeeze(Exp1_DTOF.TwoL(j_cycle,j_chan,:)), Time_ns, bg_TwoL, 2, cut_ind_mom, 3);
    
                            Mom_after = DTOF_CentralMom(Time_ns, dtof);
    
                            delMom = DTOF_DelMom(Mom_base, Mom_after);
    
                            [X_found, ~] = LMA_2_FittingDelMom(delOptProp_all, OptProp_base, n, rho, delMom,...
                                Mom_base, Time_ns, cut_ind_mom, Exp1_IRF_shifted.TwoL(:,j_chan), Options_LMA);
    
                            X_found(1) = X_found(1) + OptProp_base(1);
                            X_found(2) = X_found(2) + OptProp_base(2);
    
                            Exp1_LMA2_Mua1Mua2.TwoL(j_cycle,j_chan,:) = X_found(1:2);
                        end
        %                         figure(1000); % For a quick test
        %                         plot(Exp1_Time_min(j_cycle_all), Exp1_LMA2_Mua1Mua2.TwoL(j_cycle_all,j_chan,1) * 10,'x','Color','red')
        %                         hold on
        %                         plot(Exp1_Time_min(j_cycle_all), Exp1_LMA2_Mua1Mua2.TwoL(j_cycle_all,j_chan,2) * 10,'+','Color','blue')
        %                         grid on
        %                         plot(Exp1_Time_min, smooth(Exp1_LMA1_MuaMusp.Deep(:,j_chan,1) * 10, 10),'Color','black','LineWidth',1)
        %                         plot(Exp1_Time_min, smooth(Exp1_LMA1_MuaMusp.TwoL(:,j_chan,1) * 10, 10),'Color','green','LineWidth',1)
        %                         ylim([0.03 0.19])
        %                         xlim([-1.3 90])
        %                         plot([1 1]*Base_minutes(1,1),ylim+0.08,'Color','red')
        %                         plot([1 1]*Base_minutes(1,2),ylim+0.08,'Color','red')
        %                         plot([1 1]*Base_minutes(2,1),ylim+0.08,'Color','red')
        %                         plot([1 1]*Base_minutes(2,2),ylim+0.08,'Color','red')
        %                         drawnow
        %                         return
                    else
                    % The EXACT SAME CODE AS THE ABOVE, BUT WITH THE USE OF PARPOOL:
                        temp_res = zeros(length(j_cycle_all),2);
                        temp_dtof = squeeze(Exp1_DTOF.TwoL(j_cycle_all,j_chan,:));
                        temp_irf = Exp1_IRF_shifted.TwoL(:,j_chan);
                        parfor j = 1:length(j_cycle_all)
                            [ dtof, ~] = DTOF_filter( temp_dtof(j,:), Time_ns, bg_TwoL, 2, cut_ind_mom, 3);

                            Mom_after = DTOF_CentralMom(Time_ns, dtof);

                            delMom = DTOF_DelMom(Mom_base, Mom_after);

                            [X_found, ~] = LMA_2_FittingDelMom(delOptProp_all, OptProp_base, n, rho, delMom,...
                                Mom_base, Time_ns, cut_ind_mom, temp_irf, Options_LMA);

                            X_found(1) = X_found(1) + OptProp_base(1);
                            X_found(2) = X_found(2) + OptProp_base(2);

                            temp_res(j,:) = X_found(1:2);
                        end
                        Exp1_LMA2_Mua1Mua2.TwoL(j_cycle_all,j_chan,1) = temp_res(:,1);
                        Exp1_LMA2_Mua1Mua2.TwoL(j_cycle_all,j_chan,2) = temp_res(:,2);
                    end
                    disp([char(datetime('now','Format','HH:mm:ss')) '  Exp. ' num2str(j_Exp) '. Finished LMA2 for channel ' num2str(j_chan) ' FOR BASELINE ' num2str(j_base)])
                end
            end
            if do_save == 1
                save([folder_calc '\' file_name_1], 'Exp1_LMA2_Mua1Mua2')
                disp([char(datetime('now','Format','HH:mm:ss')) '  Saved ''Exp1_LMA2_Mua1Mua2'' into variable ''' file_name_1 ''''])
            end
        end
    

    elseif isequal(2, j_Exp) % SECOND EXPERIMENT (copy-pasted the above)
        delOptProp_all(7) = L_known(j_Exp);
        if do_load == 1
            load(file_name_2)
        else

            for j_chan = j_chan_all
                for j_base = 1:size(Base_index,1) % For each baseline
                    OptProp_base(1) = mean(Exp2_LMA1_MuaMusp.TwoL(Base_index(j_base,1,j_Exp):Base_index(j_base,2,j_Exp),j_chan, 1),1); % Mua  of Superficial layer (Measured on 2-layered)
                    OptProp_base(4) = mean(Exp2_LMA1_MuaMusp.TwoL(Base_index(j_base,1,j_Exp):Base_index(j_base,2,j_Exp),j_chan, 2),1); % Musp of Superficial layer (Measured on 2-layered)
                    OptProp_base(2:3) = mean(Exp2_LMA1_MuaMusp.Deep(Base_index(j_base,1,j_Exp):Base_index(j_base,2,j_Exp),j_chan, 1),1); % Mua  of Deep layer
                    OptProp_base(5:6) = mean(Exp2_LMA1_MuaMusp.Deep(Base_index(j_base,1,j_Exp):Base_index(j_base,2,j_Exp),j_chan, 2),1); % Musp of Deep layer
    
                    [ dtof, cut_ind_mom] = DTOF_filter( squeeze(mean(Exp2_DTOF.TwoL(Base_index(j_base,1,j_Exp):Base_index(j_base,2,j_Exp),j_chan,:),1)), Time_ns, bg_TwoL, 2, cut_lim_mom, 1); % Baseline
    
                    Mom_base = DTOF_CentralMom(Time_ns, dtof); % Baseline
    
                    if j_base == 1
                        j_cycle_all = 1:to_skip:Base_index(j_base,3,j_Exp);
                    else
                        j_cycle_all = 1+Base_index(j_base-1,3,j_Exp):to_skip:Base_index(j_base,3,j_Exp);
                    end

                    if do_use_parpool == 0
                        for j_cycle = j_cycle_all
                            [ dtof, ~] = DTOF_filter( squeeze(Exp2_DTOF.TwoL(j_cycle,j_chan,:)), Time_ns, bg_TwoL, 2, cut_ind_mom, 3);
    
                            Mom_after = DTOF_CentralMom(Time_ns, dtof);
    
                            delMom = DTOF_DelMom(Mom_base, Mom_after);
    
                            [X_found, ~] = LMA_2_FittingDelMom(delOptProp_all, OptProp_base, n, rho, delMom,...
                                Mom_base, Time_ns, cut_ind_mom, Exp2_IRF_shifted.TwoL(:,j_chan), Options_LMA);
    
                            X_found(1) = X_found(1) + OptProp_base(1);
                            X_found(2) = X_found(2) + OptProp_base(2);
    
                            Exp2_LMA2_Mua1Mua2.TwoL(j_cycle,j_chan,:) = X_found(1:2);
                        end
                    else
                    % The EXACT SAME CODE AS THE ABOVE, BUT WITH THE USE OF PARPOOL:
                        temp_res = zeros(length(j_cycle_all),2);
                        temp_dtof = squeeze(Exp2_DTOF.TwoL(j_cycle_all,j_chan,:));
                        temp_irf = Exp2_IRF_shifted.TwoL(:,j_chan);
                        parfor j = 1:length(j_cycle_all)
                            [ dtof, ~] = DTOF_filter( temp_dtof(j,:), Time_ns, bg_TwoL, 2, cut_ind_mom, 3);

                            Mom_after = DTOF_CentralMom(Time_ns, dtof);

                            delMom = DTOF_DelMom(Mom_base, Mom_after);

                            [X_found, ~] = LMA_2_FittingDelMom(delOptProp_all, OptProp_base, n, rho, delMom,...
                                Mom_base, Time_ns, cut_ind_mom, temp_irf, Options_LMA);

                            X_found(1) = X_found(1) + OptProp_base(1);
                            X_found(2) = X_found(2) + OptProp_base(2);

                            temp_res(j,:) = X_found(1:2);
                        end
                        Exp2_LMA2_Mua1Mua2.TwoL(j_cycle_all,j_chan,1) = temp_res(:,1);
                        Exp2_LMA2_Mua1Mua2.TwoL(j_cycle_all,j_chan,2) = temp_res(:,2);
                    end
                    disp([char(datetime('now','Format','HH:mm:ss')) '  Exp. ' num2str(j_Exp) '. Finished LMA2 for channel ' num2str(j_chan) ' FOR BASELINE ' num2str(j_base)])
                end
            end
            if do_save == 1
                save([folder_calc '\' file_name_2], 'Exp2_LMA2_Mua1Mua2')
                disp([char(datetime('now','Format','HH:mm:ss')) '  Saved ''Exp2_LMA2_Mua1Mua2'' into variable ''' file_name_2 ''''])
            end
        end
    
    elseif isequal(3, j_Exp) % THIRD EXPERIMENT (copy-pasted the above)
        delOptProp_all(7) = L_known(j_Exp);
        if do_load == 1
            load(file_name_3)
        else

            for j_chan = j_chan_all
                for j_base = 1:size(Base_index,1) % For each baseline
                    OptProp_base(1) = mean(Exp3_LMA1_MuaMusp.TwoL(Base_index(j_base,1,j_Exp):Base_index(j_base,2,j_Exp),j_chan, 1),1); % Mua  of Superficial layer (Measured on 2-layered)
                    OptProp_base(4) = mean(Exp3_LMA1_MuaMusp.TwoL(Base_index(j_base,1,j_Exp):Base_index(j_base,2,j_Exp),j_chan, 2),1); % Musp of Superficial layer (Measured on 2-layered)
                    OptProp_base(2:3) = mean(Exp3_LMA1_MuaMusp.Deep(Base_index(j_base,1,j_Exp):Base_index(j_base,2,j_Exp),j_chan, 1),1); % Mua  of Deep layer
                    OptProp_base(5:6) = mean(Exp3_LMA1_MuaMusp.Deep(Base_index(j_base,1,j_Exp):Base_index(j_base,2,j_Exp),j_chan, 2),1); % Musp of Deep layer
    
                    [ dtof, cut_ind_mom] = DTOF_filter( squeeze(mean(Exp3_DTOF.TwoL(Base_index(j_base,1,j_Exp):Base_index(j_base,2,j_Exp),j_chan,:),1)), Time_ns, bg_TwoL, 2, cut_lim_mom, 1); % Baseline
    
                    Mom_base = DTOF_CentralMom(Time_ns, dtof); % Baseline
    
                    if j_base == 1
                        j_cycle_all = 1:to_skip:Base_index(j_base,3,j_Exp);
                    else
                        j_cycle_all = 1+Base_index(j_base-1,3,j_Exp):to_skip:Base_index(j_base,3,j_Exp);
                    end

                    if do_use_parpool == 0
                        for j_cycle = j_cycle_all
                            [ dtof, ~] = DTOF_filter( squeeze(Exp3_DTOF.TwoL(j_cycle,j_chan,:)), Time_ns, bg_TwoL, 2, cut_ind_mom, 3);
    
                            Mom_after = DTOF_CentralMom(Time_ns, dtof);
    
                            delMom = DTOF_DelMom(Mom_base, Mom_after);
    
                            [X_found, ~] = LMA_2_FittingDelMom(delOptProp_all, OptProp_base, n, rho, delMom,...
                                Mom_base, Time_ns, cut_ind_mom, Exp3_IRF_shifted.TwoL(:,j_chan), Options_LMA);
    
                            X_found(1) = X_found(1) + OptProp_base(1);
                            X_found(2) = X_found(2) + OptProp_base(2);
    
                            Exp3_LMA2_Mua1Mua2.TwoL(j_cycle,j_chan,:) = X_found(1:2);
                        end
                    else
                    % The EXACT SAME CODE AS THE ABOVE, BUT WITH THE USE OF PARPOOL:
                        temp_res = zeros(length(j_cycle_all),2);
                        temp_dtof = squeeze(Exp3_DTOF.TwoL(j_cycle_all,j_chan,:));
                        temp_irf = Exp3_IRF_shifted.TwoL(:,j_chan);
                        parfor j = 1:length(j_cycle_all)
                            [ dtof, ~] = DTOF_filter( temp_dtof(j,:), Time_ns, bg_TwoL, 2, cut_ind_mom, 3);

                            Mom_after = DTOF_CentralMom(Time_ns, dtof);

                            delMom = DTOF_DelMom(Mom_base, Mom_after);

                            [X_found, ~] = LMA_2_FittingDelMom(delOptProp_all, OptProp_base, n, rho, delMom,...
                                Mom_base, Time_ns, cut_ind_mom, temp_irf, Options_LMA);

                            X_found(1) = X_found(1) + OptProp_base(1);
                            X_found(2) = X_found(2) + OptProp_base(2);

                            temp_res(j,:) = X_found(1:2);
                        end
                        Exp3_LMA2_Mua1Mua2.TwoL(j_cycle_all,j_chan,1) = temp_res(:,1);
                        Exp3_LMA2_Mua1Mua2.TwoL(j_cycle_all,j_chan,2) = temp_res(:,2);
                    end
                    disp([char(datetime('now','Format','HH:mm:ss')) '  Exp. ' num2str(j_Exp) '. Finished LMA2 for channel ' num2str(j_chan) ' FOR BASELINE ' num2str(j_base)])
                end
            end
            if do_save == 1
                save([folder_calc '\' file_name_3], 'Exp3_LMA2_Mua1Mua2')
                disp([char(datetime('now','Format','HH:mm:ss')) '  Saved ''Exp3_LMA2_Mua1Mua2'' into variable ''' file_name_3 ''''])
            end
        end
    end
end
if do_load == 0 && do_use_parpool == 1 % Terminate parpool
    delete(gcp)
end
if do_load == 0
    disp([char(datetime('now','Format','HH:mm:ss')) '  Finished calculating LMA 2 for all experiments'])
else
    disp([char(datetime('now','Format','HH:mm:ss')) '  LOADED LMA 2'])
end


show_fig = {'9a', '9b', '9c'};
Blood_plotting(Exp1_Time_min,Exp2_Time_min,Exp3_Time_min,    Exp1_DTOF,Exp2_DTOF,Exp3_DTOF,   Time_ns,   Exp1_IRF,Exp2_IRF,Exp3_IRF,   Base_index,...
            Exp1_Mom,Exp2_Mom,Exp3_Mom,  Exp1_LMA1_MuaMusp,Exp2_LMA1_MuaMusp,Exp3_LMA1_MuaMusp,...
            Exp1_LMA2_Mua1Mua2,Exp2_LMA2_Mua1Mua2,Exp3_LMA2_Mua1Mua2,...
            Exp1_Conc,Exp2_Conc,Exp3_Conc,  Exp1_StO2,Exp2_StO2,Exp3_StO2, set_dim, show_fig )



%% Step 4.  Calculate concentrations of Oxy and Deoxy (and water, the Mua of which is then subtracted)

% For saving/loading the calculated data:
file_name_1 = 'Blood_Exp1_Conc_StO2_vr1';
file_name_2 = 'Blood_Exp2_Conc_StO2_vr1';
file_name_3 = 'Blood_Exp3_Conc_StO2_vr1';
% file_name_1 = 'Blood_Exp1_Conc_StO2_vr2';
% file_name_2 = 'Blood_Exp2_Conc_StO2_vr2';
% file_name_3 = 'Blood_Exp3_Conc_StO2_vr2';

% do_load = 1; % 1 == load previously calculated results. 0 == Do calculation
do_load = 0;

% do_save = 1; % 1 == Save result in 'folder_calc'. Will use 'file_name' for the saved file
do_save = 0;

chan_use = 3:11; % Used for Publication 2023


% Correct Mua spectra of Water from literature,    See  Ink_Nominal_mua.m
mua_Water_all = load('matcher94_nir_water_37.txt');
mua_Water_all(:,2) = mua_Water_all(:,2) * log(10); % Absorption spectrum of Water from Matcher [OD cm-1]
mua_Water = interp1(mua_Water_all(:,1),mua_Water_all(:,2),wavelengths(chan_use)); % At chosen wavelengths

% Mua spectra of Water from another source:
% mua_Water = GetExtinctions(wavelengths(chan_use), 2); % [HbO Hb H2O lipid aa3]    [OD cm-1] for Water 
% mua_Water = mua_Water(:,3)';

% Correct Molar Absorption Coefficients for Oxy and Deoxy:
Ext_coeff = [];
Ext_coeff(:,1) = wavelengths(chan_use);
wh_source = 1; % Gratzer   % < Used for Publication 2023
% wh_source = 2; % Moaveni
% wh_source = 3; % Takatani
Ext_coeff(:,2:end+5) = GetExtinctions(Ext_coeff(:,1), wh_source);
Ext_coeff(:,2:3) = Ext_coeff(:,2:3) .* log(10);  % Molar Absorption Coefficient   [OD M-1 cm-1]   See GetExtinctions.m


[UMIX,svMIX,VMIX]=svd(Ext_coeff(:,2:3));
isv=pinv(svMIX);

mua_to_subtract = mua_Water' * (3366 + 180) / 3600; % Assuming SMOFlipid has the same Mua as Water and ignoring Yeast
% mua_to_subtract = 0;


for j_Exp = which_Exp % For which experiment to perform the calculation
% for j_Exp = 1

    if isequal(1, j_Exp) % FIRST EXPERIMENT
        if do_load == 1
            load(file_name_1)
        else
            % Calculate Concentrations:
            for j_cycle = 1:size(Exp1_DTOF.TwoL,1)
                Exp1_Conc.TwoL_LMA1(j_cycle,1:2) = VMIX*isv*UMIX'*(  squeeze(Exp1_LMA1_MuaMusp.TwoL(j_cycle,chan_use,1)*10)' -  mua_to_subtract) * 10^6; % micro Molar

                Exp1_Conc.Deep_LMA1(j_cycle,1:2) = VMIX*isv*UMIX'*(  squeeze(Exp1_LMA1_MuaMusp.Deep(j_cycle,chan_use,1)*10)' -  mua_to_subtract) * 10^6;
                
                Exp1_Conc.TwoL_LMA2(j_cycle,1:2,1) = VMIX*isv*UMIX'*(  squeeze(Exp1_LMA2_Mua1Mua2.TwoL(j_cycle,chan_use,1)*10)' -  mua_to_subtract) * 10^6;

                Exp1_Conc.TwoL_LMA2(j_cycle,1:2,2) = VMIX*isv*UMIX'*(  squeeze(Exp1_LMA2_Mua1Mua2.TwoL(j_cycle,chan_use,2)*10)' -  mua_to_subtract) * 10^6;
            end

            % Calculate StO2:
            for j_cycle = 1:size(Exp1_DTOF.TwoL,1)
                Exp1_StO2.TwoL_LMA1(j_cycle) = Exp1_Conc.TwoL_LMA1(j_cycle,1) / (Exp1_Conc.TwoL_LMA1(j_cycle,1) + Exp1_Conc.TwoL_LMA1(j_cycle,2)) * 100;

                Exp1_StO2.Deep_LMA1(j_cycle) = Exp1_Conc.Deep_LMA1(j_cycle,1) / (Exp1_Conc.Deep_LMA1(j_cycle,1) + Exp1_Conc.Deep_LMA1(j_cycle,2)) * 100;
                
                Exp1_StO2.TwoL_LMA2(j_cycle,1) = Exp1_Conc.TwoL_LMA2(j_cycle,1,1) / (Exp1_Conc.TwoL_LMA2(j_cycle,1,1) + Exp1_Conc.TwoL_LMA2(j_cycle,2,1)) * 100;

                Exp1_StO2.TwoL_LMA2(j_cycle,2) = Exp1_Conc.TwoL_LMA2(j_cycle,1,2) / (Exp1_Conc.TwoL_LMA2(j_cycle,1,2) + Exp1_Conc.TwoL_LMA2(j_cycle,2,2)) * 100;
            end

            disp([char(datetime('now','Format','HH:mm:ss')) '  Exp. ' num2str(j_Exp) '. Finished Calculating Concentrations and StO2'])
            if do_save == 1
                save([folder_calc '\' file_name_1], 'Exp1_Conc', 'Exp1_StO2')
                disp([char(datetime('now','Format','HH:mm:ss')) '  Saved ''Exp1_Conc'' and ''Exp1_StO2'' into variable ''' file_name_1 ''''])
            end
        end


    elseif isequal(2, j_Exp) % SECOND EXPERIMENT (copy-pasted the above)
        if do_load == 1
            load(file_name_2)
        else
            % Calculate Concentrations:
            for j_cycle = 1:size(Exp2_DTOF.TwoL,1)
                Exp2_Conc.TwoL_LMA1(j_cycle,1:2) = VMIX*isv*UMIX'*(  squeeze(Exp2_LMA1_MuaMusp.TwoL(j_cycle,chan_use,1)*10)' -  mua_to_subtract) * 10^6; % micro Molar

                Exp2_Conc.Deep_LMA1(j_cycle,1:2) = VMIX*isv*UMIX'*(  squeeze(Exp2_LMA1_MuaMusp.Deep(j_cycle,chan_use,1)*10)' -  mua_to_subtract) * 10^6;
                
                Exp2_Conc.TwoL_LMA2(j_cycle,1:2,1) = VMIX*isv*UMIX'*(  squeeze(Exp2_LMA2_Mua1Mua2.TwoL(j_cycle,chan_use,1)*10)' -  mua_to_subtract) * 10^6;

                Exp2_Conc.TwoL_LMA2(j_cycle,1:2,2) = VMIX*isv*UMIX'*(  squeeze(Exp2_LMA2_Mua1Mua2.TwoL(j_cycle,chan_use,2)*10)' -  mua_to_subtract) * 10^6;
            end

            % Calculate StO2:
            for j_cycle = 1:size(Exp2_DTOF.TwoL,1)
                Exp2_StO2.TwoL_LMA1(j_cycle) = Exp2_Conc.TwoL_LMA1(j_cycle,1) / (Exp2_Conc.TwoL_LMA1(j_cycle,1) + Exp2_Conc.TwoL_LMA1(j_cycle,2)) * 100;

                Exp2_StO2.Deep_LMA1(j_cycle) = Exp2_Conc.Deep_LMA1(j_cycle,1) / (Exp2_Conc.Deep_LMA1(j_cycle,1) + Exp2_Conc.Deep_LMA1(j_cycle,2)) * 100;
                
                Exp2_StO2.TwoL_LMA2(j_cycle,1) = Exp2_Conc.TwoL_LMA2(j_cycle,1,1) / (Exp2_Conc.TwoL_LMA2(j_cycle,1,1) + Exp2_Conc.TwoL_LMA2(j_cycle,2,1)) * 100;

                Exp2_StO2.TwoL_LMA2(j_cycle,2) = Exp2_Conc.TwoL_LMA2(j_cycle,1,2) / (Exp2_Conc.TwoL_LMA2(j_cycle,1,2) + Exp2_Conc.TwoL_LMA2(j_cycle,2,2)) * 100;
            end

            disp([char(datetime('now','Format','HH:mm:ss')) '  Exp. ' num2str(j_Exp) '. Finished Calculating Concentrations and StO2'])
            if do_save == 1
                save([folder_calc '\' file_name_2], 'Exp1_Conc', 'Exp1_StO2')
                disp([char(datetime('now','Format','HH:mm:ss')) '  Saved ''Exp2_Conc'' and ''Exp2_StO2'' into variable ''' file_name_2 ''''])
            end
        end

    elseif isequal(3, j_Exp) % THIRD EXPERIMENT (copy-pasted the above)
        if do_load == 1
            load(file_name_3)
        else
            % Calculate Concentrations:
            for j_cycle = 1:size(Exp3_DTOF.TwoL,1)
                Exp3_Conc.TwoL_LMA1(j_cycle,1:2) = VMIX*isv*UMIX'*(  squeeze(Exp3_LMA1_MuaMusp.TwoL(j_cycle,chan_use,1)*10)' -  mua_to_subtract) * 10^6; % micro Molar

                Exp3_Conc.Deep_LMA1(j_cycle,1:2) = VMIX*isv*UMIX'*(  squeeze(Exp3_LMA1_MuaMusp.Deep(j_cycle,chan_use,1)*10)' -  mua_to_subtract) * 10^6;
                
                Exp3_Conc.TwoL_LMA2(j_cycle,1:2,1) = VMIX*isv*UMIX'*(  squeeze(Exp3_LMA2_Mua1Mua2.TwoL(j_cycle,chan_use,1)*10)' -  mua_to_subtract) * 10^6;

                Exp3_Conc.TwoL_LMA2(j_cycle,1:2,2) = VMIX*isv*UMIX'*(  squeeze(Exp3_LMA2_Mua1Mua2.TwoL(j_cycle,chan_use,2)*10)' -  mua_to_subtract) * 10^6;
            end

            % Calculate StO2:
            for j_cycle = 1:size(Exp3_DTOF.TwoL,1)
                Exp3_StO2.TwoL_LMA1(j_cycle) = Exp3_Conc.TwoL_LMA1(j_cycle,1) / (Exp3_Conc.TwoL_LMA1(j_cycle,1) + Exp3_Conc.TwoL_LMA1(j_cycle,2)) * 100;

                Exp3_StO2.Deep_LMA1(j_cycle) = Exp3_Conc.Deep_LMA1(j_cycle,1) / (Exp3_Conc.Deep_LMA1(j_cycle,1) + Exp3_Conc.Deep_LMA1(j_cycle,2)) * 100;
                
                Exp3_StO2.TwoL_LMA2(j_cycle,1) = Exp3_Conc.TwoL_LMA2(j_cycle,1,1) / (Exp3_Conc.TwoL_LMA2(j_cycle,1,1) + Exp3_Conc.TwoL_LMA2(j_cycle,2,1)) * 100;

                Exp3_StO2.TwoL_LMA2(j_cycle,2) = Exp3_Conc.TwoL_LMA2(j_cycle,1,2) / (Exp3_Conc.TwoL_LMA2(j_cycle,1,2) + Exp3_Conc.TwoL_LMA2(j_cycle,2,2)) * 100;
            end

            disp([char(datetime('now','Format','HH:mm:ss')) '  Exp. ' num2str(j_Exp) '. Finished Calculating Concentrations and StO2'])
            if do_save == 1
                save([folder_calc '\' file_name_3], 'Exp1_Conc', 'Exp1_StO2')
                disp([char(datetime('now','Format','HH:mm:ss')) '  Saved ''Exp3_Conc'' and ''Exp3_StO2'' into variable ''' file_name_3 ''''])
            end
        end
    end
end
if do_load == 0
    disp([char(datetime('now','Format','HH:mm:ss')) '  Finished calculating Conc and StO2 for all experiments'])
else
    disp([char(datetime('now','Format','HH:mm:ss')) '  LOADED Conc and StO2'])
end


show_fig = {'12a','12b','12c','13a','13b','13c'};
Blood_plotting(Exp1_Time_min,Exp2_Time_min,Exp3_Time_min,    Exp1_DTOF,Exp2_DTOF,Exp3_DTOF,   Time_ns,   Exp1_IRF,Exp2_IRF,Exp3_IRF,   Base_index,...
            Exp1_Mom,Exp2_Mom,Exp3_Mom,  Exp1_LMA1_MuaMusp,Exp2_LMA1_MuaMusp,Exp3_LMA1_MuaMusp,...
            Exp1_LMA2_Mua1Mua2,Exp2_LMA2_Mua1Mua2,Exp3_LMA2_Mua1Mua2,...
            Exp1_Conc,Exp2_Conc,Exp3_Conc,  Exp1_StO2,Exp2_StO2,Exp3_StO2, set_dim, show_fig )