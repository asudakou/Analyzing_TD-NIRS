% Written by Aleh Sudakou (https://www.researchgate.net/profile/Aleh-Sudakou)
% Nałęcz Institute of Biocybernetics and Biomedical Engineering, Polish Academy of Sciences 
% This is one of the codes shared on  https://github.com/asudakou/Analyzing_TD-NIRS
% Last updated: 20 May 2023


%% Description:
% For given changes in Absorption and/or Scattering in up to 3 layers, 
% generate two DTOFs (using Liemert's code) and calculate the corresponding changes in moments.
% 
% This function is used by 'LMA_2_RearrangeInput.m'
% For details, see         'LMA_2_FittingDelMom.m'


function [delMom] = DelOptProp_to_DelMom(delOptProp_all, OptProp_base, n, rho, time_ns, cut_ind_mom, irf_shifted)


do_this = 0; % This must always be 0
if do_this == 1
    %% Run this section to see Example (Ctrl + Enter)
    delOptProp_all_example = zeros(8,1);
    delOptProp_all_example([1 2 7]) = [0.005 -0.005 15]; % delMua1, delMua2, and L1  -  Changes in optical properties and thicknesses
    OptProp_base_example = [0.01 0.01 0.01 1 1 1]; % [Mua1 Mua2 Mua3 Musp1 Musp2 Musp3]  -  Baseline optical properties
    time_ns_example = 0.01:0.01:7; % Time channels
    cut_ind_mom_example = [1 700]; % Moments depend on cut region, so we need to specify it
    [delMom_example] = DelOptProp_to_DelMom(delOptProp_all_example, OptProp_base_example, 1.4, 30, time_ns_example, cut_ind_mom_example, -1); % Call the current function
    disp('delMom is:'); disp(delMom_example')
    clear delOptProp_all_example OptProp_base_example time_ns_example cut_ind_mom_example delMom_example;
end


% INPUTS:
% For details, see 'LMA_2_FittingDelMom.m'

if length(delOptProp_all) ~= 8
    error('variable ""delOptProp_all"" must be length 8')
end

if length(OptProp_base) ~= 6
    error('variable ""OptProp_base"" must be length 6')
end

if rho < 5
    error('rho (source-detector distance) should be in mm')
end

OptProp_base(7:8) = delOptProp_all(7:8);

OptProp_after = OptProp_base; % Optical properties after a change in optical properties
for j = 1:6
    OptProp_after(j) = OptProp_after(j) + delOptProp_all(j);
end


%% Main body

% Generate DTOF at baseline
[R_base_conv, ~] = DTOF_generate_Liemert(OptProp_base, n, rho, time_ns, irf_shifted);

% Generate DTOF after delOptProp_all
[R_after_conv, ~] = DTOF_generate_Liemert(OptProp_after, n, rho, time_ns, irf_shifted);

% Calculate Moments
R_mom_base  = DTOF_CentralMom(time_ns(cut_ind_mom(1):cut_ind_mom(2)), R_base_conv(cut_ind_mom(1):cut_ind_mom(2)));

R_mom_after = DTOF_CentralMom(time_ns(cut_ind_mom(1):cut_ind_mom(2)), R_after_conv(cut_ind_mom(1):cut_ind_mom(2)));

% Calculate changes in Moments due to change in optical properties
delMom = DTOF_DelMom(R_mom_base, R_mom_after);


end