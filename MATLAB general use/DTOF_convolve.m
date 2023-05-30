% Written by Aleh Sudakou (https://www.researchgate.net/profile/Aleh-Sudakou)
% Nałęcz Institute of Biocybernetics and Biomedical Engineering, Polish Academy of Sciences 
% This is one of the codes shared on  https://github.com/asudakou/Analyzing_TD-NIRS
% Last updated: 20 May 2023


%% Description:
% This function convolves two DTOFs, e.g. TPSF (deconvolved or theoretical DTOF) and IRF.
% 
% For now, convolution works only if the time channels are the same for IRF and TPSF. 
% If they are not – first rescale them to make them the same, or write a separate function that will do convolution for such cases.
% To account for offset, use "interp1" built-in function
%
% The function assumes that the first time bin is time 0, such that
% colvinving with "IRF = [1 0 0 0]" does not change TPSF. Whereas
% convolving with "IRF = [0 0 0 1]" shifts TPSF by 3 time channels to the right


function [TPSF_conv] = DTOF_convolve(TPSF, IRF_shifted)


do_this = 0; % This must always be 0
if do_this == 1
    %% Run this section to see Example (Ctrl + Enter)
    % Generate first DTOF:
    OptProp_all_example = [0.01 -1 -1 1 -1 -1 0 0]; % Homogeneous optical properties, Mua=0.01 and Musp=1
    [R1_example, time_ns_example] = DTOF_generate_Liemert(OptProp_all_example, 1.4, 30, 0.1:0.1:7);
    
    % Generate second DTOF (list of zeros and 1):
    R2_example = zeros(21,1);
%     R2_example(1) = 1; % Does nothing to TPSF
    R2_example(2) = 1; % Shifts by 1 time channel to the right
%     R2_example(5) = 1; % Shifts by 4 time channels to the right
    
    R3_example = DTOF_convolve(R1_example, R2_example); % Call the current function
    
    figure(5); clf % Plot First DTOF, Second DTOF, and convolved DTOF
    [~, ~] = DTOF_filter( R1_example, time_ns_example, 0, 0, 0, 0, {5, 'red','x'});
    [~, ~] = DTOF_filter( R2_example, time_ns_example(1:length(R2_example)), 0, 0, 0, 0, {5, 'black','o'});
    [~, ~] = DTOF_filter( R3_example, time_ns_example, 0, 0, 0, 0, {5, 'blue','+'});
    xlim([-0.3 7.3])
    clear OptProp_all_example R1_example R2_example R3_example time_ns_example;
end


%% Compare two sets of time channels and rescale if needed
% Not in this version of the function. 


%% Do convolution:
% I suggest watching this YouTube tutorial: 
% "2021 How to Perform a Convolution in MATLAB | MATLAB Tutorial" https://www.youtube.com/watch?v=hcyy144Gu60

TPSF_conv = conv(IRF_shifted, TPSF);
% time_ns_conv = (0:length(TPSF_convolved)-1) .* dt_DTOF + time_ns(1);

TPSF_conv = TPSF_conv(1:1+length(TPSF)-1);
% time_ns_conv = time_ns_conv(1:1+length(TPSF)-1);


end