% Written by Aleh Sudakou (https://www.researchgate.net/profile/Aleh-Sudakou)
% Nałęcz Institute of Biocybernetics and Biomedical Engineering, Polish Academy of Sciences 
% This is one of the codes shared on  https://github.com/asudakou/Analyzing_TD-NIRS
% Last updated: 20 May 2023


%% Description:
% Using LMA algorithm (MATLAB's built-in function), search for the changes in
% optical properties that best correspond to the changes in moments. 
%
% The initial guess is calculated using the moments method based on the sensitivity factors. 
% Publication:
% A. Liebert, H. Wabnitz, J. Steinbrink, H. Obrig, M. Möller, R. Macdonald, A. Villringer, and H. Rinneberg, 
% "Time-Resolved Multidistance Near-Infrared Spectroscopy of the Adult Head: Intracerebral and Extracerebral 
% Absorption Changes from Moments of Distribution of Times of Flight of Photons," Appl. Opt. 43, 3037-3047 (2004).
% DOI:  https://doi.org/10.1364/AO.43.003037
%
% Required files: 
% Matlab Optimization Toolbox, for using in-built function 'lsqcurvefit'
% 
% Six functions:  "LMA_2_RearrangeInput.m"
%                 "DelOptProp_to_DelMom.m"
%                 "DTOF_generate_Liemert.m"
%                 "DTOF_convolve.m"
%                 "DTOF_CentralMom.m"
%                 "DTOF_delMom.m"
% 
% Inputs:
% * delOptProp_all  :  8 knowns and unknowns [delMua1 delMua2 delMua3 delMusp1 delMusp2 delMusp3 L1 L2]     [Units mm-1] and [mm for L1 and L2]
%   Value -7  means this parameter is unknown
%   Value -8  means to take the same value as for the above layer (cannot be assigned to the superficial layer or L)
%   Any other value means this is the value of the parameter.
%   (Set to 0 if we know that the optical properties didn't change in that layer)
%   (For a homogeneous medium, set to -8 both second and third layers ; values of L1 and L2 won't matter)
%   (For a two-layered medium, set to -8 the third layer; value of L2 won't matter)
%   Examples of determining different parameters:
%   delOptProp_all  = [-7 -8 -8 -7 -8 -8 0  0]  % delMua and delMusp,  Homogeneous medium
%   delOptProp_all  = [-7 -7 -8  0 -8 -8 10 0]  % delMua1 and delMua2, Two-layered medium 
%   delOptProp_all  = [-7  0 -7  0  0  0 10 10] % delMua1 and delMua3, Three-layered medium 
%   delOptProp_all  = [-7  0  0 -8 -8 -8 10 10] % delMua1 and delMusp, Three-layered medium, two unknowns
%   delOptProp_all  = [0  -7  0  0  0  0 -7 0]  % delMua2 and L,       Two-layered medium
%
% * OptProp_all     :  6 optical properties during baseline: [Mua1 Mua2 Mua3 Musp1 Musp2 Musp3]      [Units mm-1]
% 
% * n         :  refractive index, e.g. 1.4 or [1.4 1.4 1.4]
% * rho       :  source-detector distance in mm
% * delMom    :  changes in moments between baseline and after
% * Mom       :  moments during baseline (used for the first guess, if not calculating first guess then don't need this input)
% * time_ns   :  time channels of measured DTOF
% * cut_ind_mom :  indices of time channels used for calculating moments
% 
% * irf_shifted : variable containing IRF (having been moved by offset). Can be set to -1 if not used
% 
% Note, changes in first three moments don't depend on IRF, but they do depend on the time channels used,
% so we convolve with IRF to 'stretch' theoretical DTOFs, for using the same time-channels for calculating moments as for the measured DTOFs
% 
% * varargin              :  possible optional inputs
%   varargin.which_mom    :  which moments to use, e.g. [2 3] for Mean time of flight and Variance
%   varargin.weight_noise :  divide each moment by this weight. Set to [1 1 1 1] if don't want to add weights
%   varargin.init_guess   :   Use this as initial guess or if set to 'no' then use typicalX as initial guess
%   varargin.do_disp      :  [1 1 1] to display or plot: 1) Initial guess, 2) 'iter-detailed' LMA, and 3) Determined values


function [opt_prop_found, Resid] = LMA_2_FittingDelMom(delOptProp_all, OptProp_base, n, rho, delMom, Mom, time_ns, cut_ind_mom, irf_shifted, varargin)


do_this = 0; % This must always be 0
if do_this == 1
    %% Run this section to see Example (Ctrl + Enter)
    which_example = 1; % Set to 1, 2, 3, or 4

    % 1) Generate DTOFs
    OptProp_base_example = [0.01 0.01 0.01 1 1 1 5 5]; % Homogeneous, Mua=0.01 and Musp=1
    
    OptProp_new_example = OptProp_base_example;
    if which_example == 1
        OptProp_new_example([1 3]) = OptProp_new_example([1 3]) + 0.01; % Increase Mua in 1st and 3rd layers
    elseif which_example == 2
        OptProp_new_example(1) = OptProp_new_example(1) + 0.01; % Increase Mua in 1st layer
    elseif which_example == 3
        OptProp_new_example(1:3) = OptProp_new_example(1:3) + 0.01; % Increase Mua in all 3 layers
    elseif which_example == 4
        OptProp_new_example(1:3) = OptProp_new_example(1:3) + 0.01; % Increase Mua and Musp in all 3 layers
        OptProp_new_example(4:6) = OptProp_new_example(4:6) + 1;
    end
    [R1_example, time_ns_example] = DTOF_generate_Liemert(OptProp_base_example, 1.4, 30, -1);
    [R2_example, ~] =       DTOF_generate_Liemert(OptProp_new_example, 1.4, 30, -1);

    % 2) Choose cut-off region, and calculate moments and changes in moments
    [ ~, cut_ind_mom_example] = DTOF_filter( R1_example, time_ns_example, 0, 0, [0.01 0.01], 1); % At 1%
    [Mom1_example] = DTOF_CentralMom(time_ns_example(cut_ind_mom_example(1):cut_ind_mom_example(2)), R1_example(cut_ind_mom_example(1):cut_ind_mom_example(2)));
    [Mom2_example] = DTOF_CentralMom(time_ns_example(cut_ind_mom_example(1):cut_ind_mom_example(2)), R2_example(cut_ind_mom_example(1):cut_ind_mom_example(2)));
    delMom_example = DTOF_DelMom(Mom1_example, Mom2_example);

    % 3) Call the current function to determine changes in optical properties
    OptProp_base_example = OptProp_base_example(1:6); % known baseline Mua and Musp
    
    options_LMA2_example = [];
    options_LMA2_example.which_mom_use = [2 3]; % Will use Mean time of flight and Variance
    options_LMA2_example.weight = [1 1 1]; % Since we have an inverse-crime here, weight makes no difference

    options_LMA2_example.do_disp = [1 0 1]; % What information to display/plot
                                            % 1) Initial guess, 2) 'iter-detailed' LMA, and 3) Determined values

    if which_example == 1
        delOptProp_all_example = [-7  0 -7 0 0 0 5 5]; % set delMua1 and delMua3 as 2 unknowns, three-layered

    elseif which_example == 2
        delOptProp_all_example = [-7 -7 -8 0 0 0 5 5]; % set delMua1 and delMua2 as 2 unknowns, two-layered
        
    elseif which_example == 3
        delOptProp_all_example = [-7 -8 -8 0 0 0 0 0]; % set delMua1 as 1 unknown, homogeneous

    elseif which_example == 4
        delOptProp_all_example = [-7 -8 -8 -7 -8 -8 0 0]; % set delMua and delMusp as 2 unknown, homogeneous

    end
    
    [~, ~] = LMA_2_FittingDelMom(delOptProp_all_example, OptProp_base_example, 1.4, 30, delMom_example, Mom1_example, time_ns_example, cut_ind_mom_example, -1, options_LMA2_example);
    clear which_example OptProp_base_example OptProp_new_example time_ns_example R1_example R2_example Mom1_example Mom2_example...
            delMom_example options_LMA2_example delOptProp_all_example; 
    return
end


% INPUTS:
if length(delOptProp_all) ~= 8
    error('variable ""delOptProp_all"" must be length 8')
end
if length(OptProp_base) ~= 6
    error('variable ""OptProp_all"" must be length 6')
end
if rho < 5 || rho > 50
    error('Check source-detector distance (rho)')
end
if n(1) < 0.5 || n(1) > 2.5
    error('Check refractive index (n). It is too unusual for it to be below 0.5 or above 2.5 (To continue anyway, comment out this error)')
end
if size(delMom,2) ~= 1
    delMom = delMom';
end
if sum(delMom==0) == length(delMom) % Don't need to calculate anything if delMom are all zeros
    disp([char(datetime('now','Format','HH:mm:ss')) '  Returning zeros because inputted changes of moments are all zeros'])
    opt_prop_found = zeros(sum(delOptProp_all==-7),1);
    Resid = 0;
    return
end
if cut_ind_mom(1) < 0 || cut_ind_mom(2) < 0
    opt_prop_found = nan(size(sum(delOptProp_all==-7)),1);
    Resid = nan(1);
    disp([char(datetime('now','Format','HH:mm:ss')) '  WARNING: wrong cut indices, output result is set to NaN.'])
    return
end


% Default values:
temp = [0.01 0.01 0.01 1 1 1 10 10]; % (Default) [delMua1 delMua2 delMua3 delMusp1 delMusp2 delMusp3 L1 L2]
typical_delOptProp_all = temp * 0.10; % (Default) Ten percent increase
lbound = [-0.04 -0.04 -0.04 -10 -10 -10 1 1]; % (Default) Lower boundary values 
rbound = [ 0.04  0.04  0.04  10  10  10 30 30]; % (Default) Upper boundary values

typical_X = typical_delOptProp_all(delOptProp_all == -7); % Select only for the unknown properties
lbound = lbound(delOptProp_all == -7);
rbound = rbound(delOptProp_all == -7);

which_mom_use = [2 3]; % (Default) delm1 and delV
% which_mom_use = [1 2 3 4]; % For using all moments

weight_noise = [1 1 1]; % (Default)
% weight_noise = [1/10 1 1]; % Might be good to reduce the weight of attenuation
                            % The correct approach is to reduce by the uncertainties of moments

do_disp = [0 0 0]; % 1) Initial guess, 2) 'iter-detailed' LMA, and 3) Determined values

initial_guess = -10; % Will calculate the initial guess
% initial_guess = typical_X; % Uncomment this to not calculate the initial guess

% Optional inputs:
if ~isempty(varargin)
    if length(varargin) == 1
        if isfield(varargin{1},'do_disp')
            do_disp = varargin{1}.do_disp;
        end
        if isfield(varargin{1},'which_mom_use')
            which_mom_use = varargin{1}.which_mom_use;
        end
        if isfield(varargin{1},'initial_guess')
            initial_guess = varargin{1}.initial_guess;
            if isequal(initial_guess, 'no')
                initial_guess = typical_X;
            end
        end
        if isfield(varargin{1},'lbound')
            lbound = varargin{1}.lbound;
        end
        if isfield(varargin{1},'rbound')
            rbound = varargin{1}.rbound;
        end
        if isfield(varargin{1},'weight_noise')
            weight_noise = varargin{1}.weight_noise;
        end
    else
        error('Wrong number of input arguments')
    end
end


% Calculate initial guess if it wasn't specified.
% To avoid this step, set "do_this" to 0 or uncomment the above line "initial_guess = typical_X;"
do_this = 1;
% do_this = 0;
if isequal(initial_guess, -10) && do_this == 1
    
    if sum(delOptProp_all(4:6)==-7) == 0 % Wont calculate initial guess if determining changes in scattering

        which_layer = find(delOptProp_all(1:3) == -7);

        if delOptProp_all(7) == -7 % Thicknesses of both layers
            L_temp = typical_delOptProp_all(7:8);
        else
            L_temp = delOptProp_all(7:8);
        end

        K = [1, 0, 0; 0, Mom(3), Mom(4); 0, Mom(4), Mom(5) - Mom(3)^2]; % K
        K = K ./ Mom(1);

        if length(which_layer) == 3 
            % Not calculating for changes in 3 unknowns. It can be done,
            % but will take too many lines

        elseif length(which_layer) == 2 % Finding changes in 2 layers
            p = 5 / 100; % + 5 percent
            temp_base1 = OptProp_base; % Baseline
            temp_base1(which_layer(1)) = -7;
            temp_base1(delOptProp_all==-8) = -8;
            temp_base1(7:8) = L_temp;
                temp_base2 = OptProp_base;
                temp_base2(which_layer(2)) = -7;
                temp_base2(delOptProp_all==-8) = -8;
                temp_base2(7:8) = L_temp;
            temp_mom = zeros(5, 2);
            Sens_Fact = zeros(3,2); % 3 moments, 2 layers
            for j = 1:2 % Run two times, first time with +-"p" percent change and second time with more realistic changes
                if j == 1
                    temp_change = OptProp_base(which_layer(1)); % 1st DTOF (baseline)
                        [temp_R] = LMA_1_RearrangeInput(temp_change, temp_base1, n, rho, time_ns, [1 length(time_ns)], -1, 1);
                        temp_mom_base = DTOF_CentralMom(time_ns(cut_ind_mom(1):cut_ind_mom(2)),temp_R(cut_ind_mom(1):cut_ind_mom(2)));

                    temp_change = OptProp_base(which_layer(1)) * (1+p); % 2nd DTOF
                        [temp_R] = LMA_1_RearrangeInput(temp_change, temp_base1, n, rho, time_ns, [1 length(time_ns)], -1, 1);
                        temp_mom(:,1) = DTOF_CentralMom(time_ns(cut_ind_mom(1):cut_ind_mom(2)),temp_R(cut_ind_mom(1):cut_ind_mom(2)));

                    temp_change = OptProp_base(which_layer(2)) * (1+p); % 3rd DTOF
                        [temp_R] = LMA_1_RearrangeInput(temp_change, temp_base2, n, rho, time_ns, [1 length(time_ns)], -1, 1);
                        temp_mom(:,2) = DTOF_CentralMom(time_ns(cut_ind_mom(1):cut_ind_mom(2)),temp_R(cut_ind_mom(1):cut_ind_mom(2)));

                    temp = DTOF_DelMom(temp_mom_base, temp_mom(:,1));
                    Sens_Fact(:,1) = temp(1:3) / (OptProp_base(which_layer(1)) * p);
                    temp = DTOF_DelMom(temp_mom_base, temp_mom(:,2));
                    Sens_Fact(:,2) = temp(1:3) / (OptProp_base(which_layer(2)) * p);
                else
                    % For a more realistic change in optical properties
                    temp_change = OptProp_base(which_layer(1)) + delMua_guess(1); % 2nd DTOF
                        [temp_R] = LMA_1_RearrangeInput(temp_change, temp_base1, n, rho, time_ns, [1 length(time_ns)], -1, 1);
                        temp_mom(:,1) = DTOF_CentralMom(time_ns(cut_ind_mom(1):cut_ind_mom(2)),temp_R(cut_ind_mom(1):cut_ind_mom(2)));

                    temp_change = OptProp_base(which_layer(2)) + delMua_guess(2); % 3rd DTOF
                        [temp_R] = LMA_1_RearrangeInput(temp_change, temp_base2, n, rho, time_ns, [1 length(time_ns)], -1, 1);
                        temp_mom(:,2) = DTOF_CentralMom(time_ns(cut_ind_mom(1):cut_ind_mom(2)),temp_R(cut_ind_mom(1):cut_ind_mom(2)));

                    temp = DTOF_DelMom(temp_mom_base, temp_mom(:,1));
                    Sens_Fact(:,1) = temp(1:3) / delMua_guess(1);
                    temp = DTOF_DelMom(temp_mom_base, temp_mom(:,2));
                    Sens_Fact(:,2) = temp(1:3) / delMua_guess(2);
                end
                if ~isnan(Sens_Fact(1,1))
                    delMua_guess = (transpose(Sens_Fact(which_mom_use,:)) * (K(which_mom_use,which_mom_use)\Sens_Fact(which_mom_use,:)))...
                                \ (transpose(Sens_Fact(which_mom_use,:)) * (K(which_mom_use,which_mom_use)\delMom(which_mom_use)));
                end
            end
            if ~isnan(delMua_guess(1)) && max(delMua_guess)<0.02 && min(delMua_guess)>-0.02
                initial_guess(1:2) = delMua_guess(1:2);
            end

        elseif length(which_layer) == 1 % Finding changes in 1 layer
            p = 2 / 100;
            temp_base1 = OptProp_base; % Baseline
            temp_base1(which_layer) = -7;
            temp_base1(delOptProp_all==-8) = -8;
            temp_base1(7:8) = L_temp;
            Sens_Fact = zeros(3,1);
            for j = 1:2 % Run two times, first time with +-"p" percent change and second time with more realistic changes
                if j == 1 
                    temp_change = OptProp_base(which_layer); % 1st DTOF (baseline)
                        [temp_R] = LMA_1_RearrangeInput(temp_change, temp_base1, n, rho, time_ns, [1 length(time_ns)], -1, 1);
                        temp_mom_base = DTOF_CentralMom(time_ns(cut_ind_mom(1):cut_ind_mom(2)),temp_R(cut_ind_mom(1):cut_ind_mom(2)));

                    temp_change = OptProp_base(which_layer) * (1+p); % 2nd DTOF
                        [temp_R] = LMA_1_RearrangeInput(temp_change, temp_base1, n, rho, time_ns, [1 length(time_ns)], -1, 1);
                        temp_mom = DTOF_CentralMom(time_ns(cut_ind_mom(1):cut_ind_mom(2)),temp_R(cut_ind_mom(1):cut_ind_mom(2)));

                    temp = DTOF_DelMom(temp_mom_base, temp_mom);
                    Sens_Fact(:) = temp(1:3) / (OptProp_base(which_layer) * p);
                else
                    % For a more realistic change in optical properties
                    temp_change = OptProp_base(which_layer) + delMua_guess; % 2nd DTOF
                        [temp_R] = LMA_1_RearrangeInput(temp_change, temp_base1, n, rho, time_ns, [1 length(time_ns)], -1, 1);
                        temp_mom = DTOF_CentralMom(time_ns(cut_ind_mom(1):cut_ind_mom(2)),temp_R(cut_ind_mom(1):cut_ind_mom(2)));
        
                    temp = DTOF_DelMom(temp_mom_base, temp_mom);
                    Sens_Fact(:) = temp(1:3) / delMua_guess(1);
                end
                if ~isnan(Sens_Fact(1))
                    delMua_guess = (transpose(Sens_Fact(which_mom_use,:)) * (K(which_mom_use,which_mom_use)\Sens_Fact(which_mom_use,:)))...
                                \ (transpose(Sens_Fact(which_mom_use,:)) * (K(which_mom_use,which_mom_use)\delMom(which_mom_use)));
                end
            end
            if ~isnan(delMua_guess(1)) && max(delMua_guess)<0.02 && min(delMua_guess)>-0.02
                initial_guess(1) = delMua_guess;
            end
        end
    end
end
if isequal(initial_guess, -10) % If the initial guess hasn't been set
    initial_guess = typical_X; % Set initial guess to the typical X

    if do_disp(1) == 1
        disp([char(datetime('now','Format','HH:mm:ss')) '  Initial guess set to typical values (LMA 2):']); disp(initial_guess)
    end
else
    if do_disp(1) == 1
        disp([char(datetime('now','Format','HH:mm:ss')) '  Initial guess (LMA 2):']); disp(initial_guess)
    end
end


%% Prepare changes in moments and options

% Prepare IRF:
if ~isequal(irf_shifted, -1)
    irf_shifted(irf_shifted <= 0) = 1e-20; % Handle the noise, get rid of non physical values, prevent numerical problems
%     irf_shifted = irf_shifted / sum(irf_shifted); % Normalize by photon counts. Probably not needed
end

% Weight of moments, e.g. add more weight to m1 and V, and reduce the weight to N
for j = 1:3
    delMom(j) = delMom(j) / weight_noise(j);
end
% delMom = delMom * 10^3; % <This is not needed, but helps in displaying values

delMom = delMom(which_mom_use);

% Least square fitting options
options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective',...
                     'Display','off',... % 'Display','iter-detailed',...
                     'FunctionTolerance', 1e-9,...
                     'StepTolerance',1e-9,...
                     'TypicalX',typical_X,...
                     'MaxFunctionEvaluations', 200,...
                     'MaxIterations', 100);
if do_disp(2) == 1
    options.Display = 'iter-detailed'; % To output in the Command Window each Iteration of LMA
end


%% Call the LMA function

% Function handle:
fun_diff = @(delOptProp_X, time_ns)LMA_2_RearrangeInput(delOptProp_X, delOptProp_all, OptProp_base, n, rho, ...
                                                            time_ns, cut_ind_mom, irf_shifted, which_mom_use, weight_noise);

try
    [opt_prop_found, ~, Resid, ~, ~] = lsqcurvefit(fun_diff, initial_guess, time_ns, delMom, lbound, rbound, options);

    % Output result:
    if sum( opt_prop_found >= rbound*0.995 ) > 0
        disp([char(datetime('now','Format','HH:mm:ss')) '  WARNING: Reached ''rbound'' (LMA 2):']); disp(opt_prop_found)
    elseif sum( opt_prop_found <= lbound*0.995 ) > 0
        disp([char(datetime('now','Format','HH:mm:ss')) '  WARNING: Reached ''lbound'' (LMA 2):']); disp(opt_prop_found)
    else
        if do_disp(3) == 1
            disp([char(datetime('now','Format','HH:mm:ss')) '  Determined values (LMA 2):']); disp(opt_prop_found)
        end
    end

    % Manually calculate residuals to show what is being minimized by LMA. "Resid_manually" equals "Resid" outputted by the above LMA function
    do_this = 0;
    if do_this == 1
        [delMom_Found] = LMA_2_RearrangeInput(opt_prop_found, delOptProp_all, OptProp_base, n, rho, time_ns, cut_ind_mom, irf_shifted, which_mom_use, weight_noise);
        Resid_manually = delMom_Found - delMom;
        disp('Residuals from LMA and manually calculated:'); disp([Resid Resid_manually])
    end
    
catch ME
   opt_prop_found = nan(size(initial_guess));
   Resid = nan(1);
   disp([char(datetime('now','Format','HH:mm:ss')) '  WARNING: LMA function terminated, output result is set to NaN.' newline() 'Error: ' ME.message])
end

end