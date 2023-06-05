% Written by Aleh Sudakou (https://www.researchgate.net/profile/Aleh-Sudakou)
% Nałęcz Institute of Biocybernetics and Biomedical Engineering, Polish Academy of Sciences 
% This is one of the codes shared on  https://github.com/asudakou/Analyzing_TD-NIRS
% Last updated: 20 May 2023


%% Description:
% This function implements what is known as the "multi-layered curve-fitting method".
% Using LMA algorithm (MATLAB's built-in function), search for the absolute
% optical properties that best correspond to the measured DTOF. 
%
% The initial guess is calculated using the method based on moments, using the function "Mom_to_OptProp.m"
% 
% Required files: 
% Matlab Optimization Toolbox, for using in-built function "lsqcurvefit"
%
% Three functions: "LMA_1_RearrangeInput.m"
%                  "DTOF_generate_Liemert.m"
%                  "DTOF_convolve.m" (If you have DTOF)
% 
% Inputs:
% * OptProp_all  :  8 knowns and unknowns [Mua1 Mua2 Mua3 Musp1 Musp2 Musp3 L1 L2]     [Units mm-1] and [mm for L1 and L2]
%   Value -7  means this parameter is unknown
%   Value -8  means to take the same value as for the above layer (cannot be assigned to the superficial layer or L)
%   Any other value means this is the value of the parameter.
%   (For a homogeneous medium, set to -8 both second and third layers, values of L1 and L2 won't matter)
%   (For a two-layered medium, set to -8 the third layer; value of L2 won't matter)
%   Examples of determining different parameters with this function:
%   OptProp_all = [-7 -8 -8 -7 -8 -8 0  0]      % Mua and Musp,  Homogeneous medium
%   OptProp_all = [-7 -7 -8  1 -8 -8 10 0]      % Mua1 and Mua2, Two-layered medium 
%   OptProp_all = [-7  0.01 -7  0  0  0 10 10]  % Mua1 and Mua3, Three-layered medium 
%   OptProp_all = [-7  0.01 -8  1 -8 -8 10 10]  % Mua1 and Musp(all), Three-layered medium, two unknowns
%   OptProp_all = [0.01 -7  -8  1  -8  -8 -7 0] % Mua2 and L,    Two-layered medium
% 
% * n         :  refractive index, e.g. 1.4 or [1.4 1.4 1.4]
% * rho       :  source-detector distance in mm
% * time_ns   :  time channels of measured DTOF
% * cut_ind_fit :  indices of time channels used for the fitting
% 
% * irf_shifted : variable containing IRF (having been moved by offset). Can be set to -1 if not used
%                 Note, irf_shifted should correspond to time channels of DTOF for a proper convolution
%                 Therefore, to account for offset, use "interp1" built-in function
% 
% * varargin            :  possible optional inputs
%   varargin.weight_exp :  DTOF = DTOF.^weight_exp. Can accentuate the weight of later photons, e.g. 1/2 for the square root of DTOF
%   varargin.init_guess :  Use this as initial guess or if set to 'no' then use typicalX as the initial guess
%   varargin.do_disp    :  [1 1 1 1] to display or plot: 1) Initial guess, 2) "iter-detailed" LMA, 3) Determined values, and 4) Plot figures


function [opt_prop_found, Resid] = LMA_1_FittingDTOF(OptProp_all, n, rho, dtof, time_ns, cut_ind_fit, irf_shifted, varargin)


do_this = 0; % This must always be 0
if do_this == 1
    %% Run this section to see Example (Ctrl + Enter)
    which_example = 6; % Set to 1, 2, 3, 4, 5, or 6
    
    % 1) Generate DTOFs
    if which_example == 1
        OptProp_base_example = [0.01 0.01 0.01 1 1 1 5 5]; % Mua=0.01 and Musp=1, Homogeneous
    else
        OptProp_base_example = [0.01 0.015 0.015 1 1 1 5 5]; % Mua1=0.01 and Mua2=0.015
    end
    [R1_example, time_ns_example] = DTOF_generate_Liemert(OptProp_base_example, 1.4, 30, -1);

    % 2) Choose region for fitting
    [ ~, cut_ind_fit_example] = DTOF_filter( R1_example, time_ns_example, 0, 0, [0.01 0.01], 1);

    % 3) Call the current function to determine Optical properties
    options_LMA1_example = [];
    options_LMA1_example.weight = [1 1 1]; % Since we have an inverse-crime here, weight makes no difference

    options_LMA1_example.do_disp = [1 0 1 1]; % What information to display/plot
                                   % 1) Initial guess, 2) 'iter-detailed' LMA, 3) Determined values, and 4) Plot figures
    if which_example == 1
        OptProp_all_example = [-7 -8 -8 -7 -8 -8 0 0]; % Find Mua and Musp, homogeneous

    elseif which_example == 2
        OptProp_all_example = OptProp_base_example;
        OptProp_all_example(4) = -7; % Find only Musp3
        
    elseif which_example == 3
        OptProp_all_example = [-7 -7 -8 1 -8 -8 5 0]; % Find Mua1 and Mua2, two-layered
        options_LMA1_example.no_init_guess = 1;

    elseif which_example == 4
        OptProp_all_example = [-7 -7 -8 1 -8 -8 -7 0]; % Find Mua1, Mua2, and L, two-layered

    elseif which_example == 5
        OptProp_all_example = [-7 -7 -8 -7 -8 -8 5 0]; % Find Mua1, Mua2, and Musp(assuming Musp is the same in all 3 layers), two-layered

    elseif which_example == 6
        OptProp_all_example = [0.01 0.015 -7 -7 -8 -8 5 5]; % Find Mua3 and Musp(assuming Musp is the same in all 3 layers), three-layered

    end
    
    [~, ~] = LMA_1_FittingDTOF(OptProp_all_example, 1.4, 30, R1_example, time_ns_example, cut_ind_fit_example, -1, options_LMA1_example);
    clear which_example OptProp_base_example R1_example time_ns_example cut_ind_fit_example options_LMA1_example OptProp_all_example;
    return
end


% INPUTS:
if length(OptProp_all) ~= 8
    error('variable ""opt_prop_input"" must be length 8')
end
if rho < 5
    error('Check source-detector distance (rho)')
end
if n(1) < 1.2 || n(1) > 1.6
    error('Check refractive index (n), usually it should be between 1.2 and 1.6 for tissue, otherwise comment out this error message')
end
if size(dtof,2)>1
    dtof = dtof';
    if size(dtof,2)>1
        error('Dont input more than 1 DTOF. Otherwise, there may start to exist many local minima')
    end
end
if cut_ind_fit(1) < 0 || cut_ind_fit(2) < 0 || sum(dtof) == 0
    opt_prop_found = nan(sum(OptProp_all==-7),1);
    Resid = nan(1);
    disp([char(datetime('now','Format','HH:mm:ss')) '  WARNING: empty DTOF or wrong cut indices, output result is set to NaN.'])
    return
end

% Default values:
typicalX_all = [0.01 0.01 0.01 1 1 1 10 10]; % [Mua1 Mua2 Mua3 Musp1 Musp2 Musp3 L1 L2]
lbound_all = [0 0 0 0 0 0 0 0]; % Lower boundary values
rbound_all = [0.1 0.1 0.1 3 3 3 35 35]; % Upper

typical_X = typicalX_all(OptProp_all==-7); % Select only for the unknown properties
lbound = lbound_all(OptProp_all==-7);
rbound = rbound_all(OptProp_all==-7);

weight_exp = 1; % To accentuate the relative weight of later photons

do_disp = [0 0 0 0]; % 1) Initial guess, 2) 'iter-detailed' LMA, 3) Determined values, and 4) Plot figures

initial_guess = -10; % Will calculate the initial guess
% initial_guess = typical_X; % Uncomment this to not calculate the initial guess

% Optional inputs:
if ~isempty(varargin)
    if length(varargin) == 1
        if isfield(varargin{1},'do_disp')
            do_disp = varargin{1}.do_disp;
        end
        if isfield(varargin{1},'weight_exp')
            weight_exp = varargin{1}.weight_exp;
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
    else
        error('Wrong number of input arguments')
    end
end


% Calculate initial guess if it wasn't specified.
% To avoid this step, set "do_this" to 0 or uncomment the above line "initial_guess = typical_X;"
% Finding a good initial guess speeds up calculations, as can be easily verified by 
% calculating with and without this code or by comparing the initial guess and the 
% determined values (see "do_disp" for displaying the initial guess and the determined values)
do_this = 1;
% do_this = 0;
if isequal(initial_guess, -10) && do_this == 1
    [ temp, ~ ] = DTOF_filter( dtof, time_ns, 0, 0, [0.01 0.01], 1); % Use 1 percent on both sides to calc. moments
    moments_dtof = DTOF_CentralMom(time_ns, temp);

    if isequal(length(irf_shifted), 1)
        moments_irf = [0 0 0];
    else
        [ temp, ~ ] = DTOF_filter( irf_shifted, time_ns, 0, 0, [0.01 0.01], 1);
        moments_irf = DTOF_CentralMom(time_ns, temp);
    end

    temp = moments_dtof(2:3) - moments_irf(2:3);
    [mua1_ini, musp1_ini] = Mom_to_OptProp(temp(1)/10^9, temp(2)/10^18, 0, 0, rho, n(1), 0);

    if ~isnan(mua1_ini) && mua1_ini > 0.001 && mua1_ini < 0.04 && musp1_ini > 0.3 && musp1_ini < 3
        temp = [mua1_ini mua1_ini mua1_ini musp1_ini musp1_ini musp1_ini typicalX_all(7:8)];
        initial_guess = temp(OptProp_all==-7);
    end
end
if isequal(initial_guess, -10) % If the initial guess hasn't been set
    initial_guess = typical_X; % Set initial guess to the typical X

    if do_disp(1) == 1
        disp([char(datetime('now','Format','HH:mm:ss')) '  Initial guess set to typical values (LMA 1):']); disp(initial_guess)
    end
else
    if do_disp(1) == 1
        disp([char(datetime('now','Format','HH:mm:ss')) '  Initial guess (LMA 1):']); disp(initial_guess)
    end
end


%% Prepare irfs and dtofs
% irf(irf <= 0) = 1e-20; % Handle the noise, get rid of non physical values, prevent numerical problems
% irf = irf / sum(irf); % Normalize by photon counts. Probably not needed

dtof_cut = dtof(cut_ind_fit(1):cut_ind_fit(2));
dtof_cut = dtof_cut / sum(dtof_cut) * 10^6; % normalize so that dtof has Ntot = 10^6. The same will be done for theoretically calculated DTOF
% dtof_4_fit(dtof_4_fit <= 0) = 1e-20;

dtof_cut = dtof_cut.^weight_exp;

% Least square fitting options
options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective',...
                     'Display','off',... % 'Display','iter-detailed',...
                     'FunctionTolerance', 1e-9,...
                     'StepTolerance',1e-9,...
                     'TypicalX',typical_X,... 
                     'MaxFunctionEvaluations', 300,...
                     'MaxIterations',100);
if do_disp(2) == 1
    options.Display = 'iter-detailed'; % Output in the Command Window each Iteration of LMA
end


%% Call the LMA function

% Function handle:
fun_diff = @(OptProp_X, time_ns_cut)LMA_1_RearrangeInput(OptProp_X, OptProp_all, n, rho, time_ns, cut_ind_fit, irf_shifted, weight_exp);

try
    [opt_prop_found, ~, Resid, ~, ~] = lsqcurvefit(fun_diff, initial_guess, time_ns, dtof_cut, lbound, rbound, options);
    
    % Output result:
    if sum( opt_prop_found >= rbound*0.995 ) > 0
        disp([char(datetime('now','Format','HH:mm:ss')) '  WARNING: Reached ''rbound'' (LMA 1):']); disp(opt_prop_found)
    elseif sum( opt_prop_found <= lbound*0.995 ) > 0
        disp([char(datetime('now','Format','HH:mm:ss')) '  WARNING: Reached ''lbound'' (LMA 1):']); disp(opt_prop_found)
    else
        if do_disp(3) == 1
            disp([char(datetime('now','Format','HH:mm:ss')) '  Determined values (LMA 1):']); disp(opt_prop_found)
        end
    end

    % Manually calculate residuals to show what is being minimized by LMA. "Resid_manually" equals "Resid" outputted by the above LMA function
    do_this = 0;
    if do_this == 1
        R = LMA_1_RearrangeInput(opt_prop_found, OptProp_all, n, rho, time_ns, cut_ind_fit, irf_shifted, weight_exp);
        Resid_manually = R - dtof_cut;
        disp('Residuals outputted by LMA and manually calculated:'); disp([Resid Resid_manually])
    end

catch ME
   opt_prop_found = nan(size(initial_guess));
   Resid = nan(1);
   disp([char(datetime('now','Format','HH:mm:ss')) '  WARNING: LMA function terminated, output result is set to NaN.' newline() 'Error: ' ME.message])
end


%% Plot the fit if requested 
if do_disp(4) == 1
    fig_num = 101;

    figure(fig_num); clf
    
    % Measured DTOF:
%     semilogy(time_ns, dtof,'Color','blue') % Simple way of plotting
    [ ~, ~ ] = DTOF_filter( dtof, time_ns, 0, 0, 0, 0, {fig_num,'o','blue','norm'}); % Better way of plotting
    hold on
    
    % IRF:
    if length(irf_shifted)>1
        [ ~, ~ ] = DTOF_filter( irf_shifted, time_ns, 0, 0, 0, 0, {fig_num,'x','green','norm'});
    else
        plot([-1 -1],[-1 -1], 'Color','green') % For the legend
    end
    
    % DTOF for determined optical properties:
    R = LMA_1_RearrangeInput(opt_prop_found, OptProp_all, n, rho, time_ns, [1 length(time_ns)], irf_shifted, weight_exp);
    [ ~, ~ ] = DTOF_filter( R, time_ns, 0, 0, 0, 0, {fig_num,'red','norm'});

    xlabel('Time / ns','Interpreter','latex','FontSize',20)
    ylabel('Normalized counts','Interpreter','latex','FontSize',20)
    legend({'DTOF - measured','IRF - measured','DTOF - fitted'},'Interpreter','latex','FontSize',20,'Location','NorthEast')
    set(gca,'FontSize',20)
end
end