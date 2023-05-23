% Written by Aleh Sudakou (https://www.researchgate.net/profile/Aleh-Sudakou)
% Nałęcz Institute of Biocybernetics and Biomedical Engineering, Polish Academy of Sciences 
% This is one of the codes shared on  https://github.com/asudakou/Analyzing_TD-NIRS
% Last updated: 20 May 2023


%% Description:
% This function filters DTOFs by subtracting background and cutting on two sides of the maximum.
% Output: "cut_ind" are the indices for which the photon count remains nonzero


function [ dtof, cut_ind ] = DTOF_filter( dtof, dt_ns, bg_reg, bg_how, cut_lim, cut_how, varargin) % varargin is for plotting


do_this = 0; % This must always be 0
if do_this == 1
    %% Run this section to see Example (Ctrl + Enter)
    OptProp_all_example = [0.01 -1 -1 1 -1 -1 0 0]; % Homogeneous optical properties, Mua=0.01 and Musp=1
    [R_example, time_ns_example] = DTOF_generate_Liemert(OptProp_all_example, 1.4, 30, -1); % Generate DTOF
    R_example(R_example<0) = 0;
    bg_reg_example = [3 4]; bg_how_example = 2;   % For subtracting the mean photon count between 3 and 4 ns (background)
    cut_lim_example = [0.01 0.01]; cut_how_example = 1;  % For cutting below 1% on both sides of the maximum, after subtracting background
    plot_option_example = {5, 'norm', 'raw', 'blue'}; % Plotting options
    figure(5); clf
    [~, ~] = DTOF_filter( R_example, time_ns_example, bg_reg_example, bg_how_example, cut_lim_example, cut_how_example, plot_option_example); % Call the current function
    clear OptProp_all_example R_example time_ns_example bg_reg_example bg_how_example cut_lim_example cut_how_example plot_option_example;
end


% INPUTS:
% * dt_ns  :  the width of one time-channel in ns or a list of all time-channels

% * bg_reg  :  the region for background subtraction

% * bg_how == 0  :  not subtracting background, ignoring "bg_reg"
%   bg_how == 1  :  subtract constant value equal to "bg_reg"
%   bg_how == 2  :  Subtract the mean count rate between nanoseconds "bg_reg(1)" and "bg_reg(2)"
%   bg_how == 3  :  Subtract the mean count rate between indices "bg_reg(1)" and "bg_reg(2)"

% * cut_lim  :  the region where to cut dtof

% * cut_how == 0  :  not cut  -  ignoring "cut_lim"
%   cut_how == 1  :  fraction -  "cut_lim = [0.01 0.01]" will set to 0 below 1% of the maximum on both sides
%   cut_how == 2  :  ns       -  "cut_lim = [1.5 3.5]"   will set to 0 before 1.5 and after 3.5 ns
%   cut_how == 3  :  indices  -  "cut_lim = [100 200]"   will set to zero before 100 and after 200 (exclusive)

% varargin is for plotting DTOF
% * plot_option = varargin{1};   For example  "plot_option = {5, 'blue'}"
%   plot_option{1}        ==  X       :  figure number (X is any number)
%   plot_option{after 1}  ==  'norm'  :  Plot normalized
%   plot_option{after 1}  ==  'raw'   :  Also plot unfiltered DTOF
%   plot_option{after 1}  ==  'no_bg' :  Not plot lines for background region
%   plot_option{after 1}  ==  'st', '.', '-', '--', 'o', 'x', '+'
%   plot_option{after 1}  ==  'red'   :  Or any other color 
% Example:
% plot_option = {1, 'norm', 'raw', 'no_bg', '.', 'blue'}
% plot_option = {5} % Just to plot on figure 5 with default colors


if isempty(varargin)
    fig_num = 0; fig_raw = 0; % Not plot
else
    plot_option = varargin{1};
    if length(plot_option) > 6 || length(varargin) > 1
        error('Wrong number of input arguments')
    end

    if isnumeric(plot_option) == 1
        fig_num = plot_option;
    else
        fig_num = plot_option{1};
    end
    fig_norm = 0;
    fig_raw = 0;
    no_bg = 0;
    st = '.'; % Default line style
    fig_col = 0;
    
    for j = 1:5
        if length(plot_option) > j % Additional inputs
            temp = plot_option{j+1};
            if isequal(temp,'norm')
                fig_norm = 1;
            elseif isequal(temp,'raw')
                fig_raw = 1;
            elseif isequal(temp,'no_bg')
                no_bg = 1;
            elseif isequal(temp,'.') || isequal(temp,'-') || isequal(temp,'--')...
                    || isequal(temp,'o') || isequal(temp,'x') || isequal(temp,'+')
                st = temp;
            elseif ~isequal(temp,-1) && ~isequal(temp,0)
                fig_col = temp;
            end
        end
    end
end

if cut_how == 1 % fraction
    if cut_lim(1) >= 1
        error('Input as a fraction. Seems it was inputted as percentage')
%         cut1 = cut1 / 100; % convert to fraction
%         cut2 = cut2 / 100;
    end
    if cut_lim(1) == 0
        cut_lim(1) = 1e-20;
    end
    if cut_lim(2) == 0
        cut_lim(2) = 1e-20;
    end
end

if dt_ns(end) < 1e-6
    error('Time channels should be in units of ns and not s')
%     time_ns = time_ns * 1e9; % convert to ns
end

if fig_raw > 0 % If want to plot also unfiltered DTOF
    dtof_raw = dtof;
end


%% Filter background noise
if bg_how == 1 % SUBTRACT CONSTANT VALUE
    dtof = dtof - bg_reg;

elseif bg_how > 1
    if bg_how == 2 % SUBTRACT MEAN BETWEEN NANOSECONDS
        if length(dt_ns) == 1 % if dt_ns is width of time-channel
            bg_ind = round ( bg_reg / dt_ns );
        else % if dt_ns are all time-channels
            [~, bg_ind] = min(abs(dt_ns - bg_reg(1)));
            [~, bg_ind(2)] = min(abs(dt_ns - bg_reg(2)));
        end

    elseif bg_how == 3 % SUBTRACT MEAN BETWEEN INDICES
        bg_ind = bg_reg;
    end

    for j = 1:2 % Check if found indices are below 0 or above length of dtof
        if bg_ind(j) <= 0
            bg_ind(j) = 1;
        end
        if bg_ind(j) > length(dtof)
            bg_ind(J) = length(dtof);
        end
    end

    dtof = dtof - mean(dtof(bg_ind(1):bg_ind(2))); % Subtract
end

if bg_how > 0
    dtof(dtof<0) = 0; % remove negative data - assuming this data is hidden in background noise
end


%% Cut at t_lower and t_upper
cut_ind = [-1 -1];
if cut_how > 0
    [max_val, max_ind] = max(dtof);
    
    if max_val > 0 % If it is not > 0, then it is probably an empty dtof and then no need to cut
        
        if cut_how == 1 % FRACTION
            for ind = max_ind:-1:1 % Lower limit:
                if dtof(ind) < cut_lim(1) * max_val
                    cut_ind(1) = ind; % Cut everything from this index and lower
                    break;
                end
            end

            for ind = max_ind:length(dtof) % Repeat for Upper limit
                if dtof(ind) < cut_lim(2) * max_val
                    cut_ind(2) = ind;
                    break;
                end
            end

        elseif cut_how == 2 % NANOSECONDS
            if length(dt_ns) == 1 % if dt_ns is width of time-channel
                cut_ind(1) = floor ( cut_lim(1) / dt_ns ) - 1; % Cut everything from this index and lower
            else % if dt_ns are all time-channels
                [~, cut_ind(1)] = min(abs(dt_ns - cut_lim(1)));
            end

            if length(dt_ns) == 1 % Repeat for Upper limit
                cut_ind(1) = floor ( cut_lim(1) / dt_ns ) - 1;
            else
                [~, cut_ind(1)] = min(abs(dt_ns - cut_lim(1)));
            end

        elseif cut_how == 3 % INDECISE
            cut_ind(1) = cut_lim(1) - 1; % Cut everything from this index and lower
            cut_ind(2) = cut_lim(2) + 1; % Cut everything from this index and higher
        end
        
        % Cut
        if cut_ind(1) >= 1 && cut_ind(1) <= length(dtof)
            dtof(1:cut_ind(1)) = 0;
        end
        if cut_ind(2) >= 1 && cut_ind(2) <= length(dtof)
            dtof(cut_ind(2):end) = 0;
        end

        cut_ind(1) = cut_ind(1)+1; % For outputting which indices remain for the DTOF
        cut_ind(2) = cut_ind(2)-1;
    end
end


%% Plot result:
if fig_num > 0
    figure(fig_num); 
    % clf

    if length(dt_ns) == 1 % if dt_ns is width of time-channel. Otherwise it is a list of time channels
        dt_ns = dt_ns:dt_ns:dt_ns*length(dtof);
    end
    
    if fig_norm == 1
        if fig_raw > 0
            semilogy(dt_ns, dtof_raw / max(dtof_raw), st, 'Color','black')
            hold on
        end
        if ~isequal(fig_col,0)
            semilogy(dt_ns, dtof / max(dtof), st, 'Color',fig_col) % Specified color
        else
            semilogy(dt_ns, dtof / max(dtof), st) % Random color
        end
        hold on
    else
        if fig_raw > 0
            semilogy(dt_ns, dtof_raw, st, 'Color','black')
            hold on
        end
        if ~isequal(fig_col,0)
            semilogy(dt_ns, dtof, st, 'Color',fig_col) % Specified color
        else
            semilogy(dt_ns, dtof, st) % Default color
        end
        hold on
    end

    xlim([0 dt_ns(end)])
    set(gca,'YScale','log')
    grid on

    if bg_how > 1 && no_bg == 0
        plot(dt_ns(bg_ind(1))*[1 1],ylim,'--','Color','red','linewidth',1.2,'handlevisibility','off')
        plot(dt_ns(bg_ind(2))*[1 1],ylim,'--','Color','red','linewidth',1.2,'handlevisibility','off')
    end
%     if cut_how > 0
%         plot(dt_ns(cut_ind(1))*[1 1],ylim,'--','Color','red','linewidth',1.5,'handlevisibility','off')
%         plot(dt_ns(cut_ind(2))*[1 1],ylim,'--','Color','red','linewidth',1.5,'handlevisibility','off')
%     end
end
end