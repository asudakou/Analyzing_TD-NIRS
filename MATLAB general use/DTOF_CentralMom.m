% Written by Heidrun Wabnitz (https://www.researchgate.net/profile/Heidrun-Wabnitz)
% Physikalisch-Technische Bundesanstalt (PTB), Berlin
% This is one of the codes shared on  https://github.com/asudakou/Analyzing_TD-NIRS
% Last updated: 20 May 2023


%% Description:
% This function calculates central moments of DTOFs


function [mom_centr] = DTOF_CentralMom(time_ns, dtof)


do_this = 0; % This must always be 0
if do_this == 1
    %% Run this section to see Example (Ctrl + Enter)
    OptProp_all_example = [0.01 0 0 1 0 0 0 0]; % Homogeneous optical properties, Mua=0.01 and Musp=1
    [R_example, time_ns_example] = DTOF_generate_Liemert(OptProp_all_example, 1.4, 30, -1); % Generate DTOF
    [mom_centr_example] = DTOF_CentralMom(time_ns_example, R_example); % Call the current function to obtain moments
    disp('mom_centr is:'); disp(mom_centr_example(1)); disp(mom_centr_example(2:end))
    clear OptProp_all_example R_example time_ns_example mom_centr_example;
end


%% Main body
dtof = dtof(:);
time_ns = time_ns(:);

m0 = sum(dtof);
if m0 ~= 0
    % non-central moments:
    m1 = sum(dtof.*time_ns)/m0;
    m2 = sum(dtof.*time_ns.*time_ns)/m0;
    m3 = sum(dtof.*time_ns.*time_ns.*time_ns)/m0;
    m4 = sum(dtof.*time_ns.*time_ns.*time_ns.*time_ns)/m0;

    % central moments:
    mom_centr(1) = m0; % N (Total number of photons, 0th moment)
    mom_centr(2) = m1; % <t> (Mean time of flight, 1st moment)
    mom_centr(3) = m2 - m1*m1; % V (Variance, 2nd central moment)
    mom_centr(4) = m3 - 3*m2*m1 + 2*m1*m1*m1; % m3c (3rd central moment)
    mom_centr(5) = m4 - 4*m3*m1 + 6*m2*m1*m1 - 3*m1*m1*m1*m1; % m4c (4th central moment)
else
    mom_centr(1:5)=0;
end

end