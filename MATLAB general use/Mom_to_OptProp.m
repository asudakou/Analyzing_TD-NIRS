% Written by Aleh Sudakou (https://www.researchgate.net/profile/Aleh-Sudakou)
% Nałęcz Institute of Biocybernetics and Biomedical Engineering, Polish Academy of Sciences 
% This is one of the codes shared on  https://github.com/asudakou/Analyzing_TD-NIRS
% Last updated: 20 May 2023


%% Description:
% Calculate optical properties based on statistical moments of DTOFs.
% 
% PUBLICATION:
% A. Liebert, H. Wabnitz, D. Grosenick, M. Moller, R. Macdonald, and H. Rinneberg, 
% "Evaluation of optical properties of highly scattering media by moments of 
% distributions of times of flight of photons," Appl Opt 42, 5785-5792 (2003).
% DOI:  https://doi.org/10.1364/AO.42.005785


function [mua, musp] = Mom_to_OptProp(dtof_m1, dtof_V, irf_m1, irf_V, r, n, offset)


do_this = 0; % This must always be 0
if do_this == 1
    %% Run this section to see Example (Ctrl + Enter)
    OptProp_all_example = [0.01 -1 -1 1 -1 -1 0 0]; % Homogeneous optical properties
    [R_example, time_ns_example] = DTOF_generate_Liemert(OptProp_all_example, 1.4, 30, -1); % Generate DTOF
    [mom_centr_example] = DTOF_CentralMom(time_ns_example / 10^9, R_example); % Calculate Moments
    [mua_example, musp_example] = Mom_to_OptProp(mom_centr_example(2), mom_centr_example(3), 0, 0, 30, 1.4, 0); % Call the current function to retrieve optical properties
    disp('Determined Mua and Musp are:'); disp([mua_example musp_example])
    clear OptProp_all_example R_example time_ns_example mom_centr_example mua_example musp_example;
end


% INPUTS:
% dtof_m1; % == mean time of flight of DTOF (first statistical moment). Units: seconds
% irf_m1;  % == for IRF
% dtof_V;  % == variance of DTOF (second central statistical moment). Units: seconds^2
% irf_V;   % == for IRF
% r;       % == distance between source and detector fibers. Units: mm
% offset;  % == Time offset due to the distance between fibers during IRF measurement. Units: seconds
% For example: offset = 0.06/physconst('Lightspeed'); % Time offset of IRF due the 6 cm distance between fibers when measuring IRF
% n;       % == refractive index of medium


% OUTPUT : 
% mua     == absorption coefficient. Units mm^-1
% musp    == reduced scattering coefficient. Units mm^-1


%% Check units:
if dtof_m1 > 1e-6 || irf_m1 > 1e-6 % < must have units: seconds
%     dtof_m1 = dtof_m1 / 1e9;
    error('Please confirm that ''dtof_m1'' and ''irf_m1'' have units seconds and not nanoseconds')
end

if dtof_V > 1e-6 || irf_V > 1e-6 % < must have units: seconds
%     dtof_V = dtof_V / 1e18;
    error('Please confirm that ''dtof_V'' and ''irf_V'' have units seconds^2 and not nanoseconds^2')
end

if r < 10 % < must have units: mm
%     r = r * 10;
    error('Please confirm that ''r'' has units mm and not cm. If so, then comment out this line')
end

if offset > 1e-6 % < must have units: seconds
%     offset = offset / 1e9;
    error('Please confirm that ''offset'' has units seconds and not nanoseconds')
end


%% Apply formula:
C = (physconst('LightSpeed') / n) * 10^3; % units: mm/s

m1 = dtof_m1 - (irf_m1 - offset);

V = dtof_V - irf_V;

mua = m1.^3./(2.*C.*V.*(m1.^2 + V));

musp = 2.*C.*m1.*(m1.^2 + V)./(3.*r.*r.*V);


end