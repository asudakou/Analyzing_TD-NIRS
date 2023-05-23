% Written by Aleh Sudakou (https://www.researchgate.net/profile/Aleh-Sudakou)
% Nałęcz Institute of Biocybernetics and Biomedical Engineering, Polish Academy of Sciences 
% This is one of the codes shared on  https://github.com/asudakou/Analyzing_TD-NIRS
% Last updated: 20 May 2023


%% Description:
% This function calculates changes in moments for two sets of moments


function [delMom] = DTOF_DelMom(Mom_base, Mom_after)

delMom = zeros(4,1);
delMom(1) = -log( Mom_after(1) / Mom_base(1));  % 0th moment (Attenuation)
delMom(2) = Mom_after(2) - Mom_base(2);         % 1st moment (Mean time of flight)
delMom(3) = Mom_after(3) - Mom_base(3);         % 2nd central moment (Variance)
delMom(4) = Mom_after(4) - Mom_base(4);         % 3rd central moment (m3c)

end