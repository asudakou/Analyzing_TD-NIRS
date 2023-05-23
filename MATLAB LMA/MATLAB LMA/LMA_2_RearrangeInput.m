% Written by Aleh Sudakou (https://www.researchgate.net/profile/Aleh-Sudakou)
% Nałęcz Institute of Biocybernetics and Biomedical Engineering, Polish Academy of Sciences 
% This is one of the codes shared on  https://github.com/asudakou/Analyzing_TD-NIRS
% Last updated: 20 May 2023


%% Description:
% This function is an intermediate function between "LMA_2_FittingDelMom.m" and "DTOF_generate_Liemert.m"
% The purpose is to rearrange inputs for LMA, since we allow any optical property to be set as unknown. 
% For details, see "LMA_2_FittingDelMom.m", which has an example of calling the current function for manually calc. residuals


function [delMom] = LMA_2_RearrangeInput(delOptProp_X, delOptProp_all, OptProp_base, n, rho, time_ns, cut_ind_mom, irf_shifted, which_mom_use, weight_noise)


% INPUTS:
% * delOptProp_X    :  Unknown opitcal properties
% * delOptProp_all  :  All optical properties, the unknown ones are marked as -7
%                     [delMua1 delMua2 delMua3 delMusp1 delMusp2 delMusp3 L1 L2]

if length(delOptProp_all) ~= 8
    error('variable ""delOptProp_all"" must be length 8')
end
if length(OptProp_base) ~= 6
    error('variable ""OptProp_base"" must be length 6')
end
if rho < 5
    error('rho (source-detector distance) should be in mm')
end

% Move properties from 'delOptProp_X' into 'delOptProp_all' and reorganize 'delOptProp_all':
counter = 1;
if isequal(delOptProp_all(1), -7) % delMua1
    delOptProp_all(1) = delOptProp_X(counter); counter = counter + 1;
end
if isequal(delOptProp_all(2), -7) % delMua2
    delOptProp_all(2) = delOptProp_X(counter); counter = counter + 1;
elseif isequal(delOptProp_all(2), -8)
    delOptProp_all(2) = delOptProp_all(1);
end
if isequal(delOptProp_all(3), -7) % delMua3
    delOptProp_all(3) = delOptProp_X(counter); counter = counter + 1;
elseif isequal(delOptProp_all(3), -8)
    delOptProp_all(3) = delOptProp_all(2);
end
if isequal(delOptProp_all(4), -7) % delMusp1
    delOptProp_all(4) = delOptProp_X(counter); counter = counter + 1;
end
if isequal(delOptProp_all(5), -7) % delMusp2
    delOptProp_all(5) = delOptProp_X(counter); counter = counter + 1;
elseif isequal(delOptProp_all(5), -8)
    delOptProp_all(5) = delOptProp_all(4);
end
if isequal(delOptProp_all(6), -7) % delMusp3
    delOptProp_all(6) = delOptProp_X(counter); counter = counter + 1;
elseif isequal(delOptProp_all(6), -8)
    delOptProp_all(6) = delOptProp_all(5);
end
if isequal(delOptProp_all(7), -7) % L1
    delOptProp_all(7) = delOptProp_X(counter); counter = counter + 1;
end
if isequal(delOptProp_all(8), -7) % L2
    delOptProp_all(8) = delOptProp_X(counter); % counter = counter + 1;
end


% % Display the guess of the optical properties, i.e. what is the input to this function
% temp_text = cell(length(delOptProp_X));
% for j = 1:length(delOptProp_X)
%     temp_text{j} = [' X_' num2str(j) ':  ' num2str(delOptProp_X(j)) '   '];
% end
% disp([temp_text{1:end}])


%% Call the next function, as this is an intermediate function 

[delMom] = DelOptProp_to_DelMom(delOptProp_all, OptProp_base, n, rho, time_ns, cut_ind_mom, irf_shifted);

for j = 1:3
    delMom(j) = delMom(j) / weight_noise(j);
end
% delMom = delMom * 10^3; % <This is not needed, but helps in displaying values

delMom = delMom(which_mom_use); % Output only the used changes in moments


end