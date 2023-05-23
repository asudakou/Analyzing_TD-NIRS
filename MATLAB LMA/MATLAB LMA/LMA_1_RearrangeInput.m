% Written by Aleh Sudakou (https://www.researchgate.net/profile/Aleh-Sudakou)
% Nałęcz Institute of Biocybernetics and Biomedical Engineering, Polish Academy of Sciences 
% This is one of the codes shared on  https://github.com/asudakou/Analyzing_TD-NIRS
% Last updated: 20 May 2023


%% Description:
% This function is an intermediate function between "LMA_2_FittingDelMom.m" and "DTOF_generate_Liemert.m"
% The purpose is to rearrange inputs for LMA, since we allow any optical property to be set as unknown. 
% For details, see "LMA_1_FittingDTOF.m", which has an example of calling the current function for manually calc. residuals


function [R] = LMA_1_RearrangeInput(OptProp_X, OptProp_all, n, rho, time_ns, cut_ind_fit, irf_shifted, weight_exp)


% INPUTS:
% * OptProp_X    :   Unknown optical property
% * OptProp_all  :   All optical properties, the unknown ones are marked as -7

if length(OptProp_all) ~= 8
    error('variable ""OptProp_all"" must be length 8')
end
if rho < 5
    error('rho (source-detector distance) should be in mm')
end


% Move properties from 'OptProp_X' into 'OptProp_all' and reorganize 'OptProp_all'
counter = 1;
if isequal(OptProp_all(1), -7) % Mua1
    OptProp_all(1) = OptProp_X(counter); counter = counter + 1;
end
if isequal(OptProp_all(2), -7) % Mua2
    OptProp_all(2) = OptProp_X(counter); counter = counter + 1;
elseif isequal(OptProp_all(2), -8)
    OptProp_all(2) = OptProp_all(1);
end
if isequal(OptProp_all(3), -7) % Mua3
    OptProp_all(3) = OptProp_X(counter); counter = counter + 1;
elseif isequal(OptProp_all(3), -8)
    OptProp_all(3) = OptProp_all(2);
end
if isequal(OptProp_all(4), -7) % Musp1
    OptProp_all(4) = OptProp_X(counter); counter = counter + 1;
end
if isequal(OptProp_all(5), -7) % Musp2
    OptProp_all(5) = OptProp_X(counter); counter = counter + 1;
elseif isequal(OptProp_all(5), -8)
    OptProp_all(5) = OptProp_all(4);
end
if isequal(OptProp_all(6), -7) % Musp3
    OptProp_all(6) = OptProp_X(counter); counter = counter + 1;
elseif isequal(OptProp_all(6), -8)
    OptProp_all(6) = OptProp_all(5);
end
if isequal(OptProp_all(7), -7) % L1
    OptProp_all(7) = OptProp_X(counter); counter = counter + 1;
end
if isequal(OptProp_all(8), -7) % L2
    OptProp_all(8) = OptProp_X(counter); % counter = counter + 1;
end


% % Display the guess of the optical properties, i.e. what is the input to this function
% temp_text = cell(length(OptProp_X));
% for j = 1:length(OptProp_X)
%     temp_text{j} = [' X_' num2str(j) ':  ' num2str(OptProp_X(j)) '   '];
% end
% disp([temp_text{1:end}])


%% Call the next function, as this is an intermediate function 

[R, ~] = DTOF_generate_Liemert(OptProp_all, n, rho, time_ns, irf_shifted);

R = R(cut_ind_fit(1):cut_ind_fit(2));

R = R / sum(R) * 10^6; % The same normalization is done in "LMA_1_FittingDTOF.m"

R = R.^weight_exp;


end