% Written by Aleh Sudakou (https://www.researchgate.net/profile/Aleh-Sudakou)
% Nałęcz Institute of Biocybernetics and Biomedical Engineering, Polish Academy of Sciences 
% This is one of the codes shared on  https://github.com/asudakou/Analyzing_TD-NIRS
% Last updated: 20 May 2023


%% Description:
% Returns the Nominal Mua values, which are calculated based on the measured
% weights of added ingridients. These values are used only for plotting by "Ink_plotting.m"
%
% Note, this script is good for seeing how to read the Extinction coefficients 
% and the absorption of water, obtained from literature (i.e. downloaded -
% see "Where_downloaded.txt"


function [mua_nom] = Ink_Nominal_mua(wavelengths) 

% 1) Correct Molar Absorption Coefficients for Oxy and Deoxy:
% wh_source = 1; % Gratzer   % < Used for Publication 2023
% wh_source = 2; % Moaveni
% wh_source = 3; % Takatani
% Ext_coeff = [];
% Ext_coeff(:,1) = wavelengths;
% Ext_coeff(:,2:end+5) = GetExtinctions(Ext_coeff(:,1), wh_source); % [HbO Hb H2O lipid aa3]   Molar Extinction Coefficient   [OD M-1 cm-1]
% Ext_coeff(:,2:3) = Ext_coeff(:,2:3) .* log(10);  % Molar Absorption Coefficient   [OD M-1 cm-1]   See GetExtinctions.m


% 2) Correct Water's Mua:
mua_Water_all = load('matcher94_nir_water_37.txt');
mua_Water_all(:,2) = mua_Water_all(:,2) * log(10); % Absorption spectrum of Water from Matcher [OD cm-1]  (Note, not per Mole)
mua_Water = interp1(mua_Water_all(:,1),mua_Water_all(:,2),wavelengths); % At 16 wavelengths


% 3) Second source of water's Mua:
do_this = 0; % Set to 1 to see a comparison
if do_this == 1
    mua_Water_2 = GetExtinctions(640:10:950, 2); % [HbO Hb H2O lipid aa3]    [OD cm-1] for Water
    
    figure(5); clf
    plot(640:10:950, mua_Water_2(:,3),'x','Color','red','LineWidth',2,'MarkerSize',10)
    hold on
    plot(mua_Water_all(:,1),mua_Water_all(:,2),'-','Color','blue')
    grid on
    xlim([600 1000])
    title('Comparison of water''s mua from two sources')
    set(gca,'FontSize',24)
end


V_ink = zeros(21,2);
V_tot = zeros(21,2);
mua_nom = zeros(21,14,2);

for j_Exp = 1:2
    % Measured weights in gramms (Inputted from Excel file)
    if j_Exp == 1
        V_ink(:,j_Exp) = [1.493889067	2.150426133	2.798727667	3.4460336	4.095199733	4.747300267	5.4038548	6.060942067	6.721496467	7.386138067	8.052133333	8.721228933	9.3921236	10.0665378	10.7421048	11.4255842	12.1118146	12.79861267	13.48607447	14.1723136	14.8683166];
        
        V_tot(:,j_Exp) = [3617.174022	3624.721693	3632.174686	3639.616234	3647.079166	3654.575834	3662.123705	3669.6777	3677.271555	3684.912396	3692.5688	3700.260846	3707.973574	3715.726763	3723.493205	3731.350609	3739.23964	3747.135196	3755.038383	3762.927514	3770.928892];
    else
        V_ink(:,j_Exp) = [0.756402733	1.0751432	1.3971412	1.728807	2.0436	2.372680733	2.6986786	3.028589	3.359495	3.689999267	4.022913933	4.359566467	4.695162267	5.0318672	5.371951933	5.712272467	6.053859333	6.396878467	6.7413124	7.0870738	7.436206267];
        
        V_tot(:,j_Exp) = [1808.695744	1812.360043	1816.061791	1819.874682	1823.4936	1827.276772	1831.024504	1834.817214	1838.62137	1842.420908	1846.248156	1850.118375	1853.976446	1857.847267	1861.756944	1865.669331	1869.596276	1873.539687	1877.499362	1881.474299	1885.48799];
    end
end

pho = 0.009948208; % Dilution of ink
Lambda = 0.15;

for j_chan = 1:14
    % mua_water = 0.00262955;
    
    E_ink = [7300.741	7104.910	6960.671	6797.173	6538.264	6199.993	6264.458	6176.491	6182.630	6090.587	6126.242	5911.400	5624.083	5388.556];

    for j_Exp = 1:2
        mua_nom(:,j_chan,j_Exp) = mua_Water(j_chan) + V_ink(:,j_Exp) * pho ./ V_tot(:,j_Exp) * E_ink(j_chan) * (1-Lambda);
    end
end
end