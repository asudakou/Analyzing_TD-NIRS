% Written by Aleh Sudakou (https://www.researchgate.net/profile/Aleh-Sudakou)
% Nałęcz Institute of Biocybernetics and Biomedical Engineering, Polish Academy of Sciences 
% This is one of the codes shared on  https://github.com/asudakou/Analyzing_TD-NIRS
% Last updated: 20 May 2023


%% Description:
% This script can run only on a computer that has both: 
% 1. The raw data
% 2. The custom function for opening raw data ("read_DTOF_SDT.m")
% 
% This script simply reads raw data (.sdt files) and saves into .mat files in specificed variables.
% No filtering is done at this stage 
% 
% The data are only: DTOFs, IRFs, and related information (time channels and collection time) 


%% Load the measured data for 2 experiments with ink, rearrange it and save in a clear format
clear;
clc;
    DTOF = [];
    DTOF_SameND = [];
    DTOF_ConstLayer = [];
    IRF = [];

% Load all of the measured DTOFs and IRFs from .sdt files (these are described in excel file)
    DTOF_Exp1_Blu = zeros(32,16,1024);
    DTOF_Exp1_Red = zeros(32,16,1024);
    DTOF_Exp2_Blu = zeros(32,16,1024);
    DTOF_Exp2_Red = zeros(32,16,1024);

% For experiment #1
    folder = 'C:\Users\asudakou\Documents\01 Work mine\Paper 4 depth-resolved delMua\Ink phantom 2\2022-06-05\EXP1';
    filename = 'Meas_1';
    [ out_blu, out_red, Time_ns, CollectionTime_singleDTOF, where_merge ] = read_DTOF_SDT( folder, filename, 0, 0, 0, 2 ); % Returns all DTOFs
    for j=1:round(size(out_blu,1)/128)
        DTOF_Exp1_Blu(j,:,:) = squeeze(mean(out_blu(1+(j-1)*128:j*128,:,:),1));
        DTOF_Exp1_Red(j,:,:) = squeeze(mean(out_red(1+(j-1)*128:j*128,:,:),1));
    end
    filename = 'IRF_start';
    [ out_blu, out_red, Time_ns, CollectionTime_singleDTOF, where_merge ] = read_DTOF_SDT( folder, filename, 0, 0, 0, 2 );
    for j=1:16
        IRF.Exp1_Blu(:,j) = squeeze(mean(out_blu(:,j,:),1));
        IRF.Exp1_Red(:,j) = squeeze(mean(out_red(:,j,:),1));
    end

% For experiment #2
    folder = 'C:\Users\asudakou\Documents\01 Work mine\Paper 4 depth-resolved delMua\Ink phantom 2\2022-06-05\EXP2';
    filename = 'Meas_1';
    [ out_blu, out_red, Time_ns, CollectionTime_singleDTOF, where_merge ] = read_DTOF_SDT( folder, filename, 0, 0, 0, 2 );
    for j=1:round(size(out_blu,1)/128)
        DTOF_Exp2_Blu(j,:,:) = squeeze(mean(out_blu(1+(j-1)*128:j*128,:,:),1));
        DTOF_Exp2_Red(j,:,:) = squeeze(mean(out_red(1+(j-1)*128:j*128,:,:),1));
    end
    filename = 'IRF_start';
    [ out_blu, out_red, Time_ns, CollectionTime_singleDTOF, where_merge ] = read_DTOF_SDT( folder, filename, 0, 0, 0, 2 );
    for j=1:16
        IRF.Exp2_Blu(:,j) = squeeze(mean(out_blu(:,j,:),1));
        IRF.Exp2_Red(:,j) = squeeze(mean(out_red(:,j,:),1));
    end


% Now sort the loaded DTOFs into proper variables that will be used in the next scripts
    % DTOFs on a layer that we didn't change, at the start and end of experiment
        DTOF_ConstLayer.Exp1_SupL_Blu = squeeze(DTOF_Exp1_Blu([1 32],:,:));
        DTOF_ConstLayer.Exp2_Deep_Blu = squeeze(DTOF_Exp2_Blu([1 32],:,:));


    % DTOFs for Experiment 1
        POSPOS = [2 4 6 8 10 12 14 16 18:22 24:31]; % After each Mua step, we measured 2 times: before and after changing ND filter
                                                    % But only up to 24, then the ND filter is open to maximum
        DTOF.Exp1_TwoL_Blu = DTOF_Exp1_Blu(POSPOS,:,:);
        DTOF.Exp1_Deep_Red = DTOF_Exp1_Red(POSPOS,:,:);

        DTOF_SameND.Exp1_TwoL_Blu = zeros(size(DTOF.Exp1_TwoL_Blu));
        DTOF_SameND.Exp1_TwoL_Blu([2 14],:,:) = DTOF_Exp1_Blu([3 23],:,:);

        DTOF_SameND.Exp1_Deep_Red = zeros(size(DTOF.Exp1_Deep_Red));
        DTOF_SameND.Exp1_Deep_Red(2:9,:,:) = DTOF_Exp1_Red(3:2:17,:,:);
    
    % DTOFs for Experiment 2
        POSPOS = [2 4 6 8 10 12 14 16 18:1:30];% After each Mua step, we measured 2 times: before and after changing ND filter
                                               % But only up to 18, then the ND filter is open to maximum
        DTOF.Exp2_TwoL_Blu = DTOF_Exp2_Blu(POSPOS,:,:);
        DTOF.Exp2_SupL_Red = DTOF_Exp2_Red(POSPOS,:,:);

        DTOF_SameND.Exp2_TwoL_Blu  = zeros(size(DTOF.Exp2_TwoL_Blu));
        DTOF_SameND.Exp2_TwoL_Blu(2:9,:,:) = DTOF_Exp2_Blu(3:2:17,:,:);

        DTOF_SameND.Exp2_SupL_Red = zeros(size(DTOF.Exp2_SupL_Red));
        DTOF_SameND.Exp2_SupL_Red(2:7,:,:) = DTOF_Exp2_Red(3:2:13,:,:);


j_chan = 8;
figure(5); clf
temp = squeeze(DTOF.Exp2_TwoL_Blu(1,j_chan,:));
semilogy(Time_ns + (Time_ns(2)-Time_ns(1))*1, temp / max(temp))
hold on
semilogy(Time_ns, temp / max(temp))
hold on
temp = squeeze(DTOF.Exp2_TwoL_Blu(2,j_chan,:));
semilogy(Time_ns, temp / max(temp))
temp = squeeze(DTOF.Exp2_TwoL_Blu(3,j_chan,:));
semilogy(Time_ns, temp / max(temp))
temp = squeeze(DTOF.Exp2_TwoL_Blu(4,j_chan,:));
semilogy(Time_ns, temp / max(temp))
temp = squeeze(DTOF.Exp2_TwoL_Blu(10,j_chan,:));
semilogy(Time_ns, temp / max(temp))
xlim([1.5 4])

root = pwd;
cd('C:\Users\asudakou\Documents\01 Work mine\MatLab codes\Wrote by me\TR analyses\Two-layered analysis\Aleh LMA fitting\Open-access\DATA measured')
save('data_Ink_Pub2023','DTOF','DTOF_SameND','DTOF_ConstLayer','IRF','Time_ns')
cd(root)

% Loading 'data_Ink_Pub2023.mat' :
DTOF;
DTOF_SameND;
DTOF_ConstLayer;
IRF;
Time_ns;

disp('FINISHED READING AND SAVING DATA FOR EXPERIMENTS WITH INK')



%% Load the measured data for 3 experiments with blood
clear;
clc;

to_sum = 9; %  == 3 sec
CollectionTime = to_sum / 3; 

% EXP. 1
folder = 'C:\Users\asudakou\Documents\01 Work mine\Paper 4 depth-resolved delMua\Ink phantom 2\2022-06-03';
save_name = 'data_Blood_Exp1_Pub2023';

    filename = 'Meas_1';
    [ fileoutput1, fileoutput2, Time_ns, CollectionTime_singleDTOF, where_merge ] = read_DTOF_SDT( folder, filename, 0, 0, 0, 2 );
    Exp1_DTOF.TwoL = zeros(floor(size(fileoutput1,1)/to_sum),size(fileoutput1,2),size(fileoutput1,3));
    Exp1_DTOF.Deep = zeros(size(Exp1_DTOF.TwoL));
    c_2 = 1;
    for j = 1:to_sum:size(fileoutput1,1)-to_sum+1
        Exp1_DTOF.TwoL(c_2,:,:) = sum(fileoutput1(j:j+to_sum-1,:,:),1);
    
        Exp1_DTOF.Deep(c_2,:,:) = sum(fileoutput2(j:j+to_sum-1,:,:),1);
        c_2 = c_2 + 1;
    end
    
    filename = 'IRF_start';
    [ IRF_Sup_start, IRF_Deep_start, Time_ns, CollectionTime_singleDTOF, where_merge ] = read_DTOF_SDT( folder, filename, 0, 0, 0, 2 );
    Exp1_IRF.TwoL_start = squeeze(mean(IRF_Sup_start,1));
    Exp1_IRF.Deep_start = squeeze(mean(IRF_Deep_start,1));
    
    filename = 'IRF_end';
    [ IRF_Sup_end, IRF_Deep_end, Time_ns, CollectionTime_singleDTOF, where_merge ] = read_DTOF_SDT( folder, filename, 0, 0, 0, 2 );
    Exp1_IRF.TwoL_end = squeeze(mean(IRF_Sup_end,1));
    Exp1_IRF.Deep_end = squeeze(mean(IRF_Deep_end,1));

    root = pwd;
    cd('C:\Users\asudakou\Documents\01 Work mine\MatLab codes\Wrote by me\TR analyses\Two-layered analysis\Aleh LMA fitting\Open-access\DATA measured')
    save(save_name,'Exp1_DTOF','Exp1_IRF','Time_ns','CollectionTime')
    msgbox(['Saved   ' save_name])
    cd(root)

% EXP. 2
folder = 'C:\Users\asudakou\Documents\01 Work mine\Paper 4 depth-resolved delMua\Ink phantom 2\2022-06-04';
save_name = 'data_Blood_Exp2_Pub2023';

    filename = 'Meas_1';
    [ fileoutput1, fileoutput2, Time_ns, CollectionTime_singleDTOF, where_merge ] = read_DTOF_SDT( folder, filename, 0, 0, 0, 2 );
    Exp2_DTOF.TwoL = zeros(floor(size(fileoutput1,1)/to_sum),size(fileoutput1,2),size(fileoutput1,3));
    Exp2_DTOF.Deep = zeros(size(Exp2_DTOF.TwoL));
    c_2 = 1;
    for j = 1:to_sum:size(fileoutput1,1)-to_sum+1
        Exp2_DTOF.TwoL(c_2,:,:) = sum(fileoutput1(j:j+to_sum-1,:,:),1);
    
        Exp2_DTOF.Deep(c_2,:,:) = sum(fileoutput2(j:j+to_sum-1,:,:),1);
        c_2 = c_2 + 1;
    end
    
    filename = 'IRF_start';
    [ IRF_Sup_start, IRF_Deep_start, Time_ns, CollectionTime_singleDTOF, where_merge ] = read_DTOF_SDT( folder, filename, 0, 0, 0, 2 );
    Exp2_IRF.TwoL_start = squeeze(mean(IRF_Sup_start,1));
    Exp2_IRF.Deep_start = squeeze(mean(IRF_Deep_start,1));
    
    filename = 'IRF_end';
    [ IRF_Sup_end, IRF_Deep_end, Time_ns, CollectionTime_singleDTOF, where_merge ] = read_DTOF_SDT( folder, filename, 0, 0, 0, 2 );
    Exp2_IRF.TwoL_end = squeeze(mean(IRF_Sup_end,1));
    Exp2_IRF.Deep_end = squeeze(mean(IRF_Deep_end,1));

    root = pwd;
    cd('C:\Users\asudakou\Documents\01 Work mine\MatLab codes\Wrote by me\TR analyses\Two-layered analysis\Aleh LMA fitting\Open-access\DATA measured')
    save(save_name,'Exp2_DTOF','Exp2_IRF','Time_ns','CollectionTime')
    msgbox(['Saved   ' save_name])
    cd(root)

% EXP. 3
folder = 'C:\Users\asudakou\Documents\01 Work mine\Paper 4 depth-resolved delMua\Ink phantom 2\2022-06-02';
save_name = 'data_Blood_Exp3_Pub2023';

    filename = 'Meas_1';
    [ fileoutput1, fileoutput2, Time_ns, CollectionTime_singleDTOF, where_merge ] = read_DTOF_SDT( folder, filename, 0, 0, 0, 2 );
    Exp3_DTOF.TwoL = zeros(floor(size(fileoutput1,1)/to_sum),size(fileoutput1,2),size(fileoutput1,3));
    Exp3_DTOF.Deep = zeros(size(Exp3_DTOF.TwoL));
    c_2 = 1;
    for j = 1:to_sum:size(fileoutput1,1)-to_sum+1
        Exp3_DTOF.TwoL(c_2,:,:) = sum(fileoutput1(j:j+to_sum-1,:,:),1);
    
        Exp3_DTOF.Deep(c_2,:,:) = sum(fileoutput2(j:j+to_sum-1,:,:),1);
        c_2 = c_2 + 1;
    end
    
    filename = 'IRF_start';
    [ IRF_Sup_start, IRF_Deep_start, Time_ns, CollectionTime_singleDTOF, where_merge ] = read_DTOF_SDT( folder, filename, 0, 0, 0, 2 );
    Exp3_IRF.TwoL_start = squeeze(mean(IRF_Sup_start,1));
    Exp3_IRF.Deep_start = squeeze(mean(IRF_Deep_start,1));
    
    filename = 'IRF_end';
    [ IRF_Sup_end, IRF_Deep_end, Time_ns, CollectionTime_singleDTOF, where_merge ] = read_DTOF_SDT( folder, filename, 0, 0, 0, 2 );
    Exp3_IRF.TwoL_end = squeeze(mean(IRF_Sup_end,1));
    Exp3_IRF.Deep_end = squeeze(mean(IRF_Deep_end,1));

    root = pwd;
    cd('C:\Users\asudakou\Documents\01 Work mine\MatLab codes\Wrote by me\TR analyses\Two-layered analysis\Aleh LMA fitting\Open-access\DATA measured')
    save(save_name,'Exp3_DTOF','Exp3_IRF','Time_ns','CollectionTime')
    msgbox(['Saved   ' save_name])
    cd(root)

disp('FINISHED READING AND SAVING DATA FOR EXPERIMENTS WITH BLOOD')
