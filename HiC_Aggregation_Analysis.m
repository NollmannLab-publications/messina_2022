%% HiC Aggregate Peak Analysis 

% Olivier Messina from Marcelo Nollmann lab
% Centre de Biologie Structurale Montpellier FRANCE

% Execute with Matlab 2019b.

% This script required 4 files that contain the genomic coordinated of the binding sites of interest 
% separated by chromosomes.

% Example : 
% - BEAF_CHIP_Coordinateslibraries_2L_DM6.mat 
% - BEAF_CHIP_Coordinateslibraries_2R_DM6.mat 
% - BEAF_CHIP_Coordinateslibraries_3L_DM6.mat 
% - BEAF_CHIP_Coordinateslibraries_3R_DM6.mat

% Each files is a cell array divided in 4 columns

% Column 1 A number that corresponds to the identity of the bin 
% Column 2 The genomic coordinate of the bin of interest
% Column 3 The distance from the center of the bin of interest 
% Column 4 The chromosomal arm

% This script use Hi-C data from Hug et al. ArrayExpress: E-MTAB-4918

%% Compute data section 1/3

clear all 
% Complete here with HiC you want to process {'NC14','NC14','NC14','NC14_ZLD','NC14_alpha_amanitin','NC14_triptolide','NC14_water'}
HiC_files = {'NC14'};
% Chromosomes to process {'2L','2R','3L','3R','X'}
Chromosomes ={'2L','2R','3L','3R'};
% Folder containing library 
Libraries='/mnt/PALM_dataserv/DATA/Olivier/Thesis/Insulators_Project/Paper_Insulator/Code_Github';
cd(Libraries)
myObj = functionsContainer;

cd(Libraries)
fileList = dir('*.mat'); % check for all file in 
for i=1:size(fileList,1)
    for j=1:size(HiC_files,2)
        Chr=Chromosomes{i};
        NC=HiC_files{j};
        cd(Libraries)
        load(fileList(i).name);
        
        % Formating data
        Coordinateslibraries_bis=[];
        for k=1:length(Coordinateslibraries)
            Coordinateslibraries_bis(k,1)= Coordinateslibraries{k,1};
            Coordinateslibraries_bis(k,2)= Coordinateslibraries{k,2};
        end
        Coordinateslibraries=Coordinateslibraries_bis;            
   
        tic
        Matrix_all = myObj.HiC_calibration(Chr,NC,Coordinateslibraries);     
        toc
                 
        str = strcat('Matrix_all','_',Chr);
        NewNameVariable = [str,'= Matrix_all;'];
        eval(NewNameVariable);
             
    end 
end


%% Prepare for ploting section 2/3

cd(Libraries)

Coordinates_2L = load('BEAF32_old_GSM409067_top_500_Coordinateslibraries_2L_DM6_Shifted_wo_borders.mat');
Coordinates_2L = Coordinates_2L.Coordinateslibraries;
Coordinates_2R = load('BEAF32_old_GSM409067_top_500_Coordinateslibraries_2R_DM6_Shifted_wo_borders.mat');
Coordinates_2R = Coordinates_2R.Coordinateslibraries;
Coordinates_3L = load('BEAF32_old_GSM409067_top_500_Coordinateslibraries_3L_DM6_Shifted_wo_borders.mat');
Coordinates_3L = Coordinates_3L.Coordinateslibraries;
Coordinates_3R = load('BEAF32_old_GSM409067_top_500_Coordinateslibraries_3R_DM6_Shifted_wo_borders.mat');
Coordinates_3R = Coordinates_3R.Coordinateslibraries;
Coordinates_all =[Coordinates_2L;Coordinates_2R;Coordinates_3L;Coordinates_3R];

Matrix_autosome=[Matrix_all_2L;Matrix_all_2R;Matrix_all_3L;Matrix_all_3R];
Aggregation_matrix = myObj.Aggregation_matrix(Matrix_autosome,Coordinates_all);


%% Plot section 3/3

Conditions_1 =[-100;-90;-95;-85;-80;-75;-70;-65;-60;-55;-50;-45;-40;-35;-30;-25;-20;-15;-10;-5;0;5;10;15;20;25;30;35;40;45;50;55;60;65;70;75;80;85;90;95;100];
Conditions_2 = Conditions_1;
imagesc(Aggregation_matrix);
%// Cosmetic changes for the axes
set(gca, 'XTick', 1:2:41); 
set(gca, 'YTick', 1:2:41);
set(gca, 'XTickLabel', Conditions_1);
set(gca, 'YTickLabel', Conditions_2);
caxis([-0.25 0.25])
colormap(coolwarm(256))
colorbar;
set(gcf,'Color','w')
axis square

