classdef functionsContainer
	methods
        function  Matrix_all = HiC_calibration(obj,Chr,NC,Coordinateslibraries)

		X = sprintf('Normalizing chromosome %s using HiC %s be patient.',Chr,NC);
		disp(X)

		%                ================================
		%                = Chromosome X from 0 to 22.4  =
		%                = Chromosome 2L from 0 to 23   =
		%                = Chromosome 2R from 0 to 21.1 =
		%                = Chromosome 3L from 0 to 24.5 =
		%                = Chromosome 3R fron 0 to 27.9 = 
		%                ================================

		% read file
		switch NC 
			case 'NC12'
				cd('/mnt/PALM_dataserv/DATA/Commun/Hi-C/Drosophila/HiCdata_Hug')
				%disp('You Choose HiC NC12')
				fid = fopen('nuclear_cycle_12_repl_merged_5kb.tsv');
			case 'NC13'
				cd('/mnt/PALM_dataserv/DATA/Commun/Hi-C/Drosophila/HiCdata_Hug')
				%disp('You Choose HiC NC13')
				fid = fopen('nuclear_cycle_13_repl_merged_5kb.tsv');
			case 'NC14'
				cd('/mnt/PALM_dataserv/DATA/Commun/Hi-C/Drosophila/HiCdata_Hug')
				%disp('You Choose HiC NC14')
				fid = fopen('nuclear_cycle_14_repl_merged_5kb.tsv');
			case 'NC14_ZLD'
				cd('/mnt/PALM_dataserv/DATA/Commun/Hi-C/Drosophila/HiCdata_Hug')
				%disp('You Choose HiC NC14_ZLD')
				fid = fopen('nuclear_cycle_14_sh_zld.tsv');
			case 'NC14_triptolide'
				cd('/mnt/PALM_dataserv/DATA/Commun/Hi-C/Drosophila/HiCdata_Hug')
				%disp('You Choose HiC NC14_triptolide')
				fid = fopen('nuclear_cycle_14_triptolide_injected_5kb.tsv');
			case 'NC14_water'
				cd('/mnt/PALM_dataserv/DATA/Commun/Hi-C/Drosophila/HiCdata_Hug')
				%disp('You Choose HiC NC14_water')
				fid = fopen('nuclear_cycle_14_water_injected_5kb.tsv');
			case 'NC14_alpha_amanitin'
				cd('/mnt/PALM_dataserv/DATA/Commun/Hi-C/Drosophila/HiCdata_Hug')
				%disp('You Choose HiC NC14_alpha_amanitin')
				fid = fopen('nuclear_cycle_14_alpha_amanitin.tsv');
			case 'Mitosis'
				cd('/mnt/PALM_dataserv/DATA/Commun/Hi-C/Drosophila/HiCdata_Hug')
				%disp('mitotic_nuclear_cycle1-14_repl_merged.tsv')
				fid = fopen('mitotic_nuclear_cycle1-14_repl_merged.tsv');
			case '3-4h'
				cd('/mnt/PALM_dataserv/DATA/Commun/Hi-C/Drosophila/HiCdata_Hug')
				fid = fopen('3-4h_merged_5kb.tsv');
			otherwise
				disp('HiC not found')
		end


		C = textscan(fid, '%s %f %f %s %f %f %f %f', 'HeaderLines', 1);
		fclose(fid);
		Matrix_numel = cell2mat(C([2 3 5 6 7 8])); %2:start1 %3:end1 %5:start2 %6:end2 %7:count %8:balanced
		Matrix_all=[num2cell(Matrix_numel),C{:,1},C{:,4}]; %c1 start1 %c2 end %c3 start2 %c4 end2 %c5 count %c6 balanced %c7 chr1 %c8 chr2

		disp('Dataset loaded')

		%Intrachromosomal contact 

		Bin_Chr_1 = squeeze(Matrix_all(:,7));
		Bin_Chr_2 = squeeze(Matrix_all(:,8));
		idx = find(contains(Bin_Chr_1,Chr) & contains(Bin_Chr_2,Chr));
		Matrix_all=Matrix_all(idx,:);

		Distance_Bin1_Start = cell2mat(squeeze(Matrix_all(:,1)));
		Distance_Bin1_End = cell2mat(squeeze(Matrix_all(:,2)));
		Distance_Bin2_Start = cell2mat(squeeze(Matrix_all(:,3)));
		Distance_Bin2_End = cell2mat(squeeze(Matrix_all(:,4)));
		Distances = abs(((Distance_Bin1_Start+Distance_Bin1_End)/2)-((Distance_Bin2_Start+Distance_Bin2_End)/2));
		Matrix_all(:,9) = num2cell(Distances);
		Matrix_distances_counts(:,1) = (squeeze(Matrix_all(:,9)));
		Matrix_distances_counts(:,2) = (squeeze(Matrix_all(:,5)));
		Matrix_distances_counts = cell2mat(Matrix_distances_counts);
		Calibration = grpstats(Matrix_distances_counts,Matrix_distances_counts(:,1));
		Calibration = Calibration(2:end,:);
		f = fit(Calibration(:,1),Calibration(:,2),'power2');
		Calibration(:,3)=f.a*Calibration(:,1).^f.b+f.c;
		disp('Calibration Done HiC')


		for i=1:size(Matrix_all,1)
			for j=1:length(Calibration)
				if Matrix_all{i,9} == Calibration(j,1)
				    Matrix_all{i,10} = Calibration(j,2); % Change to use calibration Olivier change to 2 to 3
				    Matrix_all{i,11} = log2(Matrix_all{i,5}/Matrix_all{i,10});
				break
				end
			end
		end

		Bin1 = cell2mat(squeeze(Matrix_all(:,1)));
		Bin2 = cell2mat(squeeze(Matrix_all(:,2)));
		Bin3 = cell2mat(squeeze(Matrix_all(:,3)));
		Bin4 = cell2mat(squeeze(Matrix_all(:,4)));

		disp('Indexing Bins')

		% Indexing bin 
		Matrix_all{1,12}='';
		for i =1:size(Coordinateslibraries,1)
			idx = find(Bin1 <= Coordinateslibraries(i,2) & Bin2 >= Coordinateslibraries(i,2));
				for j=1:size(idx,1)
				    if isempty(Matrix_all{idx(j),12}) == 0
				        Matrix_all{idx(j),12} = cat(1,Matrix_all{idx(j),12},Coordinateslibraries(i,1));
				    elseif isempty(Matrix_all{idx(j),12}) == 1
				        Matrix_all{idx(j),12} = Coordinateslibraries(i,1);
				    end
				end
		end

		Matrix_all{1,13}='';
		for i =1:size(Coordinateslibraries,1)
			idx = find(Bin3 <= Coordinateslibraries(i,2) & Bin4 >= Coordinateslibraries(i,2));
				for j=1:size(idx,1)
				    if isempty(Matrix_all{idx(j),13}) == 0
				        Matrix_all{idx(j),13} = cat(1,Matrix_all{idx(j),13},Coordinateslibraries(i,1));
				    elseif isempty(Matrix_all{idx(j),13}) == 1
				        Matrix_all{idx(j),13} = Coordinateslibraries(i,1);
				    end
				end
		end

				      
		% Reorder Matrix_all /_\ THIS IS YOUR OUTPUT /_\ 
		clear A
		A(:,1)= Matrix_all(:,7); % #Colum 1 Chromosome
		A(:,2)= {NC}; % #Colum 2 NC
		A(:,3)= Matrix_all(:,1); % #Colum 3 Bin 1 Start
		A(:,4)= Matrix_all(:,2); % #Colum 4 Bin 1 End
		A(:,5)= Matrix_all(:,3); % #Colum 5 Bin 2 Start 
		A(:,6)= Matrix_all(:,4); % #Colum 6 Bin 2 End 
		A(:,7)= Matrix_all(:,12); % #Colum 7 Bin 1 (RT)
		A(:,8)= Matrix_all(:,13); % #Colum 8 Bin 2 (RT)
		A(:,9)= Matrix_all(:,9); % #Colum 9 Distance btw Bin 1 and 2 
		A(:,10)= Matrix_all(:,5); % #Colum 10 Observed number of contacts 
		A(:,11)= Matrix_all(:,10); % #Colum 11 Expected number of contacts 
		A(:,12)= Matrix_all(:,11); % #Colum 12 Log2 (Obs/Exp)  
		A(:,13)= Matrix_all(:,6); % #Colum 13 Balancing  
		Matrix_all = A;

		idx = find(~cellfun(@isempty,Matrix_all(:,11)) & ~cellfun(@isempty,Matrix_all(:,12)));
		Matrix_all = Matrix_all(idx,:);

		idx = find(~cellfun(@isempty,Matrix_all(:,7)) & ~cellfun(@isempty,Matrix_all(:,8)));
		Matrix_all = Matrix_all(idx,:);

		disp('Process finished')
		end

		function Matrix_shifted =Aggregation_matrix(obj,Matrix_autosome,Coordinates_all)

		%Remove self ligated fragement 
		deleterow = false(1,size(Matrix_autosome,1));
		for i=1:size(Matrix_autosome,1)
			if Matrix_autosome{i,9}==0
				deleterow(1,i)=1;
			end
		end
		Matrix_autosome(deleterow,:)=[];

		Bin3 = cell2mat(squeeze(Coordinates_all(:,1)));
		Matrix_autosome{1,14}='';
		Matrix_autosome{1,15}='';
		disp('Ploting 1/4')
		tic
		for i=1:size(Matrix_autosome,1)
			%i*100/size(Matrix_autosome,1)
			A=Matrix_autosome{i,7};
			B=Matrix_autosome{i,8};
			for j=1:size(A,1)
				idx1 = find(Bin3==A(j));
				Matrix_autosome{i,14}=cat(1,Matrix_autosome{i,14},Coordinates_all(idx1,3));
			end
			for k=1:size(B,1)
				idx2 = find(Bin3==B(k));
				Matrix_autosome{i,15}=cat(1,Matrix_autosome{i,15},Coordinates_all(idx2,3));
			end
		end
		toc      
		clearvars -except Coordinates_all Matrix_autosome Figures Libraries

		Bin1=cell(size(Matrix_autosome,1),1);
		Bin2=cell(size(Matrix_autosome,1),1);
		disp('Ploting 2/4')
		tic

		for i=1:size(Matrix_autosome,1)
			bin_a = (Matrix_autosome{i,14});
			bin_b = (Matrix_autosome{i,15});
			Bin1{i,1}=cell2mat(bin_a);
			Bin2{i,1}=cell2mat(bin_b);
		end

		Conditions_1 =[-100000;-95000;-90000;-85000;-80000;-75000;-70000;-65000;-60000;-55000;-50000;-45000;-40000;-35000;-30000;-25000;-20000;-15000;-10000;-5000;0;5000;10000;15000;20000;25000;30000;35000;40000;45000;50000;55000;60000;65000;70000;75000;80000;85000;90000;95000;100000];
		Conditions_2 = Conditions_1;

		disp('Ploting 3/4')

		for i=1:size(Conditions_1,1)
			i*100/size(Conditions_1,1)
			

			for j=i:size(Conditions_2,1)
			 
				if i==j 
				    a=1;
				    test1=[];
				    test2=[];
				    
				    for p=1:size(Matrix_autosome,1)
				        p;
				        test1(a)=sum(ismember(Bin1{p},Conditions_1(i)));
				        test2(a)=sum(ismember(Bin2{p},Conditions_2(j)));
				        a=a+1;
				    end

				    [idx1,B]=(find(test1 ~=0 & test2 ~=0));

				    Matrix_shifted(i,j)= mean(cellfun(@mean,Matrix_autosome(B,12)));
				
				elseif i~=j
				    test3=[];
				    test4=[];
				    test5=[];
				    test6=[];
				    a=1;
				    for p=1:size(Matrix_autosome,1)
				        p;
				        test3(a)=sum(ismember(Bin1{p},Conditions_1(i)));
				        test4(a)=sum(ismember(Bin2{p},Conditions_2(j)));
				        a=a+1;
				    end
				    
				    [idx2,C]=(find(test3 ~=0 & test4 ~=0));
				    a=1;
				    for p=1:size(Matrix_autosome,1)
				        p;
				        test5(a)=sum(ismember(Bin1{p},Conditions_1(j)));
				        test6(a)=sum(ismember(Bin2{p},Conditions_2(i)));
				        a=a+1;
				    end
				    
				    [idx3,D]=(find(test5 ~=0 & test6 ~=0));
				 
				    idx_all = [transpose(C);transpose(D)];
				    Matrix_shifted(i,j)= mean(cellfun(@mean,Matrix_autosome(idx_all,12)));
                end
            end
        end
        end
        function cmp = getPyPlot_cMap(nam,n,keepAlpha,pyCmd)
            % cmp = getPyPlot_cMap(nam [, n, keepAlpha, pyCmd])
            %
            %
            % ::INPUT::
            %
            % nam:      Colormap name from matplotlib library (case sensitve!). See
            %           below. Alternatively, '!GetNames' returns a cellstring of
            %           available colormaps.
            % n:        Number of Colors; defaults to 128
            % keepAlpha: Switch to keep the alpha channel of the colormap (4th colum);
            %           defaults to false. If true, a Nx4 matrix is returned in cmp
            %           (instead of Nx3).
            % pyCmd:    python command; defaults to 'python'
            %
            % 
            % ::OUTPUT::
            %
            % cmp       A Nx3 (Nx4 if keepAlpha is true) double array of RGB(A) values.
            %
            %
            % Colormap name can be one of the following:
            %   ###  NOTE: the set of available colormaps  ###
            %   ###  depends on your python installation!  ###
            %
            % Accent        	Accent_r      	afmhot        	afmhot_r      
            % autumn        	autumn_r      	binary        	binary_r      
            % Blues         	Blues_r       	bone          	bone_r        
            % BrBG          	BrBG_r        	brg           	brg_r         
            % BuGn          	BuGn_r        	BuPu          	BuPu_r        
            % bwr           	bwr_r         	cividis       	cividis_r     
            % CMRmap        	CMRmap_r      	cool          	cool_r        
            % coolwarm      	coolwarm_r    	copper        	copper_r      
            % cubehelix     	cubehelix_r   	Dark2         	Dark2_r       
            % flag          	flag_r        	gist_earth    	gist_earth_r  
            % gist_gray     	gist_gray_r   	gist_heat     	gist_heat_r   
            % gist_ncar     	gist_ncar_r   	gist_rainbow  	gist_rainbow_r
            % gist_stern    	gist_stern_r  	gist_yarg     	gist_yarg_r   
            % GnBu          	GnBu_r        	gnuplot       	gnuplot2      
            % gnuplot2_r    	gnuplot_r     	gray          	gray_r        
            % Greens        	Greens_r      	Greys         	Greys_r       
            % hot           	hot_r         	hsv           	hsv_r         
            % inferno       	inferno_r     	jet           	jet_r         
            % magma         	magma_r       	nipy_spectral 	nipy_spectral_r
            % ocean         	ocean_r       	Oranges       	Oranges_r     
            % OrRd          	OrRd_r        	Paired        	Paired_r      
            % Pastel1       	Pastel1_r     	Pastel2       	Pastel2_r     
            % pink          	pink_r        	PiYG          	PiYG_r        
            % plasma        	plasma_r      	PRGn          	PRGn_r        
            % prism         	prism_r       	PuBu          	PuBu_r        
            % PuBuGn        	PuBuGn_r      	PuOr          	PuOr_r        
            % PuRd          	PuRd_r        	Purples       	Purples_r     
            % rainbow       	rainbow_r     	RdBu          	RdBu_r        
            % RdGy          	RdGy_r        	RdPu          	RdPu_r        
            % RdYlBu        	RdYlBu_r      	RdYlGn        	RdYlGn_r      
            % Reds          	Reds_r        	seismic       	seismic_r     
            % Set1          	Set1_r        	Set2          	Set2_r        
            % Set3          	Set3_r        	Spectral      	Spectral_r    
            % spring        	spring_r      	summer        	summer_r      
            % tab10         	tab10_r       	tab20         	tab20_r       
            % tab20b        	tab20b_r      	tab20c        	tab20c_r      
            % terrain       	terrain_r     	viridis       	viridis_r     
            % winter        	winter_r      	Wistia        	Wistia_r      
            % YlGn          	YlGn_r        	YlGnBu        	YlGnBu_r      
            % YlOrBr        	YlOrBr_r      	YlOrRd        	YlOrRd_r 
            % 
            % V 1.4; Konrad Schumacher, 02.2021

            if strcmpi(nam,'!GetNames')
                % translate switch to retrieve colormap names into python-switch:
                nam = 'listCMapNames';
            end


            % defaults:
            if ~exist('n','var') || isempty(n)
                n = 128;
            end
            if ~exist('keepAlpha','var') || isempty(keepAlpha)
                keepAlpha = 0;
            end
            if ~exist('pyCmd','var') || isempty(pyCmd)
                pyCmd = 'python';
            end


            % check if python script is present
            pyScript = which('pyplotCMap2txt.py');
            assert(~isempty(pyScript), 'getPyPlot_cMap:PyScriptNotFound', ...
                'Could not find python script (%s).','pyplotCMap2txt.py');

            tmpf = tempname;

            % call python script
            comd = sprintf('%s "%s" %s -o "%s" -n %d', pyCmd, pyScript, nam, tmpf, n);
            [s,m] = system(comd);

            % check if system command ran w/o error
            if s~=0
                baseME = MException('getPyPlot_cMap:SystemCMDFailed', ...
                    'There was an error executing the command\n\t%s\nSystem returned:\n\t%s', ...
                    comd, m);

                errPatterns = {'(?<=ModuleNotFoundError: ).*$' ...
                               'argument cmName: invalid choice: [''\d\w+]+'};
                mtch = regexp(m,errPatterns,'match','once');

                if ~isempty(mtch{1}) % ModuleNotFoundError
                    ME = MException('getPyPlot_cMap:SystemCMDFailed:ModulNotFound', ...
                        'Python is missing a required module: %s', mtch{1});

                elseif ~isempty(mtch{2}) % cmName: invalid choice
                    ME = MException('getPyPlot_cMap:SystemCMDFailed:InvalidChoice', ...
                        'The chosen colormap name was not found in the python installation: %s', ...
                        mtch{2});

                else % UNHANDLED CASE
                    ME = MException('getPyPlot_cMap:SystemCMDFailed:Unhandled', ...
                        'There was an unexpected error while executing the python script. Sorry.');
                end

                throw(baseME.addCause(ME));
            end


            if strcmp(nam,'listCMapNames')
                % cMap names retrieved; read text file
                fid = fopen(tmpf,'r');
                cmp = textscan(fid,'%s');
                fclose(fid);
                cmp = cmp{1};

            else
                % load cMap data from text file
                cmp = load(tmpf,'-ascii');

                if keepAlpha
                else % remove 4th column of alpha values
                    cmp = cmp(:,1:3);
                end
            end


            % delete temp-file
            delete(tmpf);

        end
    end
end

