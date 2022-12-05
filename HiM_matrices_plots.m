%% Basic HiM Analysis 

% Olivier Messina from Marcelo Nollmann lab
% Centre de Biologie Structurale Montpellier FRANCE

% Run with Matlab 2019b.

% This script requires Hi-M datasets generated in this study
% See our github page for download instructions: https://github.com/NollmannLab/messina_2022

% Colormap used in the paper are generated with the getPyPlot_cMap function 
% https://fr.mathworks.com/matlabcentral/fileexchange/68239-pycolormap4matlab

%% Section 1 Drosophila melanogaster, dpp locus, nuclear cycle 14

clear all
close all
clc

% Load dataset 
Input_data = '/mnt/PALM_dataserv/DATA/Olivier/Thesis/Insulators_Project/Paper_Insulator/DATA_HiM_submited' 
cd(Input_data)
data = load('DPP_NC14_Filtered_3D.mat');
data = data.Matrix_distance_NC14_filtered;
data_median_filtered = nanmedian(data,3);
Number_cells = size(data,3);

% Plot HiM distance matrix
figure
imagesc(data_median_filtered);
axis square
caxis([350 900])
colormap(flipud(coolwarm))
colorbar
title(['PWD matrix NC14 Dpp Cells : ' num2str(Number_cells)])

% Plot HiM contact matrix 
Threshold = 250; % Threshold to consider a proximity to be a contact, in nm.
[NC14_Matrix_SC_Binary_Probrability] = Sc_analysis_HiM(data,Threshold);

figure
imagesc(NC14_Matrix_SC_Binary_Probrability);
axis square
caxis([2 15])
colormap((jet))
% colormap((PiYG))
colorbar
title(['Contact matrix NC14 Dpp Cells : ' num2str(Number_cells)])

% Plot IBPs vs Not IBPs interaction frequency

RTs_IBPs= [2,3,4,9,13,19,20,32,33,34]; %Insulator 
RTs_not_IBPs = [1,5,6,8,11,14,17,21,28,29]; %Not Insulator 
nRTs = size(RTs_IBPs,2);

data_insulators = data(RTs_IBPs,RTs_IBPs,:);
data_not_insulators = data(RTs_not_IBPs,RTs_not_IBPs,:);
Matrix_Contact_TSH = []; 

a=1;
for i=50:50:1500
    Threshold = i;
    
    [Contact_matrix_IBPs] = Sc_analysis_HiM(data_insulators,Threshold);
    [Contact_matrix_not_IBPs] = Sc_analysis_HiM(data_not_insulators,Threshold);
    
    Matrix_Contact_TSH(a,1)=Threshold;
    
    
    for j=1:size(Contact_matrix_IBPs,1)
        for k=1:size(Contact_matrix_IBPs,1)
            if j==k
                Contact_matrix_IBPs(j,k)=NaN;
                Contact_matrix_not_IBPs(j,k)=NaN;         
            end
        end
    end
    
    Matrix_Contact_TSH(a,2)=nanmean(Contact_matrix_IBPs,'all');
    Matrix_Contact_TSH(a,3)=nanmean(Contact_matrix_not_IBPs,'all');
    
    a=a+1;
end

plot(Matrix_Contact_TSH(:,1),Matrix_Contact_TSH(:,2),'LineWidth',3,'Color',[0 0.5 0])
xticks([0:100:1500])
xtickangle(45)
hold on 
plot(Matrix_Contact_TSH(:,1),Matrix_Contact_TSH(:,3),'LineWidth',3,'Color',[0.83 0.667 0])
hold off

%% Section 2 Drosophila melanogaster, dpp locus, nuclear cycle 12

% Load dataset 
data = load('DPP_NC12_Filtered_3D.mat');
data = data.Matrix_distance_NC12_filtered;
data_median_filtered = nanmedian(data,3);
Number_cells = size(data,3);

% Plot HiM distance matrix
figure
imagesc(data_median_filtered);
axis square
caxis([350 900])
colormap(flipud(coolwarm))
colorbar
title(['PWD matrix NC12 Dpp Cells : ' num2str(Number_cells)])

% Plot HiM contact matrix 
Threshold = 250; % Threshold to consider a proximity to be a contact 
[NC14_Matrix_SC_Binary_Probrability] = Sc_analysis_HiM(data,Threshold);

figure
imagesc(NC14_Matrix_SC_Binary_Probrability);
axis square
caxis([2 15])
colormap((jet))
% colormap((PiYG))
colorbar
title(['Contact matrix NC12 Dpp Cells : ' num2str(Number_cells)])

% Plot IBPs vs Not IBPs interaction frequency

RTs_IBPs= [2,3,4,9,13,19,20,32,33,34]; %Insulator 
RTs_not_IBPs = [1,5,6,8,11,14,17,21,28,29]; %Not Insulator 
nRTs = size(RTs_IBPs,2);

data_insulators = data(RTs_IBPs,RTs_IBPs,:);
data_not_insulators = data(RTs_not_IBPs,RTs_not_IBPs,:);
Matrix_Contact_TSH = []; 

a=1;
for i=50:50:1500
    Threshold = i;
    
    [Contact_matrix_IBPs] = Sc_analysis_HiM(data_insulators,Threshold);
    [Contact_matrix_not_IBPs] = Sc_analysis_HiM(data_not_insulators,Threshold);
    
    Matrix_Contact_TSH(a,1)=Threshold;
    
    
    for j=1:size(Contact_matrix_IBPs,1)
        for k=1:size(Contact_matrix_IBPs,1)
            if j==k
                Contact_matrix_IBPs(j,k)=NaN;
                Contact_matrix_not_IBPs(j,k)=NaN;         
            end
        end
    end
    
    Matrix_Contact_TSH(a,2)=nanmean(Contact_matrix_IBPs,'all');
    Matrix_Contact_TSH(a,3)=nanmean(Contact_matrix_not_IBPs,'all');
    
    a=a+1;
end

plot(Matrix_Contact_TSH(:,1),Matrix_Contact_TSH(:,2),'LineWidth',3,'Color',[0 0.5 0])
xticks([0:100:1500])
xtickangle(45)
hold on 
plot(Matrix_Contact_TSH(:,1),Matrix_Contact_TSH(:,3),'LineWidth',3,'Color',[0.83 0.667 0])
hold off

%% Functions 

function [Matrix_SC_Binary_Probrability] = Sc_analysis_HiM(data_raw,Threshold)

nRTs = 34;
data_median = nanmedian(data_raw,3);

% Diagonal = 0 
for i=1:size(data_median,1)
    for j=1:size(data_median,2)
        if i==j
            data_median(i,j)=0;
        end
    end
end

Matrix_SC_Binary = [];

for i=1:size(data_raw,3)
    % Binarize the matrix
    Matrix = NaN(nRTs,nRTs);
    Matrix(data_raw(:,:,i) < Threshold)=1;
    Matrix_SC_Binary(:,:,i) = Matrix;
    
    % Cell with combination of barcodes
    Matrix = NaN(nRTs,nRTs);
    Matrix(data_raw(:,:,i) > 0)=1;
    Matrix_SC_Binary_nb(:,:,i) = Matrix;
    
    % Fill SC with mean 
    Matrix = data_raw(:,:,i);
    Matrix(isnan(data_raw(:,:,i))) = data_median(isnan(data_raw(:,:,i)));
    data_filled(:,:,i) = Matrix;    
    
end

Matrix_SC_Binary_Sum = nansum(Matrix_SC_Binary,3);
Matrix_SC_Binary_Sum_barcodes_present = nansum(Matrix_SC_Binary_nb,3);
Matrix_SC_Binary_Probrability = (Matrix_SC_Binary_Sum./Matrix_SC_Binary_Sum_barcodes_present)*100;

for i=1:size(Matrix_SC_Binary_Probrability,1)
    for j=1:size(Matrix_SC_Binary_Probrability,1)
        if i==j
            Matrix_SC_Binary_Probrability(i,j)=100;
        end
    end
end

end

function map = coolwarm(m)
%COOLWARM cool-warm color map
%   COOLWARM(M) returns an M-by-3 matrix containing a colormap with cool-to-warm
%   colors, as commonly used in Paraview.
%   COOLWARM, by itself, is the same length as the current colormap.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(coolwarm)
%
%   Colormap is based on the colors used by the freeware program Paraview.
%   The color table used here is CoolWarmUChar33.csv, from
%   http://www.sandia.gov/~kmorel/documents/ColorMaps/
%   Reference: Moreland, Kenneth, 2009, Diverging Color Maps for Scientific 
%   Visualization, in Proceedings of the 5th International Symposium on 
%   Visual Computing.
%   The Matlab code is after haxby.m by Kelsey Jordahl, Marymount Manhattan
%   College.
%
%   See also HSV, GRAY, PINK, COOL, BONE, COPPER, FLAG, HOT
%   COLORMAP, RGBPLOT, HAXBY.

% Mark Brandon
% Yale University
% Time-stamp: <Aug 20 2012>

%% Check inputs
narginchk(0,1);

if nargin == 1
    validateattributes(m,{'numeric'},{'numel',1});
end

%% Begin Function
if nargin < 1, m = size(get(gcf,'colormap'),1); end
c=[59	76	192;
60	78	194;
61	80	195;
62	81	197;
63	83	198;
64	85	200;
66	87	201;
67	88	203;
68	90	204;
69	92	206;
70	93	207;
71	95	209;
73	97	210;
74	99	211;
75	100	213;
76	102	214;
77	104	215;
79	105	217;
80	107	218;
81	109	219;
82	110	221;
84	112	222;
85	114	223;
86	115	224;
87	117	225;
89	119	226;
90	120	228;
91	122	229;
93	123	230;
94	125	231;
95	127	232;
96	128	233;
98	130	234;
99	131	235;
100	133	236;
102	135	237;
103	136	238;
104	138	239;
106	139	239;
107	141	240;
108	142	241;
110	144	242;
111	145	243;
112	147	243;
114	148	244;
115	150	245;
116	151	246;
118	153	246;
119	154	247;
120	156	247;
122	157	248;
123	158	249;
124	160	249;
126	161	250;
127	163	250;
129	164	251;
130	165	251;
131	167	252;
133	168	252;
134	169	252;
135	171	253;
137	172	253;
138	173	253;
140	174	254;
141	176	254;
142	177	254;
144	178	254;
145	179	254;
147	181	255;
148	182	255;
149	183	255;
151	184	255;
152	185	255;
153	186	255;
155	187	255;
156	188	255;
158	190	255;
159	191	255;
160	192	255;
162	193	255;
163	194	255;
164	195	254;
166	196	254;
167	197	254;
168	198	254;
170	199	253;
171	199	253;
172	200	253;
174	201	253;
175	202	252;
176	203	252;
178	204	251;
179	205	251;
180	205	251;
182	206	250;
183	207	250;
184	208	249;
185	208	248;
187	209	248;
188	210	247;
189	210	247;
190	211	246;
192	212	245;
193	212	245;
194	213	244;
195	213	243;
197	214	243;
198	214	242;
199	215	241;
200	215	240;
201	216	239;
203	216	238;
204	217	238;
205	217	237;
206	217	236;
207	218	235;
208	218	234;
209	219	233;
210	219	232;
211	219	231;
213	219	230;
214	220	229;
215	220	228;
216	220	227;
217	220	225;
218	220	224;
219	220	223;
220	221	222;
221	221	221;
222	220	219;
223	220	218;
224	219	216;
225	219	215;
226	218	214;
227	218	212;
228	217	211;
229	216	209;
230	216	208;
231	215	206;
232	215	205;
232	214	203;
233	213	202;
234	212	200;
235	212	199;
236	211	197;
236	210	196;
237	209	194;
238	209	193;
238	208	191;
239	207	190;
240	206	188;
240	205	187;
241	204	185;
241	203	184;
242	202	182;
242	201	181;
243	200	179;
243	199	178;
244	198	176;
244	197	174;
245	196	173;
245	195	171;
245	194	170;
245	193	168;
246	192	167;
246	191	165;
246	190	163;
246	188	162;
247	187	160;
247	186	159;
247	185	157;
247	184	156;
247	182	154;
247	181	152;
247	180	151;
247	178	149;
247	177	148;
247	176	146;
247	174	145;
247	173	143;
247	172	141;
247	170	140;
247	169	138;
247	167	137;
247	166	135;
246	164	134;
246	163	132;
246	161	131;
246	160	129;
245	158	127;
245	157	126;
245	155	124;
244	154	123;
244	152	121;
244	151	120;
243	149	118;
243	147	117;
242	146	115;
242	144	114;
241	142	112;
241	141	111;
240	139	109;
240	137	108;
239	136	106;
238	134	105;
238	132	103;
237	130	102;
236	129	100;
236	127	99;
235	125	97;
234	123	96;
233	121	95;
233	120	93;
232	118	92;
231	116	90;
230	114	89;
229	112	88;
228	110	86;
227	108	85;
227	106	83;
226	104	82;
225	102	81;
224	100	79;
223	98	78;
222	96	77;
221	94	75;
220	92	74;
218	90	73;
217	88	71;
216	86	70;
215	84	69;
214	82	67;
213	80	66;
212	78	65;
210	75	64;
209	73	62;
208	71	61;
207	69	60;
205	66	59;
204	64	57;
203	62	56;
202	59	55;
200	57	54;
199	54	53;
198	51	52;
196	49	50;
195	46	49;
193	43	48;
192	40	47;
190	37	46;
189	34	45;
188	30	44;
186	26	43;
185	22	41;
183	17	40;
181	11	39;
180	4	38];
%... Interpolate get requested size for color table
pp=1:(m-1)/(size(c,1)-1):m;
r=interp1(pp,c(:,1),1:m);
g=interp1(pp,c(:,2),1:m);
b=interp1(pp,c(:,3),1:m);
%... Normalize to range [0,1], and divide again by maximum value
% to correct for round-off errors associated with the interpolation.
map=[r' g' b']/255;
map = map/max(map(:));
end




