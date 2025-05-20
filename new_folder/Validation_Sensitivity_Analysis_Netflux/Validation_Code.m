clc;
clear all;
close all;
%% it is different from original version, in a way that this new version 
% reads the data from an excel file

tic
%% address to read the data
It2reachSS = 5;
thr = 0.001;

% it shows the order in the netwrok inputs (model excel file)
netInputs = {'AngII', 'TGFB', 'mechanical' ...
            ,'IL6','IL1','TNFa','NE','PDGF','ET1', 'NP', 'Forskolin'};
inputPrefered_order = {'TGFB','AngII',  'mechanical' ...
    'ET1','PDGF','TNFa','IL1','IL6','NE','Forskolin','NP'};
outputPrefered_order = {'CI', 'CIII', 'aSMA','PAI1',...
    'periostin','fibronectin', 'proliferation','CTGF',...
    'migration', 'EDAFN',  'CImRNA', 'CIIImRNA',...
    'TIMP1','TIMP2','MMP1', 'MMP2', 'MMP9','MMP14','latentTGFB','IL6'};


fname = ['Khorasani_Atrial_Fibroblast_Model.xlsx'];
path = '/Users/najmekhorasani/Library/CloudStorage/Dropbox/Najme Folder/article/Atrial-fibrosis-Model/Github/Validation_Sensitivity_Analysis_Netflux/';
pathname = [path,'models/'];

%% read validation data
flg = 0;
while flg==0
    flg = input('1: atrial Fbs Validation, 2: Cardiac Fbs Validation: ');
    if flg ==1
        fnameVal = 'atrial-like_Validation_Data_TGFB_AngII.xlsx';
    elseif flg==2
        fnameVal = 'cardiacFb_Validation_Data_preOrderInp';
    else
        disp('Not Valid Input!');
        flg=0;
    end
end

pathnameVal = [path,'Validation_Data/'];
file = [pathnameVal,fnameVal];
Val_mat = readtable(file, 'Range', 'A:G');

%% run netflux for baseline condition
ind_equal_to_0_in_w = [];
%% define Wmat
indices_in_orig_w = ind_equal_to_0_in_w;
W_Mat = zeros(1,size(indices_in_orig_w,2));
%% run netflux
[specID, lastValuesCtl] = Netflux_NJ_for_predefined_inp(indices_in_orig_w, ...
            W_Mat, fname,pathname, It2reachSS);

%% determine inputs and their corresponding data
[~, NetNodesInd] = ismember(outputPrefered_order, specID);

sz1 = size(netInputs,2);
sz2 = size(outputPrefered_order,2);

predMat = zeros(sz1,sz2);
predMarks = ones(sz1,sz2)*(-1);

for i = 1:length(netInputs)
    key = netInputs{i};  % e.g., 'TGFB', 'IL1'
    
    inPut_indeces = find(ismember(netInputs, key));
    ind = inPut_indeces;   
    W_Mat_sti = [1,W_Mat];
    indices_in_orig_w_sti = [ind,indices_in_orig_w];

    [~,lastValuesSti] = Netflux_NJ_for_predefined_inp(indices_in_orig_w_sti, ...
            W_Mat_sti, fname,pathname, It2reachSS);

    diffTMP = lastValuesSti(1,NetNodesInd)-lastValuesCtl(1,NetNodesInd);
    predMat(i,:) = diffTMP;

    % Save NetNode and Measurement into struct
    idx = strcmp(Val_mat.Input, key);  % logical index for matching rows
    if sum(idx)>=1
        TMP = Val_mat.NetNode(idx);
        [~, outPut_indeces] = ismember(TMP, specID);
        predTMP = lastValuesSti(1,outPut_indeces)'-lastValuesCtl(1,outPut_indeces)';
        predTMP = (predTMP > thr)*1 + (predTMP < -1*thr)*(-1) + (predTMP<=thr & predTMP >= -1*thr)*0;
        Val_mat.Prediction(idx) = predTMP;
        mesTMP = Val_mat.Measurement_V(idx);
    
        matchTMP = (abs(predTMP-mesTMP)==0)*1 + (abs(predTMP-mesTMP)>0)*0;
        Val_mat.Match(idx) = matchTMP;
        [~, matchIndTMP] = ismember(TMP, outputPrefered_order);
        matchTMP(matchIndTMP==0) =[];
        matchIndTMP(matchIndTMP==0) =[];
        predMarks(i,matchIndTMP) = matchTMP;
    end
    
end

%% to change the mat in a order I want
predMat_ordered = zeros(sz1,sz2);
predMarks_ordered = ones(sz1,sz2)*(-1);

for i = 1:length(inputPrefered_order)
    key = inputPrefered_order{i};  % e.g., 'TGFB', 'IL1'
    
    indTMP = find(ismember(netInputs, key));
    
    predMat_ordered(i,:) = predMat(indTMP,:);
    predMarks_ordered(i,:) = predMarks(indTMP,:);
    
end

%% Define the labels and plot the heatmap
% change the order of inputs
% Define the range of possible values in deltaActivity

% to create mymap
no = 100;
stp = 1/no;
c1 = 0:stp:1;
c2 = 0.99:-1*stp:0;
mymap = [[c1';ones(no,1)] [c1';c2'] [ones(no+1,1);c2'] ];

val1 = -1;  % Lower bound
val2 = 1;   % Upper bound
stp = (val2-val1)*stp/2;
potInterval = val1:stp:val2;

yLabels = inputPrefered_order;
columnLabels = outputPrefered_order;

% to Modify diff_Mat
diff_Mat_mod = predMat_ordered;

diff_Mat_mod(abs(diff_Mat_mod)<.001) = 0;
diff_Mat_mod(diff_Mat_mod>0) = val2;
diff_Mat_mod(diff_Mat_mod<0) = val1;

%% determine the colors:
% Plot heatmap
data = diff_Mat_mod;
figure;
imagesc(data); % Display the matrix
ax = gca;  
% define one tick per row
% Define custom colormap: -1 (green), 0 (white), 1 (red)
customColormap = [21/255 100/255 190/255; % Blue for -1
                  1 1 1; % White for 0
                  1 60/255 60/255]; % Red for 1

colormap(customColormap); % Apply colormap
caxis([-1 1]); % Set color limits to match data values
colorbar; % Add color legend

% grids
xticks([]);
yticks([]);
nc = size(data,2);
nr = size(data,1);
hold on;
% full extents for lines
xL = [0.5, nc+0.5];
yL = [0.5, nr+0.5];
% horizontal lines
for r = 1:nr-1
    yy = r + 0.5;
    plot(xL, [yy yy], 'k-', 'LineWidth', .5);
end
% vertical lines
for c = 1:nc-1
    xx = c + 0.5;
    plot([xx xx], yL, 'k-', 'LineWidth', .5);
end
hold off;
ax.YTick = 1:numel(yLabels);
ax.YTickLabel = yLabels;  % set your input‐node names
ylabel('Input Nodes');    % axis label
% Set the x-axis labels
ax.XTick = 1:numel(columnLabels);
ax.XTickLabel = columnLabels; 


%% to create marks
A = predMarks_ordered;
% Create the figure
figure;
axis off;

[nRows, nCols] = size(A);

% Loop through the matrix and place symbols
for i = 1:nRows
    for j = 1:nCols
        x = j;
        y = nRows - i + 1;  % Flip y-axis to match matrix visual
        
        % Draw the grid cell (optional for visual)
        rectangle('Position', [x-1, y-1, 1, 1], 'EdgeColor', 'k');
        
        % Add symbols
        if A(i, j) == 1
            text(x - 0.5, y - 0.5, '✔', 'FontSize', 20, 'HorizontalAlignment', 'center');
        elseif A(i, j) == 0
            text(x - 0.5, y - 0.5, '✘', 'FontSize', 20, 'HorizontalAlignment', 'center');
        end
    end
end

xlim([0 nCols])
ylim([0 nRows])

toc
%% defining functions
function lastValues = readTXTfile(pathname4Data,fName,specIDToExtract,specID)
    % Specify the file path
    nfname = fullfile(pathname4Data,fName);
    nfilename = nfname; 

    % Read the data using importdata
    data = importdata(nfilename);


    % Initialize variables to store the last values
    lastValues = struct();

    % Extract the last values of specific factors
    for i = 1:numel(specIDToExtract)
        factorIndex = find(strcmp(specID, specIDToExtract{i}));
        factorData = data.data(factorIndex,:);
        lastValue = factorData(end);
        lastValues.(specIDToExtract{i}) = lastValue;
    end


end

function outPut_indeces = findIndicesInSpecID(output_nodes,specID)
    outPut_indeces = zeros(1, numel(output_nodes));  % Create an array to store indices
    for i = 1:numel(output_nodes)
        node = output_nodes{i};
        index = find(ismember(specID, node));  % Find the index for the current node
        outPut_indeces(i) = index;  % Store the index in the array
    end
end

function perf = calculate_Performance(diff_Mat_mod)
    matrix = ones(size(diff_Mat_mod,1),size(diff_Mat_mod,2))*1000;
    
    matrix(1,:) =[1 1 1 1000 1 1000 1 1 1 1 1 1 1 1 1000 1000  1 1];
    matrix(2,:) =[1 1 1 1000 1 1 1000 1000 1 1 1000 1000 1 1000 1 1000  1000 1];
    matrix(3,:) =[1 1 1 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 -1 1000  1000 1];
    matrix(4,:) =[1 1 1 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1 1000  1000 1000];
    matrix(5,:) =[1 1 1 1 1000 1000 1000 1000 1000 1 1000 1 1000 1000 1 1000  1000 1];
    matrix(6,:) =[1 1 1000 1000 1 1 1000 1000 1000 1000 1000 1 1000 1000 1000 1  1000 1];
    matrix(7,:) =[-1 -1 0 1 0 1 0 1000 1000 1000 1000 1 0 1000 -1 1000  1 1000];
    matrix(8,:) =[-1 -1 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1 1000  1000 1000];
    matrix(9,:) =[-1 1000 1000 1000 1 1000 1 1 1000 1000 1000 1000 1000 1000 1 1000  1 1];
    matrix(10,:) =[-1 1000 -1 1000 1000 1000 1000 1000 1000 1000 1000 1 1000 1000 -1 1000  1000 1000];
    matrix(11,:) =[0 0 0 1000 0 0 1000 1000 0 0 1000 1000 0 1000 0 1000  0 1000];
    % matrix(12,:) =[1000 1000 1000 1000 -1 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000  1000 1000 1000];



    perf = 0;
    sz1 = size(matrix,1)*size(matrix,2);
    sz = sz1 - size(find(matrix==1000),1);
    diff = diff_Mat_mod-matrix;
    sz_True = size(find(diff==0),1);
    perf = sz_True/sz;


end
