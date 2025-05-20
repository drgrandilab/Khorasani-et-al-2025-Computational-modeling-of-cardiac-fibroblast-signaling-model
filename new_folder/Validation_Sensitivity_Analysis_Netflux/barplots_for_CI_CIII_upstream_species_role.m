close all;
clear all;
clc;
%% u need to take care of these values
% parametrs and variables:
path = '/Users/najmekhorasani/Library/CloudStorage/Dropbox/Najme Folder/article/Atrial-fibrosis-Model/Github/Validation_Sensitivity_Analysis_Netflux/';

pathName = [path, 'barplots_for_CI_CIII_upstream_species_role_files/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind1 = 3;
ind2 = 4;
flg = 0;
while flg==0
    flg = input('1: AngII, 2: TGFB: ');
    if flg ==1
        fileNames = {'base_24h_angII_24h_AngII_72h_T1.txt',...
            'base_24h_angII_24h_AngII_CaKD_72h_T1.txt',...
            'base_24h_angII_24h_AngII_CaKD_MMP1KD_72h_T1.txt',...
        'base_24h_angII_24h_AngII_CaKD_proMMP1KD_72h_T1.txt',...
        'base_24h_angII_24h_AngII_CaKD_CImRNAOE_72h_T1.txt',...
        'base_24h_angII_24h_AngII_CaKD_CIIImRNAOE_72h_T1.txt',...
        };
    elseif flg==2
        fileNames = {'base_24h_TGFB_24h_TGFB_72h_T1.txt',...
            'base_24h_TGFB_24h_TGFB_CaKD_72h_T1.txt',...
            'base_24h_TGFB_24h_TGFB_CaKD_MMP1KD_72h_T1.txt',...
        'base_24h_TGFB_24h_TGFB_CaKD_proMMP1KD_72h_T1.txt',...
        'base_24h_TGFB_24h_TGFB_CaKD_CImRNAOE_72h_T1.txt',...
        'base_24h_TGFB_24h_TGFB_CaKD_CIIImRNAOE_72h_T1.txt'};
    else
        disp('wrong input!');
        flg=0;
    end
end

nodes4activity = {'CI','CIII',...
 'CImRNA','CIIImRNA', 'MMP1', 'proMMP1','Ca'};

labels = {'CTL', 'Ca-KD'...
    'Ca-KD-MMP1-KD', 'Ca-KD-proMMP1-KD',...
    'Ca-KD-CImRNA-OE', 'Ca-KD-CIIImRNA-OE'};
tend_ind = [480, 720,720,720,720,720];
colors_ = {[0.5 0.5 0.5],'r-.','g','k--','b:','b:'};
ln_width = [3, 6, 2, 2, 4, 4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sz1 = size(nodes4activity,2);
sz2 = size(fileNames,2);

load([pathName,'specID.mat']);

Vals = zeros(sz1,sz2);

for i=1:sz2
    name1 = fileNames{i};
    fName2Read1 = [name1];
    wholeValues1 = readTXTfile(fName2Read1,nodes4activity,pathName,specID);
    Vals(:,i) = wholeValues1(:,tend_ind(i));
    wholeValues{i} =wholeValues1;
end

%% Dynamics of Collagens and their upstreams species 
%% CIII, CIIImRNA : [2,4]
x=[1:10:1200];
arr1 = [2,4];
arr2 = [1,2,3,4,6];
for ind1=1:2
    i = arr1(ind1);
    figure;
    for ind2=1:5
        co = arr2(ind2);
        traceTMP = wholeValues{co};
        if isnumeric(colors_{co})
            plot(x,traceTMP(i,x), '-', 'Color', colors_{co}, 'LineWidth', ln_width(co))
        else
            plot(x,traceTMP(i,x), colors_{co}, 'LineWidth', ln_width(co))
        end
        hold on;
    end
    ylabel('Activity')
    title(nodes4activity{i})
    ylim([0,1.2]);
    yticks(0:.2:1.2);
    xlim([0,1200]);
    xticks(0:240:1200);  
    xticklabels([0,24,48,72,96,120]);
    xlabel('Time (h)')
    legend(labels(arr2));
    set(gca,'FontSize',20);
end

%% CI, CImRNA, MMP1, proMMP1: [1,3,5,6]
x=[1:10:1200];
arr1 = [1,3,5,6];
arr2 = [1,2,3,4,5];
for ind1=1:4
    i = arr1(ind1);
    figure;
    for ind2=1:5
        co = arr2(ind2);
        traceTMP = wholeValues{co};
        if isnumeric(colors_{co})
            plot(x,traceTMP(i,x), '-', 'Color', colors_{co}, 'LineWidth', ln_width(co))
        else
            plot(x,traceTMP(i,x), colors_{co}, 'LineWidth', ln_width(co))
        end
        hold on;
    end
    ylabel('Activity')
    title(nodes4activity{i})
    ylim([0,1.2]);
    yticks(0:.2:1.2);
    xlim([0,1200]);
    xticks(0:240:1200);  
    xticklabels([0,24,48,72,96,120]);
    xlabel('Time (h)')
    legend(labels);
    set(gca,'FontSize',20);
end

%% barplots of CI and CIII
% Define custom colors
customColors = [
    0.5 0.5 0.5;  % gray
    1   0   0;    % red
    0   1   0;    % green
    0   0   0;    % black
    0   0   1;    % solid blue
    0.6 0.8 1     % blue again (for striped, same base color for now)
];

for co=1:2
    figure;

    b = bar(Vals(co,:), 'FaceColor', 'flat');  % Use 'flat' to enable CData coloring
    
    % Apply custom colors to individual bars
    b.CData = customColors;

    ylim([0,1])
    yticks(0:.5:1)
    xticks(1:sz2);                % Set x-tick positions (match number of bars)
    xticklabels(labels);
    set(gca,'FontSize',20);
    xtickangle(45);
    ylabel('Activity')

    title(nodes4activity{co});
    axis square
end

%% functions
function lastValues = readTXTfile(fName,specIDToExtract,pathName,specID)
    % specIDToExtract = {'CI','CIII','aSMA', 'MMP1', 'MMP2','MMP9'};
    global w n K tau ymax y0 specID reactionIDs reactionRules paramList ODElist CNAmodel tstep tUnit tUnitLabel myAxes
    global tEnd tNow tCum yCum tspan options count indPlot pathname2saveData
    % Specify the file path
    %     default = 'Data Sheet1.txt';
    nfname = fullfile(pathName,fName);
    nfilename = nfname; 

    % Read the data using importdata
    data = importdata(nfilename);
    

    % Initialize variables to store the last values
    lastValues = [];

    % Extract the last values of specific factors
    for i = 1:numel(specIDToExtract)
        factorIndex = find(strcmp(specID, specIDToExtract{i}));
        factorData = data.data(factorIndex,:);
        lastValues(i,:) = factorData;
    end

end
