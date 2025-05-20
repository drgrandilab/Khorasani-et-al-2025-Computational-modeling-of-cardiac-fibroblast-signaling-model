close all;
clear all;
clc;
%% u need to take care of these values
% parametrs and variables:
path = '/Users/najmekhorasani/Library/CloudStorage/Dropbox/Najme Folder/article/Atrial-fibrosis-Model/Github/Validation_Sensitivity_Analysis_Netflux/';

pathName = [path, 'barplots_for_AF_TRPC_KD_CaMKII_KD_files/'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nodes4activity = {'fibronectin','CI','CIII',...
 'aSMA','periostin'};

fileNames = {'base_120h_angII_72h_base_120h_angII_TRPCKD_72h_T1',...
    'base_120h_angII_72h_base_120h_angII_TRPCKD_72h_T1',...
'base_120h_angII_72h_base_120h_angII_CaMKIIn_72h_T1',...
'base_120h_angII_72h_base_120h_angII_CaMKIIn_72h_T1',...
'base_120h_angII_48h_base_120h_angII_TRPCKD_48h_base_120h_TRPCKD_48h_T1',...
'base_120h_angII_48h_base_120h_angII_TRPCKD_48h_base_120h_TRPCKD_48h_T1',...
'base_120h_angII_48h_base_120h_angII_TRPCKD_48h_base_120h_TRPCKD_48h_T1',...
'base_120h_angII_48h_base_120h_angII_TRPCKD_48h_base_120h_TRPCKD_48h_T1'};

labels = {'AF', 'AF+pyr3', 'AF', 'AF+CaMKIIn', ...
    'CTL', 'AngII', 'AngII+pyr3', 'pyr3'};
tend_ind = [1200+240, 1200+720+1200+240, ...
    1200+720, 1200+720+1200+720, ...
    1200, 1200+240, 1200+480+1200+240, (1200+48)*2+1200+240];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sz1 = size(nodes4activity,2);
sz2 = size(fileNames,2);

load([pathName,'specID.mat']);

Vals = zeros(sz1,sz2);

for i=1:sz2
    name1 = fileNames{i};
    fName2Read1 = [name1,'.txt'];
    wholeValues1 = readTXTfile(fName2Read1,nodes4activity,pathName,specID);
    Vals(:,i) = wholeValues1(:,tend_ind(i));
end


%% plot aSMA
for co=4:4
    figure;
    
    bar(Vals(co,1:2),'FaceColor', 'k');

    ylim([0,1])
    yticks(0:.5:1)
    xticks(1:sz2);                % Set x-tick positions (match number of bars)
    xticklabels(labels(1:2));
    set(gca,'FontSize',20);
    xtickangle(45);
    ylabel('Activity')
    
    title(nodes4activity{co});
    axis square
end


%% plot fibrosis level based on the average value over fibronectin, CI , and CIII
V = 0;
for co=1:3
    V = V + Vals(co,3:4);
end
V = V/3;
figure;

bar(V,'FaceColor', 'k');

ylim([0,1])
yticks(0:.5:1)
xticks(1:sz2);                % Set x-tick positions (match number of bars)
xticklabels(labels(3:4));
set(gca,'FontSize',20);
xtickangle(45);
ylabel('Activity')

title('Fibrosis');
axis square

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
