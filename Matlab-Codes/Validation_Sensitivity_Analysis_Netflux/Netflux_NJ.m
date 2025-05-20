% this code do the KD after reaching SS in base line

clear all;
clc;
close all;
global w n K tau ymax y0 specID reactionIDs reactionRules paramList ODElist CNAmodel tstep tUnit tUnitLabel myAxes
global tEnd tNow tCum yCum tspan options count indPlot pathname2saveData
tEnd = 10;
tNow = 0;
tCum = []; yCum = [];
tspan = [tNow,tEnd];
options = [];
count = 0; 

%% I wanna ask if we need to set an input value high in both CTL and KNO 
% Press "1" here if you want to set two input signals to 1. 
% Your selection at this step determines the first input signal.
flgFixedHigh = input('Do u wanna set one of the inputs equal to 1? Yes: 1 No: o.w.? ');
extented_name = '';
if flgFixedHigh ==1
    indTmp = input(['Which signal input? 1:AngII, 2:TGFB, 3:mechanical, \n' ...
            ,'4: IL6, 5: IL1, 6: TNFa, 7: NE, 8: PDGF, 9: ET1, 10: NP, \n' ...
            '11: Forskolin, 12: Mel, 13: PDE10A, 14: NR4A1, ... ']);
    extented_name = ['high_',num2str(indTmp)];
end

%% Main Code
% the model file you wanna upload
fname = ['Khorasani_Atrial_Fibroblast_Model.xlsx'];

% the folder you wanna save the model dynamics over time that you ran in
% different conditions
folderName = ['dataNJ_',extented_name];

% the path for loading model
path = '/Users/najmekhorasani/Library/CloudStorage/Dropbox/Najme Folder/article/Atrial-fibrosis-Model/Github/Validation_Sensitivity_Analysis_Netflux/';
pathname = [path,'models/'];
% the path for saving result files
pathname2saveData = [path,folderName,'/'];

% open the model
openModel(fname,pathname);

%% to save nodes names
save([pathname2saveData,'specID.mat'],'specID');
% iteration no to reach SS
It2reachSS = 15;

%% to determine what we wanna do with Ctl data
% if it is your first time running this code, you need to save everything.
% therefore, you need to press 1.
flg = input('Do u wanna save data for Ctl condition, Yes: 1 No: o.w.? ');

% here you need to choose if you want all inputs as baseline (press any key), or you need
% to set one of the input signals to 1 (press 1). 
if flg ==1
    % to simulate high expression of TGFB, or high mechanichal
    flgC = input('Do u have any specific signalling context, Yes: 1 No: o.w.? ');
    if flgC ==1
        % which input u wanna consider high in val.
        ind = input(['Which signal input? 1:AngII, 2:TGFB, 3:mechanical, \n' ...
            ,'4: IL6, 5: IL1, 6: TNFa, 7: NE, 8: PDGF, 9: ET1, 10: NP, \n' ...
            '11: Forskolin, 12: Mel, 13: PDE10A, 14: NR4A1, ... ']);
        w(ind)=.9;
        
        if ind ==1
            fName2w = 'Data Sheet_Ctl_highAng.txt';
        elseif ind==2
            fName2w = 'Data Sheet_Ctl_highTGFB.txt';
        elseif ind==3
            fName2w = 'Data Sheet_Ctl_highMech.txt';
        elseif ind==4
            fName2w = 'Data Sheet_Ctl_highIL6.txt';
        elseif ind==5
            fName2w = 'Data Sheet_Ctl_highIL1.txt';
        elseif ind==6
            fName2w = 'Data Sheet_Ctl_highTNFa.txt';
        elseif ind==7
            fName2w = 'Data Sheet_Ctl_highNE.txt';
        elseif ind==8
            fName2w = 'Data Sheet_Ctl_highPDGF.txt';
        elseif ind==9
            fName2w = 'Data Sheet_Ctl_highET1.txt';
        elseif ind==10
            fName2w = 'Data Sheet_Ctl_highNP.txt';
        elseif ind==11
            fName2w = 'Data Sheet_Ctl_highForskolin.txt';
        elseif ind==12
            fName2w = 'Data Sheet_Ctl_highMel.txt';
        elseif ind==13
            fName2w = 'Data Sheet_Ctl_highPDE10A.txt';
        elseif ind==14
            fName2w = 'Data Sheet_Ctl_highNR4A1.txt';
        end
    else
            fName2w = 'Data Sheet_Ctl.txt';
    end

    %% fixed high input
    if flgFixedHigh ==1
        w(indTmp) = .9;
    end

    % calculate network output in baseline
    for i=1:It2reachSS
        indPlot = 1;
        start();
    end

    % save the Ctl data
    dataSheet(fName2w);
end

%% to determine what we wanna do with KNO data
flg = input('Do u wanna save data for KD condition, Yes: 1 No: o.w.? ');
if flg ==1
    % to simulate high expression of TGFB, or high mechanichal
    flgC = input('Do u have any specific signalling context, Yes: 1 No: o.w.? ');
    if flgC ==1
        % which input u wanna consider high in val.
        ind = input(['Which signal input? 1:AngII, 2:TGFB, 3:mechanical, \n' ...
            ,'4: IL6, 5: IL1, 6: TNFa, 7: NE, 8: PDGF, 9: ET1, 10: NP, \n' ...
            '11: Forskolin, 12: Mel, 13: PDE10A, 14: NR4A1, ... ']);
        w(ind)=.9;
        
        if ind ==1
            fName2wEND = '_highAng.txt';
        elseif ind==2
            fName2wEND = '_highTGFB.txt';
        elseif ind==3
            fName2wEND = '_highMech.txt';
        elseif ind==4
            fName2wEND = '_highIL6.txt';
        elseif ind==5
            fName2wEND = '_highIL1.txt';
        elseif ind==6
            fName2wEND = '_highTNFa.txt';
        elseif ind==7
            fName2wEND = '_highNE.txt';
        elseif ind==8
            fName2wEND = '_highPDGF.txt';
        elseif ind==9
            fName2wEND = '_highET1.txt';
        elseif ind==10
            fName2wEND = '_highNP.txt';
        elseif ind==11
            fName2wEND = '_highForskolin.txt';
        elseif ind==12
            fName2wEND = '_highMel.txt';
        elseif ind==13
            fName2wEND = '_highPDE10A.txt';
        elseif ind==14
            fName2wEND = '_highNR4A1.txt';
        end
    else
            fName2wEND = '.txt';
    end
    
    % determine specIDs to KNO
    specIDToKNO = specID;
    for i = 1:numel(specIDToKNO)
        factorIndex = find(strcmp(specID, specIDToKNO{i}));
        % we want to plot only the last iteration
        indPlot = -1;
        % reset the sys. to start over
        resetParameters();
        resetSimulation();
        
        %% read the Ctl data
        fName2Read = ['Data Sheet_Ctl',fName2wEND];
        specIDToExtract = specID;
        lastValuesCtl = readTXTfile(fName2Read,specIDToExtract);
        valuesCtl = structfun(@(x) x, lastValuesCtl);
        
        % to set yo equal to yend of Ctl
        y0 = real(valuesCtl);

        %  to KNO
        ymax(factorIndex) = ymax(factorIndex) * .1;  
        
        % to simulate high expression of TGFB
        if flgC ==1
            w(ind)=.9;
        end

        %% fixed high input
        if flgFixedHigh ==1
            w(indTmp) = .9;
        end

        for j=1:It2reachSS
            % to plot the last trace
            if j == It2reachSS
                indPlot = 1;
            end

            start();
        end

        %to save the values
        fName2w = ['Data Sheet KNO_',specIDToKNO{i},fName2wEND];
        dataSheet(fName2w);  
    end
end



%% read the data
% to simulate high expression of TGFB, or high mechanichal
flgC = input('Do u have any specific signalling context, Yes: 1 No: o.w.? ');
if flgC ==1
    % which input u wanna consider high in val.
    ind = input('Which signal input? 1:AngII, 2:TGFB, 3:mechanical, 12:Mel... ');
    
     if ind ==1
            fName2wEND = '_highAng.txt';
        elseif ind==2
            fName2wEND = '_highTGFB.txt';
        elseif ind==3
            fName2wEND = '_highMech.txt';
        elseif ind==4
            fName2wEND = '_highIL6.txt';
        elseif ind==5
            fName2wEND = '_highIL1.txt';
        elseif ind==6
            fName2wEND = '_highTNFa.txt';
        elseif ind==7
            fName2wEND = '_highNE.txt';
        elseif ind==8
            fName2wEND = '_highPDGF.txt';
        elseif ind==9
            fName2wEND = '_highET1.txt';
        elseif ind==10
            fName2wEND = '_highNP.txt';
        elseif ind==11
            fName2wEND = '_highForskolin.txt';
        elseif ind==12
            fName2wEND = '_highMel.txt';
        elseif ind==13
            fName2wEND = '_highPDE10A.txt';
        elseif ind==14
            fName2wEND = '_highNR4A1.txt';
     end
else
        fName2wEND = '.txt';
end

%% read the Ctl data
fName2Read = ['Data Sheet_Ctl',fName2wEND];
specIDToExtract = specID;
lastValuesCtl = readTXTfile(fName2Read,specIDToExtract);

%% to read data after KNO
% to determine KNO-specIDs to read from their data
KNOspecIDToRead = specID;

for i = 1:numel(KNOspecIDToRead)
    fName2Read = ['Data Sheet KNO_',KNOspecIDToRead{i},fName2wEND];
    lastValues{i} = readTXTfile(fName2Read,specIDToExtract);
end

%% comparing Ctl versus KNO
valuesCtl = structfun(@(x) x, lastValuesCtl);
for i = 1:numel(KNOspecIDToRead)
    values(:,i) = structfun(@(x) x, lastValues{i});
    deltaActivity(:,i) =  values(:,i) - valuesCtl;

    % Replace zero values in valuesCtl with a small non-zero value
    eps = 0;%0.0000001;
    valuesCtl(valuesCtl == 0) = eps;

    deltaActivityR(:,i) = deltaActivity(:,i)./valuesCtl;
end

%% Define the labels and plot the heatmap
yLabels = KNOspecIDToRead;
columnLabels = specIDToExtract;

% to create mymap
no = 100;
stp = 1/no;
c1 = 0:stp:1;
c2 = 0.99:-1*stp:0;
mymap = [[c1';ones(no,1)] [c1';c2'] [ones(no+1,1);c2'] ];

% Define the range of possible values in deltaActivity
val1 = -.25;  % Lower bound
val2 = .25;   % Upper bound
stp = (val2-val1)*stp/2;
potInterval = val1:stp:val2;

%% to calculate sensitivity and influemce of each node through all
% knockdowns
number_of_top_vals = 10;
[sensitivityArr,influenceArr,...
    top_n_max_sens_val,top_n_max_infl_val,...
    indicesSens, indicesInfl,...
    top_n_max_sens_nodes, top_n_max_infl_nodes] = sens_inf_values(deltaActivity, specID,number_of_top_vals);

figure;
[sval,sind] = sort(top_n_max_sens_val);
hbh = barh(sval,'FaceColor',[.95 .95 .95],'EdgeColor',[0 0 0],'LineWidth',2.5);
axis tight;
ln = length(top_n_max_sens_nodes);
lbl = top_n_max_sens_nodes(sind);
xlim([0,ceil(max(sval))])
set(gca, 'YTick', [1:ln], 'YTickLabel', lbl)%,'XTickLabelRotation', 90);
set(gca, 'XTick',[0:ceil(max(sval)):ceil(max(sval))]);
set(gca,'fontsize',30);

title('Bar Plot of sensitivity values.');
exportgraphics(gcf, [pathname2saveData,'sen_barPlot_',fName2wEND,'.pdf']);

figure;
[sval,sind] = sort(top_n_max_infl_val);
hbh = barh(sval,'FaceColor',[.95 .95 .95],'EdgeColor',[0 0 0],'LineWidth',2.5);
axis tight;
ln = length(top_n_max_infl_nodes);
lbl = top_n_max_infl_nodes(sind);
xlim([0,ceil(max(sval))])
set(gca, 'YTick', [1:ln], 'YTickLabel', lbl);%, 'XTickLabelRotation', 90);
set(gca, 'XTick',[0:ceil(max(sval)):ceil(max(sval))]);
set(gca,'fontsize',30);

% Set the title
title('Bar Plot of influence values.');
exportgraphics(gcf, [pathname2saveData,'infl_barPlot_',fName2wEND,'.pdf']);


%% to Modify deltaActivity
deltaActivity_mod = deltaActivity;

deltaActivity_mod(deltaActivity_mod>=val2) = val2;
deltaActivity_mod(abs(deltaActivity_mod)<.001) = 0;
deltaActivity_mod(deltaActivity_mod<=val1) = val1;

%% prepare data to be visualized
% define range of data in deltaActivity
mn = min(min(deltaActivity_mod));
mx = max(max(deltaActivity_mod));

%scale mn and mx to  two-floor-digit values: (bec. no =100)
mn = floor(mn*no)/no;
mx = floor(mx*no)/no;
potInterval = floor(potInterval*no)/no;

% Determine the number of colors to select from mymap
ind1 = find(potInterval==mn);
ind2 = find(potInterval==mx);

% DEfine a new mymap
figure;
mymap2 = mymap(ind1:ind2,:);

h= heatmap(yLabels, columnLabels, deltaActivity_mod);
h.GridVisible = 'off';
colormap(mymap2);
ax = gca;
ax.FontSize = 05; 
% saveas(gcf, [pathname2saveData,'sensAna',fName2wEND,'.pdf']);
exportgraphics(gcf, [pathname2saveData,'sensAna',fName2wEND,'.pdf']);


%% define my functions
function [sensitivityArr,influenceArr,...
    top_n_max_sens_val,top_n_max_infl_val,...
    indicesSens, indicesInfl,...
    top_n_max_sens_nodes, top_n_max_infl_nodes] = sens_inf_values(deltaActivity, specID,number_of_top_vals)

    sensitivityArr = sum(abs(deltaActivity),2);
    sensitivityArr = sensitivityArr';
    influenceArr = sum(abs(deltaActivity),1);

    [top_n_max_sens_val, indicesSens] = maxk(sensitivityArr, number_of_top_vals);
    [top_n_max_infl_val, indicesInfl] = maxk(influenceArr, number_of_top_vals);
    top_n_max_sens_nodes = specID(indicesSens);
    top_n_max_infl_nodes = specID(indicesInfl);

end

%%%%%%%%%%%%%%%%%%%%%%%
%% end of part 1.
%%%%%%%%%%%%%%%%%%%%%%%

%% GUI subfunctions
function updateTend()%checked
    global tEnd tNow tspan;
    tspan = [tNow tNow+tEnd];
end
function updateTstep()%checked
    global tstep;
    tstep = .1;
end
function updateTunit()
    tUnit = tUnitLabel.Value;
    tEnd.Label = sprintf('Simulation time (%s):',tUnit);
    tstepLabel.Label = sprintf('Time step (%s):',tUnit);
    xlabel(sprintf('Time (%s)', tUnit))
    %insert a test 
end

function resetParameters()%checked
   global w n K tau ymax y0  paramList
   [w,n,K,tau,ymax,y0] = paramList{:};   
   disp('Parameters have been reset!');
end  
function start()%checked
    global w n K tau ymax y0 specID reactionIDs reactionRules paramList ODElist CNAmodel tstep tUnit tUnitLabel myAxes
    global tEnd tNow tCum yCum tspan options count
    if (K.^n < .5)
        disp( 'Warning: EC50 and n combinations are negative');
    end
    runSimulation;
end
function runSimulation()
    global w n K tau ymax y0 specID reactionIDs reactionRules paramList ODElist CNAmodel tstep tUnit tUnitLabel myAxes
    global tEnd tNow tCum yCum tspan options count indPlot
    signal = 1;
    if isequal(signal,1) %shows start of simulation
        disp('Simulation running...');
    end
    rpar = [w;n;K];         %Reaction Parameters
    params = {rpar,tau,ymax,specID};
    [t,y]=ode23(@util.ODE,tspan,y0,options,params,ODElist);
    y = real(y);

    disp('Interpolating...');
    % interpolate to equally spaced time points
    x = (tNow:tstep:tspan(end))';
    yi = interp1(t,y,x,'pchip');
    t = x;
    y = yi;

    tNow = t(end);
    updateTend;
    y0 = y(end,:);         
    if count < 1
         tCum = [tCum;t];
         yCum = [yCum;y];
    else
        tCum = [tCum; t(2:length(t))];   % I want to remove the 1st entry
        sizeY = size(y);
        sizeY = sizeY(1);
        yCum = [yCum; y(2:sizeY,:)]; % then concat to the end of the cumulative 
    end
    count = count + 1; 
    % plot results
    if indPlot ~= -1
        if indPlot == 0
            selectedVars = input('which specID u wanna plot? ');
        elseif indPlot == 1
                selectedVars = 90;%CI
        end
        
        set(gcf,'DefaultLineLineWidth',1.7); % set the width of the lines in the plot
        plot(tCum,yCum(:,selectedVars),'.'); 
        xlabel(sprintf('Time (%s)',tUnit)); ylabel('Activity');
        legend(specID{selectedVars});
    end
    signal = 2;
    if ~isequal(signal,1) 
        disp('Simulation Successful!'); 
    end
%     updateDisplayedSpeciesParams
end
function plotData(hWidget)
    selectedVars = ismember(specID,speciesListbox.Value);
    if(isempty(yCum))
        statusLabel.MenuItems = 'No Plot: please ''Simulate'' first'; 
        statusLabel.Value = [];
    else
        plot(tCum,yCum(:,selectedVars));        
        xlabel(sprintf('Time (%s)',tUnit)); ylabel('Activity');
        legend(speciesListbox.Value);
        statusLabel.MenuItems = 'Plot Successful!'; 
        statusLabel.Value = [];
    end
end
function resetSimulation()%checked
    global w n K tau ymax y0 specID reactionIDs reactionRules paramList ODElist CNAmodel tstep tUnit tUnitLabel myAxes
    global tEnd tNow tCum yCum tspan options count
    tNow = 0;
    updateTend;   
    tCum = [];
    yCum = [];
    count = 0; 
    y0 = paramList{6};
    plot(0,0); axis([0 1 0 1]); 
    xlabel(sprintf('Time (%s)',tUnit));ylabel('Activity');
    disp('Simulation has been reset');
end


function openModel(fname,pathname)
    global w n K tau ymax y0 specID reactionIDs  paramList ODElist CNAmodel tstep tUnit tUnitLabel myAxes
    try
        if ~isequal(fname,0)
            xlsfilename = [pathname fname];

            [specID,reactionIDs,reactionRules,paramList,ODElist,CNAmodel, error] = util.xls2Netflux(strrep(fname,'.xls',''),xlsfilename);
            [w,n,K,tau,ymax,y0] = paramList{:};

            tauLabel.ValueRange = [1e-9 1e9];
            yLabel.ValueRange = [0 10];
            y0Label.ValueRange = [0 10];
            nLabel.ValueRange = [.5 5];
            kLabel.ValueRange = [0 5];

            errsignal = false;
            if ~isempty(error{1}) || ~isempty(error{2}) %check to see if error passed from xls2Netflux exists
                statusLabel.MenuItems = error; % display the errors in the status window
                statusLabel.Value = [];
                errsignal = true;
                updateDisplayedSpeciesParams;
                updateDisplayedReactionParams;
                updateTstep;
            end

            if ~errsignal
                resetParameters;
                resetSimulation;
                updateTstep;
            end
            errsignal = false;
        else
            disp('Excel import canceled');
        end

    catch e
        errorstr = {'Error reading file, try resaving as new .xlsx',...
            e.identifier, e.message,'line:', e.stack.name,num2str(e.stack.line)};
        disp( errorstr);
    end
        
end

function dataSheet(fName)
    global w n K tau ymax y0 specID reactionIDs reactionRules paramList ODElist CNAmodel tstep tUnit tUnitLabel myAxes
    global tEnd tNow tCum yCum tspan options count indPlot pathname2saveData
    %check to see if simulation has run
    if (isempty(tCum))
        disp('Need to simulate the model first');  
    else        
        %prompt for filename
%         default = 'Data Sheet1.txt';
        nfname = fullfile(pathname2saveData,fName);
        nfilename = nfname; 

        % create labeled time values
        timeLabeli = tCum;                   
        timeLabeli = cellfun(@num2str,num2cell(timeLabeli),'UniformOutput',false);
        tHeader = {'T'};
        timeLabel = strcat(tHeader,timeLabeli);
        header1 = {'Species_Name'};
        write = cat(2,header1,timeLabel');  %concatenates "species name" and time points        
        tDCum = tCum;  %D stands for Data
        yDCum = yCum;         
        for j = 1:length(specID)                               
            yi = real(yCum); % convert activity values to real numbers           
            %yi = round(yi*10^4)./10^4; %rounds the output to 4 decimal places
            %yi = num2cell(yi);
            label = specID(j);
            for iCol = 1:length(yi(:,j))
                label = horzcat(label,sprintf('%0.4f',yi(iCol,j)')); % concats species header to the y values from t0 to tend
            end
            write = cat(1,write,label);
        end
        util.textwrite(nfilename,write);
        a = 'Data written to ';
        z = horzcat(a,strrep(nfilename,pathname2saveData,''));
        statusLabel.MenuItems = z;
        statusLabel.Value = [];
    end
end


function lastValues = readTXTfile(fName,specIDToExtract)
    % specIDToExtract = {'CI','CIII','aSMA', 'MMP1', 'MMP2','MMP9'};
    global w n K tau ymax y0 specID reactionIDs reactionRules paramList ODElist CNAmodel tstep tUnit tUnitLabel myAxes
    global tEnd tNow tCum yCum tspan options count indPlot pathname2saveData
    % Specify the file path
    %     default = 'Data Sheet1.txt';
    nfname = fullfile(pathname2saveData,fName);
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