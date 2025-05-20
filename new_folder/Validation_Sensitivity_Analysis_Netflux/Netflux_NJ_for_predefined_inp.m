% clear all;
% clc;
function [specID,y0] = Netflux_NJ_for_predefined_inp(indices_in_orig_w, ...
            values_for_w_in_netFlux, fname, pathname,It2reachSS)
    global w n K tau ymax y0 specID reactionIDs reactionRules paramList ODElist CNAmodel tstep tUnit tUnitLabel myAxes
    global tEnd tNow tCum yCum tspan options count indPlot
    tEnd = 10;
    tNow = 0;
    tCum = []; yCum = [];
    tspan = [tNow,tEnd];
    options = [];
    count = 0; 
    
    %% Main Code
    % y0 keeps the final values of the nodes
    
    
    %% to load the data from file,
    % the nodes, connections and also initial weights
    openModel(fname,pathname);
    
    %to modify w based on ga's W_Mat
    % 1:AngII, 2:TGFB, 3:mechanical, \n' ...
    % ,'4: IL6, 5: IL1, 6: TNFa, 7: NE, 8: PDGF, 9: ET1, 10: NP, \n' ...
    % '11: Forskolin, 12: Mel, 13: PDE10A, 14: NR4A1
    % missed inputs: 3: mechanical, 5:IL1, 7:NE, 11:Forskolin, 12:Mel => W_Mat[1:5]
    w(indices_in_orig_w) = values_for_w_in_netFlux;
    % w(missed_inp_in_ExcelFile) = W_Mat(missed_inp_in_W);
    % tmp = first_non_inp_in_W;
    % w(1,last_inp_ind_in_exFile+1:end) = W_Mat(1,tmp:end);

    
    %% to run for predefined W
    % calculate network output in baseline
    for i=1:It2reachSS
        indPlot = 1;
        start();
    end
    disp('NetFlux Done!')
end

%% GUI subfunctions
function updateTend()%checked
    global tEnd tNow tspan;
    tspan = [tNow tNow+tEnd];
end
function updateTstep()%checked
    global tstep;
    tstep = .1;
end

function resetParameters()%checked
   global w n K tau ymax y0  paramList
   [w,n,K,tau,ymax,y0] = paramList{:};   
   % disp('Parameters have been reset!');
end  
function start()%checked
    global w n K tau ymax y0 specID reactionIDs reactionRules paramList ODElist CNAmodel tstep tUnit tUnitLabel myAxes
    global tEnd tNow tCum yCum tspan options count
    % if (K.^n < .5)
    %     disp( 'Warning: EC50 and n combinations are negative');
    % end
    runSimulation;
end
function runSimulation()
    global w n K tau ymax y0 specID reactionIDs reactionRules paramList ODElist CNAmodel tstep tUnit tUnitLabel myAxes
    global tEnd tNow tCum yCum tspan options count indPlot
    signal = 1;
    % if isequal(signal,1) %shows start of simulation
    %     disp('Simulation running...');
    % end
    rpar = [w;n;K];         %Reaction Parameters
    params = {rpar,tau,ymax,specID};
    [t,y]=ode23(@util.ODE,tspan,y0,options,params,ODElist);
    y = real(y);

    % disp('Interpolating...');
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
    % if indPlot ~= -1
    %     if indPlot == 0
    %         selectedVars = input('which specID u wanna plot? ');
    %     elseif indPlot == 1
    %             selectedVars = 90;%CI
    %     end
        
        % set(gcf,'DefaultLineLineWidth',1.7); % set the width of the lines in the plot
        % plot(tCum,yCum(:,selectedVars),'.'); 
        % xlabel(sprintf('Time (%s)',tUnit)); ylabel('Activity');
        % legend(specID{selectedVars});
    % end
    signal = 2;
    % if ~isequal(signal,1) 
    %     disp('Simulation Successful!'); 
    % end
%     updateDisplayedSpeciesParams
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
    % plot(0,0); axis([0 1 0 1]); 
    % xlabel(sprintf('Time (%s)',tUnit));ylabel('Activity');
    % disp('Simulation has been reset');
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