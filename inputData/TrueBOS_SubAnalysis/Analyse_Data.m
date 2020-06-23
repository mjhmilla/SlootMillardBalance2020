% Functional BOS estimation
% 
% Based on Pilot BOS data
% 
% LST, 16-06-2020

addpath(genpath('C:\Users\Lizeth\Desktop\HeiAGE\PROJECTS_HEIAGE\_GENERAL ANALYSIS\1 Matlab\btk-0.3.0_Win7_MatlabR2009b_64bit'))
addpath(genpath('C:\Users\Lizeth\Desktop\HeiAGE\PROJECTS_HEIAGE\2b. LEA_Gait_STS\- STS\Analysis\3 Matlab Analysis'))

warning off
% clear all, clc

DataPath = ('C:\Users\Lizeth\Desktop\HeiAGE\PROJECTS_HEIAGE\2b. LEA_Gait_STS\- STS\Data\Pilot_BOS\');

ThresholdForce = 50;    % Force Threshold: remove low values (leaning to the other side, more inaccurate FP readings)

% Pilot1
for k=1
DataPathSub = 'Pilot1\';

clear Trials
Trials.Bare_stnd_nar = [1:3];      % narrow, barefoot, stance (both feet)
Trials.Bare_stnd_med = [6];      % medium : not included 7 as feet are moved between 6 and 7 trials
Trials.Bare_stnd_wid = [4 5];      % wide
Trials.Bare_stnd_L   = [8:10];     % one foot (L or R)
Trials.Bare_sts_med  = [11];        % not included 12 as feet are moved between trials
Trials.Shoe_stnd_nar = [15 16];    % with shoe conditions
Trials.Shoe_stnd_med = [17];
Trials.Shoe_stnd_wid = [13 14];
Trials.Shoe_stnd1_L   = [18];           % moved feet in between trials
Trials.Shoe_stnd2_L   = [19];
Trials.Shoe_sts_med  = [20];
end

% Pilot2
for k=1
DataPathSub = 'Pilot2\';

clear Trials
Trials.Bare_stnd_nar = [20];      % narrow, barefoot, stance (both feet) - no trial 21
Trials.Bare_stnd_med = [18 19];      % medium
Trials.Bare_stnd_wid = [16 17];      % wide
Trials.Bare_stnd_L   = [22 23];     % Non-dominant
% Trials.Bare_stnd_R   = [];          % Dominant
Trials.Bare_sts_nar  = [13];        % note: some lower GRF data when sitting!
Trials.Bare_sts_med  = [14];
Trials.Bare_sts_wide  = [15];
Trials.Shoe_stnd_nar = [8 9];         % with shoe conditions
Trials.Shoe_stnd_med = [7];         % 6 is missing dile
Trials.Shoe_stnd_wid = [4 5];
Trials.Shoe_stnd_L   = [10 11];
Trials.Shoe_stnd_R   = [12];
Trials.Shoe_sts_nar  = [1];
Trials.Shoe_sts_med  = [2];
Trials.Shoe_sts_wid  = [3];
Trials.Bare_fall_forw = [35:39];          % Falling forward
Trials.Bare_fall_backw = [46:48];         % Falling backward (trust fall)
Trials.Bare_stnd_toe = [24];           % Standing Toes
Trials.Bare_stnd_toe_L = [27];           % Standing Toes one foot
end


%% For-loop Analysis
% ConsNames = {'Shoe_stnd_nar'};
ConsNames = fieldnames(Trials);

PlotOn = 1;

for iCon=1:size(ConsNames,1)
    
    %% Readin data
    for k=1
        COP=[]; FOR=[]; MOM=[]; MarkerData=[];
        clear COP_filt FOR_filt MOM_filt
        
        for iTrial = 1:length( Trials.(ConsNames{iCon}) )
            
            %%% C3D from Qualisys: (some reason without COP data)
            [POINTdat,VideoFrameRate,ANALOGdat,AnalogFrameRate,~,ParameterGroup,~,~] = ...
                readC3D_mhs([DataPath,DataPathSub,'MAT\','test00',...
                sprintf('%.2d',Trials.(ConsNames{iCon})(iTrial)),'_done.c3d']);
            
            % POINTdat = 3dim x markers x time     --> Dim: (1) ML - (2) AP - (3) VT   %
            % ANALOGdat (chan x time)
            MarkerLabels = ParameterGroup(1).Parameter(6).data;
            
            % ForceLabels = ParameterGroup(2).Parameter(2).data;
            
            %%% MAT file from Qualisys: (incorrect COP data) --> OLD
            for k=1
                % Filename = 'test0001_done';
                % DataF = load([Filename,'.mat']);
                % DataF = DataF.(Filename);
                % MarkerLabels = DataF.Trajectories.Labeled.Labels;
                % MarkerData = DataF.Trajectories.Labeled.Data/1000;          % [markers x 4dim x time] -> (1)ML - (2)AP - (3)VT
                % COPData1 = DataF.Force(1).COP./1000;                      % [dim x time] --> (1)AP - (2)ML - (3)VT
                % COPData2 = DataF.Force(2).COP./1000;
                % FP2 = DataF.Force(2).Force;                               % [dim x time] --> (1)AP - (2)ML - (3)VT
                % MOM2 = DataF.Force(2).Moment;
                % ax = -MOM2(2,:) ./ FP2(3,:) ;    % ax = -My / Fz
                % ay =  MOM2(1,:) ./ FP2(3,:) ;    % ay = Mx / Fz
            end
            
            %%% TXT file from export V3D (see special pipeline in data folder)
            % COP from V3D because Qualisys exports something weird with offset (also
            % if you calculate from the forces and moments)
            Temp = dlmread([DataPath,DataPathSub,'MAT\','test00',...
                sprintf('%.2d',Trials.(ConsNames{iCon})(iTrial)),'_done_Forces.txt'],'\t',5);
            
            % Throw away data if one foot trial
            for iMark = 1:size(MarkerLabels,2)
                if strcmp(MarkerLabels{iMark}(1),'L')
                    MarkerSide(iMark)=1;
                elseif strcmp(MarkerLabels{iMark}(1),'R')
                    MarkerSide(iMark)=2;
                else
                    MarkerSide(iMark)=0;
                end
            end
            
            if strcmp( ConsNames{iCon}(end-1:end) ,'_L')
                POINTdat(:, find(MarkerSide==2), :) = NaN;
            elseif strcmp( ConsNames{iCon}(end-1:end) ,'_R')
                POINTdat(:, find(MarkerSide==1), :) = NaN;
            end
            
            if strcmp( ConsNames{iCon}(end-1:end) ,'_L') || strcmp( ConsNames{iCon}(end-1:end) ,'_R')
                if nanmean( Temp(:,13) )>350
                    Temp(:,[2:4 8:10 14:16])=NaN; 
                else
                    Temp(:,[5:7 11:13 17:19])=NaN; 
                end
            end
            
            % COP: time, COP11, COP1y, COP1z, COP2x, COP2y, COP2z
            COP = [COP; Temp(:,2:7)];
            FOR = [FOR; Temp(:,8:13)];
            MOM = [MOM; Temp(:,14:19)];
            MarkerData = cat(3, MarkerData, (POINTdat./ 1000)) ;
            clear Temp
            
        end, clear iTrial
        
    end % k
    
    %% Filter Data
    for k=1   
    FS_ana = 900;
    FS_mar = 150;
    FS_fac = FS_ana / FS_mar;
    
    FS_filt = 10;
    [b, a] = butter(4, FS_filt/(FS_ana/2) );
    
    % Filter COP:
    COP_filt = COP;  FOR_filt = FOR; MOM_filt = MOM;
    for ichan=1:6
        Ind = find(~isnan(COP(:,ichan)));
        if ~isempty(Ind)
        COP_filt(Ind,ichan) = filtfilt(b,a,COP(Ind,ichan));
        FOR_filt(Ind,ichan) = filtfilt(b,a,FOR(Ind,ichan));
        MOM_filt(Ind,ichan) = filtfilt(b,a,MOM(Ind,ichan));
        end
    end
    % figure, hold on, plot(COP), plot(COP_filt,'--')
    % figure, hold on, plot(FOR), plot(FOR_filt,'--')
    Weight = nanmean(FOR_filt(:,6)+FOR_filt(:,3));
    
    % Filter Kinem:
    [b, a] = butter(4, FS_filt/(FS_mar/2) );
    Marker_filt = MarkerData;
    for ichan=1:size(MarkerData,2)
        Ind = find(~isnan(MarkerData(1,ichan,:)));
        for idim=1:3
            Marker_filt(idim,ichan,Ind) = filtfilt(b,a, squeeze(MarkerData(idim,ichan,Ind)) );
        end
    end
    % figure, hold on, plot(squeeze(MarkerData(1,:,:))'), plot(squeeze(Marker_filt(1,:,:))','--')
    clear a b FS_filt
    
    end
    
    %% Threshold Force data (vert)
    Ind = find(FOR_filt(:,3) < ThresholdForce);
    FOR_filt(Ind,1:3)=NaN;  MOM_filt(Ind,1:3)=NaN;  COP_filt(Ind,1:3)=NaN;
    
    Ind = find(FOR_filt(:,6) < ThresholdForce);
    FOR_filt(Ind,4:6)=NaN;  MOM_filt(Ind,4:6)=NaN;  COP_filt(Ind,4:6)=NaN;
    clear Temp Ind
    
    %% Plot 2D & Save Fig:
    for k=1
        if PlotOn==1
            figure, hold on, set(gcf,'color',[1 1 1],'position',[20 500 1000 500])
            ConNamePlot = ConsNames{iCon}; ConNamePlot( regexp(ConsNames{iCon},'_') )=' ';  title( ConNamePlot ), clear ConNamePlot
            xlabel('ML (m)'), xlabel('AP (m)')
            plot(COP_filt(:,1),COP_filt(:,2),'k.')
            plot(COP_filt(:,4),COP_filt(:,5),'k.')
            plot(squeeze(Marker_filt(1,:,:))',squeeze(Marker_filt(2,:,:))','.c')
            plot( nanmean(squeeze(Marker_filt(1,[1 3 2 4 5 6 1],:)),2), nanmean(squeeze(Marker_filt(2,[1 3 2 4 5 6 1],:)),2) ,'-b','linewidth',2)
            plot( nanmean(squeeze(Marker_filt(1,[7 9 8 10 11 12 7],:)),2), nanmean(squeeze(Marker_filt(2,[7 9 8 10 11 12 7],:)),2) ,'-b','linewidth',2)
            xlim([ 0.25 0.95]), ylim([0 0.35])
            
            figname = [DataPath,DataPathSub,'Figs\',ConsNames{iCon}];
            saveas(gcf,figname,'fig'), saveas(gcf,figname,'png')
            
            % ==> Could also plot rotated data?!
        end
    end
    
    %% Calculate Limits
    for k=1
        Heel_L = squeeze(Marker_filt([1:2], find(strcmp(MarkerLabels,'L_FCC')==1), :));
        Toe_L =  squeeze(Marker_filt([1:2], find(strcmp(MarkerLabels,'L_FM2')==1), :));
        Heel_R = squeeze(Marker_filt([1:2], find(strcmp(MarkerLabels,'R_FCC')==1), :));
        Toe_R =  squeeze(Marker_filt([1:2], find(strcmp(MarkerLabels,'R_FM2')==1), :));
        
        COP_filt2 =   COP_filt(1:FS_fac:end,:);  COP_rot=[]; Marker_rot = [];
        
        POSrel_R_MLheel = []; POSrel_R_MLtoe = []; POSrel_R_AP = [];
        POSrel_L_MLheel = []; POSrel_L_MLtoe = []; POSrel_L_AP = [];
        
        for iT = 1:length(COP_filt2)
            % FP1 = second FP (under Right Foot)
            
            %%% Left: AP
            % get rotation angle:
            u = [Toe_L(:,iT) - Heel_L(:,iT); 0]';   % u = [nanmean(Toe_L,2) - nanmean(Heel_L,2); 0]';
            v = [0, u(2), 0];
            t = atan2 ( norm(cross(u,v)) , dot(u,v) ); %  deg: t*(180/pi)  % dot: |a| |b| cos(t)
            Rm = [cos(t), -sin(t); sin(t), cos(t)];
            COP_rot(iT,[4 5]) = ( (COP_filt2(iT,[4 5]) - Heel_L(:,iT)') *Rm) + Heel_L(:,iT)';
            MarkerL = find(MarkerSide==1);
            for iM=1:6
                Marker_rot(1:2,MarkerL(iM),iT) = ((squeeze(Marker_filt(1:2,MarkerL(iM),iT)) - ...
                    Heel_L(:,iT))' * Rm) + Heel_L(:,iT)';
            end
            
            % AP position:
            DistAP_L = sqrt( (Toe_L(1,iT)-Heel_L(1,iT)).^2 + (Toe_L(2,iT)-Heel_L(2,iT)).^2 );
            POSrel_L_AP(iT) = (COP_rot(iT,[5]) - Heel_L(2,iT)') ./ DistAP_L *100;      % [ml ap]
            
            % ML position:
            Ank_Lmed = squeeze(Marker_rot([1:2], find(strcmp(MarkerLabels,'L_TAM')==1), :));
            Ank_Llat = squeeze(Marker_rot([1:2], find(strcmp(MarkerLabels,'L_FAL')==1), :));
            Toe_Lmed = squeeze(Marker_rot([1:2], find(strcmp(MarkerLabels,'L_FM1')==1), :));
            Toe_Llat = squeeze(Marker_rot([1:2], find(strcmp(MarkerLabels,'L_FM5')==1), :));

            MidFootAP = (nanmean([Toe_Lmed(2,iT),Toe_Llat(2,iT)]) - nanmean([Ank_Lmed(2,iT),Ank_Llat(2,iT)]))/2 + ...
                nanmean([Ank_Lmed(2,iT),Ank_Llat(2,iT)]);
            if COP_rot(iT,5)<MidFootAP  % lower part: compare to ML ankles
                POSrel_L_MLtoe(iT)  = NaN;
                POSrel_L_MLheel(iT) = 100 - ( COP_rot(iT,4)-Ank_Llat(1,iT)) ./ ...
                    (Ank_Lmed(1,iT)-Ank_Llat(1,iT)) *100;
            else                        % upper part: compare to ML MTP (toes)
                POSrel_L_MLtoe(iT)  = 100 - ( COP_rot(iT,4)-Toe_Llat(1,iT)) ./ ...
                    (Toe_Lmed(1,iT)-Toe_Llat(1,iT)) *100;
                POSrel_L_MLheel(iT) = NaN;
            end, clear MidFootAP Ank_Lmed Ank_Llat Toe_Lmed Toe_Llat
            
            %%% Right: AP
            % get rotation angle:
            u = [Toe_R(:,iT) - Heel_R(:,iT); 0]';   % u = [nanmean(Toe_L,2) - nanmean(Heel_L,2); 0]';
            v = [0, u(2), 0];
            t = atan2 ( norm(cross(u,v)) , dot(u,v) ); %  deg: t*(180/pi)  % dot: |a| |b| cos(t)
            Rm = [cos(-t), -sin(-t); sin(-t), cos(-t)];
            COP_rot(iT,[1 2]) = ( (COP_filt2(iT,[1 2]) - Heel_R(:,iT)') *Rm) + Heel_R(:,iT)';
            MarkerR = find(MarkerSide==2);
            if strcmp(DataPathSub, 'Pilot2\') && strcmp(ConsNames{iCon},'Bare_fall_backw') && iT==1
               MarkerR(7)=[]; for iTemp = 8:13, MarkerLabels{iTemp-1} = MarkerLabels{iTemp}; end, MarkerLabels{13} = [];
            end
            for iM=1:6
                Marker_rot(1:2,MarkerR(iM),iT) = ((squeeze(Marker_filt(1:2,MarkerR(iM),iT)) - ...
                    Heel_R(:,iT))' * Rm) + Heel_R(:,iT)';
            end
            
            % AP position:
            DistAP_R = sqrt( (Toe_R(1,iT)-Heel_R(1,iT)).^2 + (Toe_R(2,iT)-Heel_R(2,iT)).^2 );
            POSrel_R_AP(iT) = (COP_rot(iT,[2]) - Heel_R(2,iT)') ./ DistAP_R *100;      % [ml ap]
            
            % ML position:
            Ank_Rmed = squeeze(Marker_rot([1:2], find(strcmp(MarkerLabels,'R_TAM')==1), :));
            Ank_Rlat = squeeze(Marker_rot([1:2], find(strcmp(MarkerLabels,'R_FAL')==1), :));
            Toe_Rmed = squeeze(Marker_rot([1:2], find(strcmp(MarkerLabels,'R_FM1')==1), :));
            Toe_Rlat = squeeze(Marker_rot([1:2], find(strcmp(MarkerLabels,'R_FM5')==1), :));
            
            MidFootAP = (nanmean([Toe_Rmed(2,iT),Toe_Rlat(2,iT)]) - nanmean([Ank_Rmed(2,iT),Ank_Rlat(2,iT)]))/2 + ...
                nanmean([Ank_Rmed(2,iT),Ank_Rlat(2,iT)]);
            if COP_rot(iT,2)<MidFootAP  % lower part: compare to ML ankles
                POSrel_R_MLtoe(iT)  = NaN;
                POSrel_R_MLheel(iT) = ( COP_rot(iT,1)-Ank_Rmed(1,iT)) ./ ...
                    (Ank_Rlat(1,iT)-Ank_Rmed(1,iT)) *100;
            else                        % upper part: compare to ML MTP (toes)
                POSrel_R_MLtoe(iT)  = ( COP_rot(iT,1)-Toe_Rmed(1,iT)) ./ ...
                    (Toe_Rlat(1,iT)-Toe_Rmed(1,iT)) *100;
                POSrel_R_MLheel(iT) = NaN;
            end, clear MidFootAP  Ank_Rmed Ank_Rlat Toe_Rmed Toe_Rlat
            
%                     figure, hold on,
%                         plot(COP_filt(iT,1),COP_filt(iT,2),'k*'), plot(COP_filt2(iT,4),COP_filt(iT,5),'k*')
%                         plot(squeeze(Marker_filt(1,:,iT))',squeeze(Marker_filt(2,:,iT))','*b')
%                         plot(COP_rot(iT,1),COP_filt(iT,2),'r*'), plot(COP_rot(iT,4),COP_filt(iT,5),'r*')
%                         plot(squeeze(Marker_rot(1,:,iT))',squeeze(Marker_rot(2,:,iT))','*c')
            
        end, clear iT
    end
    
    %% Save Data Cons
    for k=1
        Data.(ConsNames{iCon}).FOR = FOR_filt;         Data.(ConsNames{iCon}).MOM = MOM_filt;
        Data.(ConsNames{iCon}).COP = COP_filt;         Data.(ConsNames{iCon}).MARK = Marker_filt;
        Data.(ConsNames{iCon}).MLabels = MarkerLabels;
        Data.(ConsNames{iCon}).COProt = COP_rot;         Data.(ConsNames{iCon}).MARKrot = Marker_rot;
        
        Data.(ConsNames{iCon}).POSrel_R_MLheel = POSrel_R_MLheel;
        Data.(ConsNames{iCon}).POSrel_R_MLtoe = POSrel_R_MLtoe;
        Data.(ConsNames{iCon}).POSrel_R_AP = POSrel_R_AP;
        Data.(ConsNames{iCon}).POSrel_L_MLheel = POSrel_L_MLheel;
        Data.(ConsNames{iCon}).POSrel_L_MLtoe = POSrel_L_MLtoe;
        Data.(ConsNames{iCon}).POSrel_L_AP = POSrel_L_AP;
    end
    
end

%%% Save Data

% save([DataPath,DataPathSub,'Data.mat'], 'Data')

 

%% Plot Results
        
figure, set(gcf,'color',[1 1 1])
        
for iCon = 1:size(ConsNames,1)
    subplot(ceil(size(ConsNames,1)/4),4,iCon), hold on, 
    ConNamePlot = ConsNames{iCon}; ConNamePlot( regexp(ConsNames{iCon},'_') )=' ';  title( ConNamePlot ), clear ConNamePlot
    % L
    if ~isempty( unique(round(Data.(ConsNames{iCon}).POSrel_L_AP( find(~isnan(Data.(ConsNames{iCon}).POSrel_L_AP))))) )
    plot(47, unique(round(Data.(ConsNames{iCon}).POSrel_L_AP( find(~isnan(Data.(ConsNames{iCon}).POSrel_L_AP))))), '.k')
    plot(unique(round(Data.(ConsNames{iCon}).POSrel_L_MLtoe( find(~isnan(Data.(ConsNames{iCon}).POSrel_L_MLtoe))))), 74, '.k')
    plot(unique(round(Data.(ConsNames{iCon}).POSrel_L_MLheel( find(~isnan(Data.(ConsNames{iCon}).POSrel_L_MLheel))) )), 34, '.k')
    end
    % R
    if ~isempty( unique(round(Data.(ConsNames{iCon}).POSrel_R_AP( find(~isnan(Data.(ConsNames{iCon}).POSrel_R_AP))))) ) && ...
           ~isempty( unique(round(Data.(ConsNames{iCon}).POSrel_R_MLheel( find(~isnan(Data.(ConsNames{iCon}).POSrel_R_MLheel)) ))) )
    plot(52, unique(round(Data.(ConsNames{iCon}).POSrel_R_AP( find(~isnan(Data.(ConsNames{iCon}).POSrel_R_AP))))), '.b')
    plot(unique(round(Data.(ConsNames{iCon}).POSrel_R_MLtoe( find(~isnan(Data.(ConsNames{iCon}).POSrel_R_MLtoe))))), 76, '.b')
    plot(unique(round(Data.(ConsNames{iCon}).POSrel_R_MLheel( find(~isnan(Data.(ConsNames{iCon}).POSrel_R_MLheel)) ))), 36, '.b')
    end
    
    ylim([0 120]), ylabel('Percentage Toe-Heel marker')
    xlim([0 110]), xlabel('Med      -    Percentage ML med-lat marker    -     Lat')
    
    Markers = [50 100; 0 75; 0 35; 50 0; 100 35; 100 75; 50 100];
    plot(Markers(:,1),Markers(:,2),'-c')

end 

% save fig
figname = [DataPath,DataPathSub,'Figs\ALLCONS'];
saveas(gcf,figname,'fig'), saveas(gcf,figname,'png')


