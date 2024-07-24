% extract trace from diff tasks
%MotionCorrectionForOmeTiff
clear,clc;
%gcp; 
addpath(genpath('utilities')); addpath(genpath('deconvolution'));
foldername = pwd; Temp=strsplit(foldername,'\');DataID=[Temp{1,end}];DataID=strrep(DataID,'_','-');%% folder where all the files are located.
MetaFiles=subdir(fullfile(foldername,'Experiment.xml'));% list raw data meta info location
subFolder=subdir(fullfile(foldername,['*','MotionCond.mat'])); % list subfolder of mat
% delete ChannelB folder
for iSub=1:length(subFolder) % delete ChannelB folder
    ChannelBFile=subdir(fullfile(subFolder(iSub).folder,'ChanB_Preview.tif'));
    if ~isempty(ChannelBFile) 
        subFolder(iSub)=[];
        MetaFiles(iSub)=[];
        break;
    end
end
SyncFiles=subdir(fullfile(foldername,'Episode001.h5'));% list sync events;
JavaFiles=subdir(fullfile(foldername,['*','ser']));% list java files(behavioral data);
% Task Events in Java file
SerialData=[];
for iitr=1:size(JavaFiles,1)
    tempFile=JavaFiles(iitr,1).name;
    tempSerialData=ser2mat(tempFile);
    SerialData=[SerialData; tempSerialData];
end 

TaskEvents=GetTaskEventsFromJavaFile(SerialData); % java file task events
TaskType=TaskEvents.TaskType;

TraceMatFile=subdir(fullfile(foldername,['*FilePart1','*.mat']));
numFiles = length(SyncFiles);
mkdir([foldername,'\Figure\ROITracesFolder\']);
mkdir([foldername,'\Figure\ROISpikeFolder\']);

PreROINum=[800 1000 1000 900]; % PP16 pre determine ROI number for each plane

for i = 1:numFiles  % following files except first part

    [tree, ~, ~]=xml_read(MetaFiles(i).name); % read raw data meta info
    frameRate=tree.LSM.ATTRIBUTE.frameRate; pixelX=tree.LSM.ATTRIBUTE.pixelX; pixelY=tree.LSM.ATTRIBUTE.pixelY;
    ChannelNum=length(tree.Wavelengths.Wavelength); pixelSizeUM=tree.LSM.ATTRIBUTE.pixelSizeUM; % width/pixel
    zFastEnable=tree.Streaming(1).ATTRIBUTE.zFastEnable; % determine zFast(multiple planes imaging), 1(Yes), 0(No)
    if zFastEnable==1; Planes=tree.ZStage(1).ATTRIBUTE.steps; elseif zFastEnable==0; Planes=1; end % planes is imaging planes(multiple planes imaging)
    flybackFrames=tree.Streaming(1).ATTRIBUTE.flybackFrames; % flybackFrames came from fastZ imaging(multiple planes imaging)
    if Planes ~= 1 && zFastEnable==1
        frameRate = frameRate/(Planes + flybackFrames);
    end
        
     %% load Motion Correction File    
     
    AllTraceFile=struct([]); %preallocating a struct
       
    for itr = 1:Planes   % go through each plane  
        %% get ROI template,First downsampling data
        tic;
        Files=subdir(fullfile(subFolder(i).folder,['*', '_nr.h5'])); % list h5 file(motion correction data) of each planes in subfolder
        [folder_name,file_name,~] = fileparts(Files(itr).name); % list folder path and file name
        index=1;
        if frameRate <=10; tsub = 3;
        elseif frameRate >10&&frameRate <=20; tsub = 5;
        else; tsub = 10;
        end  % degree of downsampling (for 30Hz imaging rate you can try also larger, e.g. 8-10)
        FOV=[pixelX pixelY];
        [data,Ts,F_dark]=DownsamplingData(Files,itr,tsub,FOV,index); % downsampling data
        
        %% Estimates the level of percentile filtering to be used in DF/F extraction
        try
            [~,cdf_val] = estimate_percentile_level(data.Y(:,:,1:round(length(data.Y)/5):end)); % choose five sample to calculate percentile of cdf
            TotalPertile=zeros; 
            itrr=round(length(cdf_val)/10);
            for itrp=1:9
                TotalPertile(itrp)=mean(cdf_val(((itrp-1)*itrr+1):itrp*itrr));
            end
            percentile=mean(TotalPertile);
        catch
            percentile=30;
        end
        %figure; imagesc(data.Y(:,:,405*2))
        %% downsampling is over, next set parameters
        
        sizY = data.sizY; d1=sizY(1); d2=sizY(2); T=sizY(3);
        K = PreROINum(itr);                             % number of components to be found
        tau = round(7/pixelSizeUM);                     % std of gaussian kernel (half size of neuron(default neural size:16um))
        p = 2;                                          % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
        
        options = CNMFSetParms(...
            'd1',d1,'d2',d2,...                         % dimensionality of the FOV
            'deconv_method','constrained_foopsi',...    % neural activity deconvolution method
            'p',p,...                                   % order of AR dynamics
            'gSig',[tau,tau],...                        % half size of neuron
            'merge_thr',0.80,...                        % merging threshold
            'max_size_thr',20,'min_size_thr',8,...      % max/min acceptable size for each component
            'spatial_method','regularized',...          % method for updating spatial components
            'df_prctile',percentile,...                 % take the 50% of background fluorescence to compute baseline fluorescence
            'fr',frameRate/tsub,...                     % frame rate
            'nb',1,...                                  % number of background components
            'gnb',1,...                                 % number of global background components
            'min_SNR',2.0,...                           % minimum SNR threshold, the component is rejected if SNR is below than this value,
            'space_thresh',0.35,...                     % space correlation threshold,the component is rejected if 'space_thresh' is below than this value,
            'cnn_thr',0.2,...                           % threshold for CNN classifier,the component is rejected if 'cnn_thr' is below than this value,
            'decay_time',1.2,...                        % indicator decay time, GCaMP6f:0.4s; GCaMP6s:1.28s, based on paper
            'make_avi',1 ...                            % make and save video
            );
        %% Data pre-processing
        [P,Y] = preprocess_data(data.Y,p);
        
        if i~=1
            tempFile=subdir(fullfile(foldername,['*','FilePart','*.mat']));
            tempFile=load(tempFile(1).name,'AllTraceFile');  
            center_or=tempFile(1).AllTraceFile(itr).center_or; 
            P.ROI_list=center_or;
%             center_or=[]; 
%             for itri=1:length(json_file)
%                 center_or{itri,1}=json_file(itri).centroid;
%             end
%             center_or=cell2mat(center_or);           
%             P.ROI_list=ROISite{itr,1};           
        end
       
        [Ain,Cin,bin,fin,center] = initialize_components(Y,K,tau,options,P);  % initialize
        
        %% display centers of found components
        Cn = correlation_image(Y,8);  Cn(isnan(Cn)==1)=0; Cn(isinf(Cn)==1)=0;% image statistic (only for display purposes)
        %Cn = correlation_image_max(Y,8);
        figure;
        imagesc(Cn,[min(Cn(:)),max(Cn(:))]); %colormap(gray);% Cn=mean(ChAMotionCo{itr,1},3);
        axis equal; axis tight; hold all;
        scatter(center(:,2),center(:,1),'mo');
        title('Center of ROIs found from initialization algorithm');
        drawnow;
        toc;
        namefn1 = ['ImagingPlanes' num2str(itr) '.fig'];namefn2 = ['ImagingPlanes' num2str(itr) '.png'];
        saveas(gcf,[folder_name,'\',namefn1]);
        saveas(gcf,[folder_name,'\',namefn2]);
        close all;
        
        %% update spatial components
        %Yr = reshape(Y,d,T);
        [A,b,Cin] = update_spatial_components(data.Yr,Cin,fin,[Ain,bin],P,options); %Cin ChA,fin,ChA,bin ChA,P ChA, options ChA; Ain ChB
        %% update temporal components
        P.p = 0;    % set AR temporarily to zero for speed
        [C,f,P,S,YrA] = update_temporal_components(data.Yr,A,b,Cin,fin,P,options);   % YrA: residual
        
        if isempty(TraceMatFile)&&i==1
            %% classify components ()
            rval_space = classify_comp_corr(Y,A,C,b,f,options);
            ind_corr = rval_space > options.space_thresh;           % components that pass the correlation test
            %% further classification with cnn_classifier
            try  % matlab 2017b or later is needed
                [ind_cnn,value] = cnn_classifier(A,[d1,d2],'cnn_model',options.cnn_thr);
            catch
                ind_cnn = true(size(A,2),1);                        % components that pass the CNN classifier
            end
            %% event exceptionality
            fitness = compute_event_exceptionality(C+YrA,options.N_samples_exc,options.robust_std);
            ind_exc = (fitness < options.min_fitness);
            
            %% select components
            keep = (ind_corr | ind_cnn) & ind_exc;
            %             keep = ind_corr & ind_cnn;
             % manually select ROI
             A_keep = A(:,keep);
             C_keep = C(keep,:);
             
             Ain=A_keep; Cin=C_keep; center = options.ssub*com(Ain,options.d1,options.d2); %recreated center
             %% manually refine components (optional)
             refine_components = true;  % flag for manual refinement
             if refine_components %&& i==1
                 [Ain,Cin,center] = manually_refine_components(Y,Ain,Cin,center,Cn,tau,options);
             end
            close all;
            %% re-update & reject ROI
             
            [A,b,Cin] = update_spatial_components(data.Yr,Cin,fin,[Ain,bin],P,options); %Cin ChA,fin,ChA,bin ChA,P ChA, options ChA; Ain ChB
            P.p = 0;    % set AR temporarily to zero for speed
            [C,f,P,S,YrA] = update_temporal_components(data.Yr,A,b,Cin,fin,P,options);   % YrA: residual
            rval_space = classify_comp_corr(Y,A,C,b,f,options);
            ind_corr = rval_space > options.space_thresh;    
            try  % matlab 2017b or later is needed
                [ind_cnn,value] = cnn_classifier(A,[d1,d2],'cnn_model',options.cnn_thr);
            catch
                ind_cnn = true(size(A,2),1);                        % components that pass the CNN classifier
            end
             keep = ind_corr & ind_cnn;
            
            %% view contour plots of selected and rejected components (optional)
            throw = ~keep;
            A_keep = A(:,keep);
            C_keep = C(keep,:);
            Coor_k = [];
            Coor_t = [];
            figure;
            ax1 = subplot(121); plot_contours(A(:,keep),Cn,options,0,[],Coor_k,[],1,find(keep)); title('Selected components','fontweight','bold','fontsize',14);
            ax2 = subplot(122); plot_contours(A(:,throw),Cn,options,0,[],Coor_t,[],1,find(throw));title('Rejected components','fontweight','bold','fontsize',14);
            linkaxes([ax1,ax2],'xy')
            namefn1 = ['Keep and Discard ROI_Plane' num2str(itr) '.fig'];namefn2 = ['Keep and Discard ROI_Plane' num2str(itr) '.png'];
            saveas(gcf,namefn1); %['planes',itr,'Keep and Discard ROI'],'fig')
            saveas(gcf,namefn2);  
            close all;
                       
            %% merge found components
            [Am,Cm,K_m,merged_ROIs,Pm,Sm] = merge_components(data.Yr,A_keep,b,C_keep,f,P,S,options);
            %%
            display_merging = 1; % flag for displaying merging example
            if and(display_merging, ~isempty(merged_ROIs))
                %for i = 1:randi(length(merged_ROIs))
                imer = randi(length(merged_ROIs));
                ln = length(merged_ROIs{imer});
                figure;
                set(gcf,'Position',[300,300,(ln+2)*300,300]);
                for j = 1:ln
                    subplot(1,ln+2,j); imagesc(reshape(A_keep(:,merged_ROIs{imer}(j)),d1,d2));
                    title(sprintf('Component %i',j),'fontsize',16,'fontweight','bold'); axis equal; axis tight;
                end
                subplot(1,ln+2,ln+1); imagesc(reshape(Am(:,K_m-length(merged_ROIs)+imer),d1,d2));
                title('Merged Component','fontsize',16,'fontweight','bold');axis equal; axis tight;
                subplot(1,ln+2,ln+2);
                plot(1:T,(diag(max(C_keep(merged_ROIs{imer},:),[],2))\C_keep(merged_ROIs{imer},:))');
                hold all; 
                plot(1:T,Cm(K_m-length(merged_ROIs)+imer,:)/max(Cm(K_m-length(merged_ROIs)+imer,:)),'--k')
                title('Temporal Components','fontsize',16,'fontweight','bold')
                drawnow;
                namefn1 = ['Merged ROI_Plane' num2str(itr) '.fig'];namefn2 = ['Merged ROI_Plane' num2str(itr) '.png'];
                saveas(gcf,namefn1); %['planes',itr,'Keep and Discard ROI'],'fig')
                saveas(gcf,namefn2); 
                close all;                
            end
            %% refine estimates excluding rejected components
            Pm.p = p;    % restore AR value
            [A2,b2,C2] = update_spatial_components(data.Yr,Cm,f,[Am,b],Pm,options);
            [C2,f2,P2,S2,YrA2] = update_temporal_components(data.Yr,A2,b2,C2,f,Pm,options);
                        
            [A_or,C_or,S_or,P_or] = order_ROIs(A2,C2,S2,P2); % order components
            K_m = size(C_or,1);
            %view_components(Y,A,C,b,f,Cn,options)
             
        end
        
        if isempty(TraceMatFile)&&i==1
            %% plot ROI map
            % Cn=AllTraceFile(itr).Cn; A_or=AllTraceFile(itr).A_or; options=AllTraceFile(itr).options;
            figure;
            [Coor,json_file] = plot_contours(A_or,Cn,options,1);
            save(fullfile(folder_name,[file_name,'_json_file']),'json_file','Coor');  % optional save json file with component coordinates (requires matlab json library)          
            namefn1 = ['ROI Order_Plane' num2str(itr) '.fig'];
            namefn2 = ['ROI Order_Plane' num2str(itr) '.png'];
            saveas(gcf,namefn1); saveas(gcf,namefn2);
            close all;
            %% ROI location site
            center_or = options.ssub*com(A_or,options.d1,options.d2); % ROI Order
            ROISite{itr,1}=center_or;
            save tempROISite.mat ROISite
            %[C_df1,Df1] = extract_DF_F(Y2,A_or,C_full,P_or,options);
            % display components
            %         plot_components_GUI(data.Yr,A_or,C_or,b2,f2,Cn,options);
            % show a video
            % 
            % view_components.m
            % extract residual signals for each trace
            if exist('YrA','var'); R_keep = YrA(keep,:); else; R_keep = compute_residuals(data,A_or,b,C_or,f); end
            %% extract fluorescence on native temporal resolution
            options.fr = options.fr*tsub;                   % revert to origingal frame rate
            N = size(C_or,1);                               % total number of components
            T = sum(Ts);                                    % total number of timesteps
            %T = sum(Ts)+round(frameRate*6);
            C_full = imresize(C_or,[N,T]);                  % upsample to original frame rate
            R_full = imresize(R_keep,[N,T]);                % upsample to original frame rate
            F_full = C_full + R_full;                       % full fluorescence
            f_full = imresize(f,[size(f,1),T]);             % upsample temporal background
            
            S_full = zeros(N,T);
            P_or.p = 0; % P_or.p = 0;
            ind_T = [0;cumsum(Ts(:))];
            %ind_T = [0;cumsum(T(:))];
            options.nb = options.gnb;
            inds = ind_T(1)+1:ind_T(1+1);   % indeces of file i to be updated
            ProcessData=read_file(Files(itr).name); %ProcessData=cat(3,ProcessData, ProcessData(:,:,(end-26:end)));
            [~,ProcessY] = preprocess_data(ProcessData,p);
            [C_full(:,inds),f_full(:,inds),~,S_full(:,inds),R_full(:,inds)] = update_temporal_components_fast(ProcessY,A_or,b,C_full(:,inds),f_full(:,inds),P_or,options);
            disp(['Extracting raw fluorescence at native frame rate. File ',num2str(itr),' out of ',num2str(Planes),' finished processing.'])
            % extract DF/F
            [C_df,Df] = extract_DF_F(ProcessY,A_or,C_full,P_or,options);
            [F_dff,F0] = detrend_df_f(A_or,[b,ones(prod(FOV),1)],C_full,[f_full;-double(F_dark)*ones(1,T)],R_full,options);
            
        else
            %% plot ROI map
            figure('Visible','Off');
            [Coor,json_file] = plot_contours(A,Cn,options,1);
            save(fullfile(folder_name,[file_name,'_json_file']),'json_file','Coor');  % optional save json file with component coordinates (requires matlab json library)
            %savejson('jmesh',json_file,'.mat');        % optional save json file with component coordinates (requires matlab json library)
            namefn1 = ['ROI Order_Plane' num2str(itr) '.fig'];namefn2 = ['ROI Order_Plane' num2str(itr) '.png'];
            saveas(gcf,[folder_name,'\',namefn1]);
            saveas(gcf,[folder_name,'\',namefn2]);
            close all;
            % extract residual signals for each trace
            if exist('YrA','var');  R_keep = YrA; else; R_keep = compute_residuals(data,A,b,C,f); end
            %% extract fluorescence on native temporal resolution
            options.fr = options.fr*tsub;                   % revert to origingal frame rate
            N = size(C,1);                                  % total number of components
            T = sum(Ts);                                    % total number of timesteps
            %T = sum(Ts)+round(frameRate*12);
            C_full = imresize(C,[N,T]);                     % upsample to original frame rate
            R_full = imresize(R_keep,[N,T]);                % upsample to original frame rate
            F_full = C_full + R_full;                       % full fluorescence
            f_full = imresize(f,[size(f,1),T]);             % upsample temporal background
            
            S_full = zeros(N,T);
            P.p = 0; % P_or.p = 0;
            ind_T = [0;cumsum(Ts(:))];
            %ind_T = [0;cumsum(T(:))];
            options.nb = options.gnb;
            inds = ind_T(1)+1:ind_T(1+1);   % indeces of file i to be updated
            ProcessData=read_file(Files(itr).name); %ProcessData=cat(3,ProcessData, ProcessData(:,:,((end-round(frameRate*12)+1):end)));
            [~,ProcessY] = preprocess_data(ProcessData,p);
            [C_full(:,inds),f_full(:,inds),~,S_full(:,inds),R_full(:,inds)] = update_temporal_components_fast(ProcessY,A,b,C_full(:,inds),f_full(:,inds),P,options);
            disp(['Extracting raw fluorescence at native frame rate. File ',num2str(itr),' out of ',num2str(Planes),' finished processing.'])
            % extract DF/F
            [C_df,Df] = extract_DF_F(ProcessY,A,C_full,P,options);
            [F_dff,F0] = detrend_df_f(A,[b,ones(prod(FOV),1)],C_full,[f_full;-double(F_dark)*ones(1,T)],R_full,options);
        end
        F_dff(isinf(F_dff)==1)=0; F_dff(isnan(F_dff)==1)=0;
        C_df(isinf(C_df)==1)=0; C_df(isnan(C_df)==1)=0;
     
        %% spike deconvolution      
        C_dec = zeros(N,T);         % deconvolved DF/F traces
        S_dec = zeros(N,T);         % deconvolved neural activity
        bl = zeros(N,1);            % baseline for each trace (should be close to zero since traces are DF/F)
        neuron_sn = zeros(N,1);     % noise level at each trace
        g = cell(N,1);              % discrete time constants for each trace
        if p == 1; model_ar = 'ar1'; elseif p == 2; model_ar = 'ar2'; else; error('This order of dynamics is not supported'); end
        
        for itr0 = 1:N
            spkmin = options.spk_SNR*GetSn(F_dff(itr0,:));
            lam = choose_lambda(exp(-1/(options.fr*options.decay_time)),GetSn(F_dff(itr0,:)),options.lam_pr);
            [cc,spk,opts_oasis] = deconvolveCa(F_dff(itr0,:),model_ar,'method','thresholded','optimize_pars',true,'maxIter',20,...
                'window',150,'lambda',lam,'smin',spkmin);
            bl(itr0) = opts_oasis.b;
            C_dec(itr0,:) = cc(:)' + bl(itr0);
            S_dec(itr0,:) = spk(:);
            neuron_sn(itr0) = opts_oasis.sn;
            g{itr0} = opts_oasis.pars(:)';
            disp(['Performing deconvolution. Trace ',num2str(itr0),' out of ',num2str(N),' finished processing.'])
        end
        C_dec(isinf(C_dec)==1)=0; C_dec(isnan(C_dec)==1)=0; %delete inf/nan
        S_dec(isinf(S_dec)==1)=0; S_dec(isnan(S_dec)==1)=0;
          %% plot trace
        %% plot trace
        if i==1; figure; else; figure('Visible','Off'); end%colormap('bone');figure;
        set(gcf,'position',[50 0 1920 1000]);
        for ip=1:size(C_df,1)%:10%:2
            plot(C_df(ip,:)+ip)
            hold on
        end
        namefn1 = ['ExampleTrace-FilePart' num2str(i) '-Plane' num2str(itr) '.fig'];
        namefn2 = ['ExampleTrace-FilePart' num2str(i) '-Plane' num2str(itr) '.png'];
        saveas(gcf,[foldername '\Figure\ROITracesFolder\', namefn1]);
        saveas(gcf,[foldername '\Figure\ROITracesFolder\', namefn2]);
        close all;
        % plot detrend_df_f
        if i==1; figure; else; figure('Visible','Off'); end%colormap('bone');figure;
        set(gcf,'position',[50 0 1920 1000]);
        for ip=1:size(F_dff,1)%:10%:2
            plot(F_dff(ip,:)+ip)
            hold on
        end
        namefn1 = ['ExampleF_dff-FilePart' num2str(i) '-Plane' num2str(itr) '.fig'];
        namefn2 = ['ExampleF_dff-FilePart' num2str(i) '-Plane' num2str(itr) '.png'];
        saveas(gcf,[foldername '\Figure\ROITracesFolder\', namefn1]);
        saveas(gcf,[foldername '\Figure\ROITracesFolder\', namefn2]);
        close all;
        %% plot spike
        figure%('Visible','Off'); %colormap('bone');figure;
        set(gcf,'position',[50 0 1920 1000]);
        for ip=1:size(S_dec,1)%:10%:2
            plot(S_dec(ip,:)+ip)
            hold on
        end
        namefn1 = ['ExampleSpike-FilePart' num2str(i) '-Plane' num2str(itr) '.fig'];
        namefn2 = ['ExampleSpike-FilePart' num2str(i) '-Plane' num2str(itr) '.png'];
        saveas(gcf,[foldername '\Figure\ROISpikeFolder\', namefn1]);
        saveas(gcf,[foldername '\Figure\ROISpikeFolder\', namefn2]);
        close all;
        %% plot C_dec
        figure%('Visible','Off'); %colormap('bone');figure;
        set(gcf,'position',[50 0 1920 1000]);
        for ip=1:size(C_dec,1)%:10%:2
            plot(C_dec(ip,:)+ip)
            hold on
        end
        namefn1 = ['ExampleC_dec-FilePart' num2str(i) '-Plane' num2str(itr) '.fig'];
        namefn2 = ['ExampleC_dec-FilePart' num2str(i) '-Plane' num2str(itr) '.png'];
        saveas(gcf,[foldername '\Figure\ROISpikeFolder\', namefn1]);
        saveas(gcf,[foldername '\Figure\ROISpikeFolder\', namefn2]);
        close all;
        
        
        %%
        if TaskType==1 % DPA Task
            [SampleClip,TestClip,AllOdorClip,TrialStart,TrialStartClip,TrialEndClip,trialLength,~,~,~] = GetEventsFromDiffTaskin2P(SyncFiles,i,itr,Planes,zFastEnable,frameRate,TaskType,flybackFrames); % take diff events from diff task
        else % DRT or dual task
            [SampleClip,TestClip,AllOdorClip,TrialStart,TrialStartClip,TrialEndClip,trialLength,AbolishClip,AbolishStart,...
                AbolishEnd] = GetEventsFromDiffTaskin2P(SyncFiles,i,itr,Planes,zFastEnable,frameRate,TaskType,flybackFrames); % take diff events from diff task
        end
%         tempF_dff=F_dff; tempC_df=C_df; tempC_dec=C_dec; tempS_dec=S_dec;
%         % tempF_dff(:,TrialEndClip(1,46)+1:end)=[];  tempC_df(:,TrialEndClip(1,46)+1:end)=[];
%         % tempC_dec(:,TrialEndClip(1,46)+1:end)=[]; tempS_dec(:,TrialEndClip(1,46)+1:end)=[];
%         tempF_dff=[tempF_dff tempF_dff(:,5548:5568)];
%         tempC_df=[tempC_df tempC_df(:,5548:5568)];
%         tempC_dec=[tempC_dec tempC_dec(:,5548:5568)];
%         tempS_dec=[tempS_dec tempS_dec(:,5548:5568)];
%         F_dff=tempF_dff; C_df=tempC_df; C_dec=tempC_dec; S_dec=tempS_dec;
        % assign trace to diff task type
        if TaskType==1 || TaskType==5 || TaskType==6% DPA or DRT-DPA dual task
            AbolishStart=[]; AbolishEnd=[];
            [AllCdf,deconSpike,AllDdf,RawTrace,dFF0,~,~,~,~,~] = AssignFluoTraceToDiffTask(F_dff,...
                C_df,C_dec,S_dec,TrialStartClip,trialLength,AbolishStart,AbolishEnd,frameRate,N,TaskType);
        else % DRT/Dual task
            [AllCdf,deconSpike,AllDdf,RawTrace,dFF0,AboDff,AboTrace,AboDec,AboSpike,AbodFF0] = AssignFluoTraceToDiffTask(F_dff,...
                C_df,C_dec,S_dec,TrialStartClip,trialLength,AbolishStart,AbolishEnd,frameRate,N,TaskType);
        end
                
        %% save in a struct
                       
        AllTraceFile(itr).dFF0=dFF0; AllTraceFile(itr).F_dff=AllCdf;  AllTraceFile(itr).C_df=RawTrace;  AllTraceFile(itr).C_dec=AllDdf;  AllTraceFile(itr).S_dec=deconSpike;
        AllTraceFile(itr).Cn=Cn; AllTraceFile(itr).options=options;
        %         if i==1 AllTraceFile(itr).A_or=A_or; AllTraceFile(itr).center_or=center_or; else AllTraceFile(itr).A_or=A; AllTraceFile(itr).center_or=center1; end
        if isempty(TraceMatFile)&&i==1; AllTraceFile(itr).A_or=A_or; else;  AllTraceFile(itr).A_or=A; end
        AllTraceFile(itr).center_or=center_or; AllTraceFile(itr).Df=Df;
        AllTraceFile(itr).AllOdorClip=AllOdorClip; AllTraceFile(itr).SampleClip=SampleClip; AllTraceFile(itr).TestClip=TestClip;
        AllTraceFile(itr).TrialStartClip=TrialStartClip; AllTraceFile(itr).TrialEndClip=TrialEndClip;
        % abolished trials
        if  exist('AbolishClip', 'var')
            if ~isempty(AbolishClip)
                AllTraceFile(itr).AbodFF0=AbodFF0; AllTraceFile(itr).AboDff=AboDff;  AllTraceFile(itr).AboTrace=AboTrace;  AllTraceFile(itr).AboDec=AboDec;  AllTraceFile(itr).AboSpike=AboSpike;
                AllTraceFile(itr).AbolishClip=AbolishClip; AllTraceFile(itr).AbolishStart=AbolishStart; AllTraceFile(itr).AbolishEnd=AbolishEnd;
            end
        end
        disp(['Performing Extract Trace In Plane ',num2str(itr),' out of ',num2str(Planes),' finished processing.'])
        
%         save([DataID,'FilePart',num2str(i)],'AllTraceFile')
        
    end
    
    %% save in mat
    %AA={}; BB={}; CC={}; DD={}; EE={};
    AA=cell(length(AllTraceFile),1); 
    BB=cell(length(AllTraceFile),1); 
    CC=cell(length(AllTraceFile),1); 
    DD=cell(length(AllTraceFile),1); 
    EE=cell(length(AllTraceFile),1);
 
    for n=1:length(AllTraceFile)
        temp1=AllTraceFile(n).F_dff;
        AA{n,1}=temp1;
        temp2=AllTraceFile(n).C_df;
        BB{n,1}=temp2;
        temp3=AllTraceFile(n).C_dec;
        CC{n,1}=temp3;
        temp4=AllTraceFile(n).S_dec;
        DD{n,1}=temp4;
        temp5=AllTraceFile(n).dFF0;
        EE{n,1}=temp5;
    end
    %SumTraceFile=[];
    SumTraceFile.AllTraceFile=AllTraceFile;
    SumTraceFile.TaskEvents=TaskEvents;
    SumTraceFile.SumF_dff= cat(1,AA{:,1});
    SumTraceFile.SumC_df= cat(1,BB{:,1});
    SumTraceFile.SumC_dec= cat(1,CC{:,1});
    SumTraceFile.SumS_dec= cat(1,DD{:,1});
    SumTraceFile.SumdFF0= cat(1,EE{:,1});
    
    save([DataID,'FilePart',num2str(i)],'SumTraceFile','AllTraceFile','TaskEvents','options','DataID','frameRate','trialLength','AllOdorClip','-v7.3') 
    
    disp(['Performing Extract Trace In SubFile ',num2str(i),' out of ',num2str(numFiles),' finished processing.'])
    
    clearvars -except DataID foldername MetaFiles subFolder SyncFiles JavaFiles TraceMatFile numFiles i PreROINum TaskEvents TaskType
    
end
%end

CombineFiles=true;
if CombineFiles
    clear;
    % combine all subfiles to one
    foldername = pwd; Temp=strsplit(foldername,'\');DataID=[Temp{1,end}];DataID=strrep(DataID,'_','-');%% folder where all the files are located.
    TraceMatFile=subdir(fullfile(foldername,['*NewTraceSubFile','*.mat']));
    if ~isempty(TraceMatFile)
        subFile=subdir(fullfile(foldername,'*NewTraceSubFile*.mat'));
    else
        subFile=subdir(fullfile(foldername,'*FilePart*.mat'));
    end
    Plane=load(subFile(1).name,'AllTraceFile'); Plane=length(Plane.AllTraceFile);
    ImageParameter=load(subFile(1).name,'frameRate');
    %AllFileF_dff{length(subFile),1}=[];  
    AllFileF_dff=cell(length(subFile),1); 
    AllFileC_df=cell(length(subFile),1);   
    AllFileC_dec=cell(length(subFile),1); 
    AllFileS_dec=cell(length(subFile),1);   
    AllFiledFF0=cell(length(subFile),1);  
    AllFileTrace=struct([]); %preallocating a struct 
    
    %AllFileF_dff=[]; AllFileC_df=[]; AllFileC_dec=[]; AllFileS_dec=[]; AllFiledFF0=[]; AllFileTrace=[];
    tic;
    for nn = 1:Plane % go through each planes
        for n=1:length(subFile)
            tempfile=load(subFile(n).name,'AllTraceFile');
            AllFileF_dff{n,1}=tempfile.AllTraceFile(nn).F_dff;
            AllFileC_df{n,1}=tempfile.AllTraceFile(nn).C_df;
            AllFileC_dec{n,1}=tempfile.AllTraceFile(nn).C_dec;
            AllFileS_dec{n,1}=tempfile.AllTraceFile(nn).S_dec;
            AllFiledFF0{n,1}=tempfile.AllTraceFile(nn).dFF0;
        end
        AllFileTrace(nn).F_dff=cat(3,AllFileF_dff{:,1});
        AllFileTrace(nn).C_df=cat(3,AllFileC_df{:,1});
        AllFileTrace(nn).C_dec=cat(3,AllFileC_dec{:,1});
        AllFileTrace(nn).S_dec=cat(3,AllFileS_dec{:,1});
        AllFileTrace(nn).dFF0=cat(3,AllFiledFF0{:,1});
        
    end
    toc;
%% 20231225
%     Trials=cell(length(subFile),1); Sample=cell(length(subFile),1); Test=cell(length(subFile),1); Laser=cell(length(subFile),1);
%     AllFileEvents=struct([]);
% 
%     for in=1:Plane
%         tempfile1=load(subFile(in).name,'TaskEvents');
%         Trials{in,1}=tempfile1.TaskEvents.Trials;  
%         Sample{in,1}=tempfile1.TaskEvents.Sample;  
%         Test{in,1}=tempfile1.TaskEvents.Test;  
%         Laser{in,1}=tempfile1.TaskEvents.Laser;          
%     end
% 
%    AllFileEvents(1).Trials=cat(1,Trials{:,1});
%    AllFileEvents(1).Sample=cat(1,Sample{:,1});
%    AllFileEvents(1).Test=cat(1,Test{:,1});
%    AllFileEvents(1).Laser=cat(1,Laser{:,1});

   %% 20231225
    
    %AllFileTrials=[]; AllFileSample=[]; AllFileTest=[]; AllFileLaser=[];
    %AllFileEvents=struct([]);
    tempfile1=load(subFile(1).name,'TaskEvents');
    AllFileEvents=tempfile1.TaskEvents; % AllFileEvents=TaskEvents;

    
    if ~isempty(TraceMatFile) % already exist trace mat files, so creat a new folder
        mkdir(foldername,[DataID,'-NewTrace','-SumFluoTraceFile']);
        NewTracePath=fullfile(foldername,[DataID,'-NewTrace','-SumFluoTraceFile']);
        cd(NewTracePath);
        save([DataID,'-SumFluoTraceFile'],'AllFileTrace','AllFileEvents','ImageParameter')
        cd(foldername);
    else
        save([DataID,'-SumFluoTraceFile'],'AllFileTrace','AllFileEvents','ImageParameter','-v7.3')
    end
    
end

