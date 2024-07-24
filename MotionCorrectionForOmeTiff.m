%% Motion Correction For OME-Tiff file
%SimplifyingExtractTraceCode
clear,clc;
addpath(genpath('../NoRMCorre'));               % add the NoRMCorre motion correction package to MATLAB path
%gcp;                                            % start a parallel engine
CurrentPath=pwd;
% AllPath=subdir(fullfile(CurrentPath,'Flag.mat')); % list across day data
AllPath=subdir(fullfile(CurrentPath,'Image_scan_1_region_0_0.tif')); % list data of every day

for iAllPath=1:length(AllPath)
    
    [pathstr,~,~]=fileparts(AllPath(iAllPath).folder);
    cd(pathstr);
    %     cd(AllPath(iAllPath).folder);
    foldername = pwd;% folder where all the files are located.
    Temp=strsplit(foldername,'\');DataID=[Temp{1,end}];DataID=strrep(DataID,'_','-');
    MetaFiles=subdir(fullfile(foldername,'Experiment.xml'));% list raw data meta info location
    files = subdir(fullfile(foldername,'Image_scan_1_region_0_0.tif'));   % list of filenames (search all subdirectories)
    TemplateFile=subdir(fullfile(foldername,['*template' '*.tif'])); % list template files
    numFiles = length(files); % subfolder for 2P files
    MotionCorrFile=subdir(fullfile(foldername,'*MotionCond.mat')); %determine whether motion correction file existence
     
    for i=1:numFiles % go through each subfolder
        
        if isempty(MotionCorrFile)  % determine whether exist motion correction file, if not, continue
            
            [tree, ~, ~]=xml_read(MetaFiles(i).name); % read raw data meta info
            Planes=tree.ZStage(1).ATTRIBUTE.steps;% planes is imaging planes(multiple planes imaging)
            frameRate=tree.LSM.ATTRIBUTE.frameRate; pixelX=tree.LSM.ATTRIBUTE.pixelX; pixelY=tree.LSM.ATTRIBUTE.pixelY;
            ChannelNum=length(tree.Wavelengths.Wavelength); pixelSizeUM=tree.LSM.ATTRIBUTE.pixelSizeUM; % width/pixel
            frameNUM=tree.Streaming(1).ATTRIBUTE.frames;
            zFastEnable=tree.Streaming(1).ATTRIBUTE.zFastEnable; % determine zFast, 1(Yes), 0(No)
            flybackFrames=tree.Streaming(1).ATTRIBUTE.flybackFrames; % flybackFrames came from fastZ imaging(multiple planes imaging)
            if Planes ~= 1 && zFastEnable==1;  frameRate = frameRate/(Planes + flybackFrames); end
            
            %% ome-tiff file info
            Info=imfinfo(files(i).name);
            tif='tif';
            format=Info.Format;
            if (strcmp(format ,tif)==0); disp('file is not tif format, please load correct file'); end
            %% clear invalid tiff image
            % StripByteCounts in Info could show which image is invalid
            % if StripByteCounts=0 mean that it's invalid
            for temp1= 1:length(Info)
                if Info(temp1).StripByteCounts==0 && Info(temp1+5).StripByteCounts==0 && Info(temp1+10).StripByteCounts==0
                    break;
                end
            end
            %% load tiff file
            if Info(temp1).StripByteCounts==0 %&& ChannelNum==2 && tree.PMT.ATTRIBUTE.enableB==1
                Slice=temp1-1;
            elseif Info(temp1).StripByteCounts~=0 && ChannelNum==2 && tree.PMT.ATTRIBUTE.enableB==1 % 2 channel and imaged all frames
                Slice=temp1/2;
            elseif Info(temp1).StripByteCounts~=0 && ChannelNum==1% 1 channel and imaged all frames
                Slice=temp1;
            end
            Width=Info.Width;  Height=Info.Height;
            Image=zeros(Height,Width,Slice);
            tic;
            for itr=1:Slice
                try
                    Image(:,:,itr)=imread(files(i).name,itr);
                catch
                    Image(:,:,itr)=imread(files(i).name,itr+Planes); % if unexpectedly lost some frame, using next frame replace it
                end
            end
            toc;
            %figure; imagesc(Image(:,:,3)); colormap(gray);axis tight; axis equal;
            ChannelA=Image; clearvars Image;
            
            if zFastEnable==1 % multiple planes imaging
                
                for itrPlane = 1:Planes % load each planes data, ChannelA
                    
                    motion_correct = true;                            % perform motion correction
                    non_rigid = true;                                 % flag for non-rigid motion correction
                    output_type = 'h5';                               % format to save registered files
                    
                    if non_rigid; append0 = ['-Plane' num2str(itrPlane)]; append = '_nr'; else; append0 = ['-Plane' num2str(itrPlane) ]; append = '_rig'; end  % use this to save motion corrected files
                    % change diff planes append name
                    options_mc = NoRMCorreSetParms('d1',pixelX,'d2',pixelY,'grid_size',[128,128],'init_batch',200,...
                        'overlap_pre',32,'mot_uf',4,'bin_width',200,'max_shift',24,'max_dev',8,'us_fac',50,...
                        'output_type',output_type,'fr',frameRate);
                    
                    Planei=single(ChannelA(:,:,itrPlane:Planes:end)); %read single plane images
                    %% determine template
                    if i==1
                        template=[];
                    else
                        MotionCondFile=subdir(fullfile(foldername,'*MotionCond.mat'));
                        tempT=load(MotionCondFile(1).name);
                        template=tempT.MotionCond(itrPlane).template;
                    end
                    %template=single(imread(TemplateFile(itrPlane).name));
                    %figure;imagesc(template);colormap(gray);axis tight; axis equal;
                    %% start motion correct
                    col_shift = []; shifts=[];
                    fullname = files(i).name; [folder_name,file_name,ext] = fileparts(fullname);
                    output_filename = fullfile(folder_name,[file_name,append0,append,'.',output_type]);
                    options_mc = NoRMCorreSetParms(options_mc,'output_filename',output_filename,'h5_filename','','tiff_filename',''); % update output file name
                    tic;
                    if motion_correct
                        [M,shifts,template,options_mc,col_shift] = normcorre_batch_even(Planei,options_mc,template);
                            [~, ~, frmsAvgdORGNL, ~] = normcorre_batch_Rigid_4GPU(tempH5Name, NoRMCorreSetParms(options4motCrct), frmsAvgdORGNL, MaxCPUcoreNum4SPMD1, 1, 0, iterNum4step1, 0, 0);
                    else    % if files are already motion corrected convert them to h5
                        convert_file(fullname,'h5',fullfile(folder_name,[file_name,'_mc.h5']));
                    end
                    toc;
                    MotionCond(itrPlane).shifts=shifts; MotionCond(itrPlane).template=template; MotionCond(itrPlane).options_mc=options_mc; MotionCond(itrPlane).col_shift=col_shift;
                    
                    disp(['Performing MotionCorrection In Plane ',num2str(itrPlane),' out of ',num2str(Planes),' finished processing.'])
                    
                end
                               
                clearvars ChannelA M Planei; %release space
                save(fullfile(folder_name,[file_name,'_MotionCond','.mat']),'MotionCond','-v7.3');           % save shifts of each file at the respective folder

                %%  ChannelB
                if ChannelNum==2 && tree.PMT.ATTRIBUTE.enableB==1 % two channel imaging
                    
                    realFrameNum=frameNUM-(frameNUM/(Planes + flybackFrames))*flybackFrames; % real saved frame number, 2*realFrameNum included two channel
                    if 2*realFrameNum > length(Info) % max length(Info)=65536, so exist two Image_scan_1_region*.tif
                       
                        if (length(Info)-realFrameNum)>=Slice %if remained ChannelB >= valid ChannelA frame, then using Valid A Frame number
                            tempimage1=zeros(Height,Width,Slice);
                            for itr=1:Slice
                                tempimage1(:,:,itr)=imread(files(i).name,itr+realFrameNum);
                            end
                            tempimage2=[];
                        else
                            files1 = subdir(fullfile(foldername,'Image_scan_1_region_0_1.tif'));
                            %Info1=imfinfo(files1(i).name);
                            tempimage1=zeros(Height,Width,(length(Info)-realFrameNum));
                            for itr=1:(length(Info)-realFrameNum)
                                tempimage1(:,:,itr)=imread(files(i).name,itr+realFrameNum);
                            end
                            tempimage2=zeros(Height,Width,(Slice-(length(Info)-realFrameNum)));
                            for itr=1:(Slice-(length(Info)-realFrameNum))
                                tempimage2(:,:,itr)=imread(files1(i).name,itr);
                            end
                        end
                                                                     
                        ChannelB=cat(3,tempimage1,tempimage2);
                        clearvars tempimage1 tempimage2;                       
                    else %just one Image_scan_1_region*.tif file
                        Image=zeros(Height,Width,Slice);
                        tic;
                        for itr=1:Slice
                            try
                                Image(:,:,itr)=imread(files(i).name,itr+length(Info)/2);
                            catch
                                Image(:,:,itr)=imread(files(i).name,itr+length(Info)/2+Planes);
                            end
                        end
                        toc;
                        ChannelB=Image; clearvars Image;
                    end
                    %% perform motion correction in channelB
                    for itrPlane = 1:Planes % load each planes data
                        motion_correct = true;                            % perform motion correction
                        non_rigid = true;                                 % flag for non-rigid motion correction
                        output_type = 'h5';                               % format to save registered files
                        
                        if non_rigid; append0 = ['-Plane' num2str(itrPlane)]; append = '_Bnr'; else; append0 = ['-Plane' num2str(itrPlane)]; append = '_Brig'; end  % use this to save motion corrected files
                        % change diff planes append name
                        options_mc_B = NoRMCorreSetParms('d1',pixelX,'d2',pixelY,'grid_size',[128,128],'init_batch',200,...
                            'overlap_pre',32,'mot_uf',4,'bin_width',200,'max_shift',24,'max_dev',8,'us_fac',50,...
                            'output_type',output_type,'fr',frameRate);
                        
                        PlaneBi=single(ChannelB(:,:,itrPlane:Planes:end)); %read single plane images
                        %% determine template
%                         if i==1
                            template_B=[];
%                         else
%                             MotionCondFile=subdir(fullfile(foldername,'*MotionCond.mat'));
%                             load(MotionCondFile(1).name);
%                             template_B=MotionCond(itrPlane).template_B;
%                         end
                        %% start motion correct
                        col_shift_B = [];
                        fullname = files(i).name; [folder_name,file_name,ext] = fileparts(fullname);
                        output_filename = fullfile(folder_name,[file_name,append0,append,'.',output_type]);
                        options_mc_B = NoRMCorreSetParms(options_mc_B,'output_filename',output_filename,'h5_filename','','tiff_filename',''); % update output file name
                        tic;
                        if motion_correct
                            [M,shifts_B,template_B,options_mc_B,col_shift_B] = normcorre_batch_even(PlaneBi,options_mc_B,template_B);
                        else    % if files are already motion corrected convert them to h5
                            convert_file(fullname,'h5',fullfile(folder_name,[file_name,'_mc.h5']));
                        end
                        toc;
                        MotionCond(itrPlane).shifts_B=shifts_B; MotionCond(itrPlane).template_B=template_B; MotionCond(itrPlane).options_mc_B=options_mc_B; MotionCond(itrPlane).col_shift_B=col_shift_B;
                        
                        disp(['Performing MotionCorrection ChannelB In Plane ',num2str(itrPlane),' out of ',num2str(Planes),' finished processing.'])
                                               
                    end
                    
                    clearvars ChannelB M PlaneBi; %release space
                    save(fullfile(folder_name,[file_name,'_MotionCond','.mat']),'MotionCond','-v7.3');
                end
            else % one plane imaging
                motion_correct = true;                            % perform motion correction
                non_rigid = true;                                 % flag for non-rigid motion correction
                output_type = 'h5';                               % format to save registered files
                if non_rigid; append = '_nr'; else; append = '_rig'; end  % use this to save motion corrected files
                options_mc = NoRMCorreSetParms('d1',pixelX,'d2',pixelY,'grid_size',[128,128],'init_batch',200,...
                    'overlap_pre',32,'mot_uf',4,'bin_width',200,'max_shift',24,'max_dev',8,'us_fac',50,...
                    'output_type',output_type,'fr',frameRate);
                
                Planei=single(ChannelA);
                % plane1=mean(Planei(:,:,2000:2200),3);
                %figure; imagesc(plane1);colormap(gray); axis tight; axis equal;
                %corrcoef(template,plane1)
                %% determine template
                template=[];
%                 if i==1 && ~isempty(TemplateFile)
%                     template=single(imread(TemplateFile(1).name));
%                 elseif i==1 && isempty(TemplateFile)
%                     template=[];
%                 else
%                     template=MotionCond.template;
%                 end
                %figure;imagesc(template);colormap(gray);axis tight; axis equal;
                %% start motion correct
                col_shift = [];
                fullname = files(i).name; [folder_name,file_name,ext] = fileparts(fullname);
                output_filename = fullfile(folder_name,[file_name,append,'.',output_type]);
                options_mc = NoRMCorreSetParms(options_mc,'output_filename',output_filename,'h5_filename','','tiff_filename',''); % update output file name
                tic;
                if motion_correct
                    [M,shifts,template,options_mc,col_shift] = normcorre_batch_even(Planei,options_mc,template);
                else    % if files are already motion corrected convert them to h5
                    convert_file(fullname,'h5',fullfile(folder_name,[file_name,'_mc.h5']));
                end
                toc;
                MotionCond.shifts=shifts; MotionCond.template=template; MotionCond.options_mc=options_mc; MotionCond.col_shift=col_shift;
                
                %% channelB
                Image=zeros(Height,Width,Slice);
                tic;
                for itr=1:Slice
                    try
                        Image(:,:,itr)=imread(files(i).name,itr+length(Info)/2);
                    catch
                        Image(:,:,itr)=imread(files(i).name,itr+length(Info)/2+Planes);
                    end
                end
                toc;
                ChannelB=Image; clearvars Image;
                
                if ChannelNum == 1
                    clearvars ChannelA M Planei; %release space
                    save(fullfile(folder_name,[file_name,'_MotionCond','.mat']),'MotionCond','-v7.3');           % save shifts of each file at the respective folder
                else
                    clearvars ChannelA M Planei; %release space
                    save(fullfile(folder_name,[file_name,'_MotionCond','.mat']),'MotionCond','ChannelB','-v7.3');
                end
            end
            disp(['Performing MotionCorrection In SubFile ',num2str(i),' out of ',num2str(numFiles),' finished processing.'])
            
        end
    end
    disp([DataID,' Motion Correction Finished'])
end
cd(CurrentPath);


%% show some motion correction movie
ShowCorrResults=0;
if ShowCorrResults==1
    
    foldername = pwd;% folder where all the files are located.
    MFiles=subdir(fullfile(foldername,['Image_scan','*.h5']));
    
    Y=read_file(MFiles(2).name);
    Planei=read_file(MFiles(4).name); 
    M2=Planei;
    %Y=read_file(M);
    
    % Example plot planes1 and planes2
    % T=length(ChAMotionCo{1,1})/10;
    T=length(Y);
    nnY = quantile(Y(:),0.005); mmY = quantile(Y(:),0.995);
    
    figure;
    for t = 1:1:500 % first 1 is initiate, second 1 is increase value, T is terminal value
        subplot(121);imagesc(Y(:,:,t),[nnY,mmY]); xlabel('corrected data','fontsize',14,'fontweight','bold'); axis equal; axis tight;
        title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); %colormap('bone')
        %set(gca,'XTick',[],'YTick',[]);
        subplot(122);imagesc(M2(:,:,t),[nnY,mmY]); xlabel('raw data','fontsize',14,'fontweight','bold'); axis equal; axis tight;
        title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); %colormap('bone')
        %set(gca,'XTick',[],'YTick',[]);
        drawnow;
        pause(0.01);
    end
    %% preview each plane's image
    figure;
    for it=1:100
        imagesc(Planei(:,:,it));colormap('bone');
        axis equal; axis tight;
        drawnow;
        pause(0.05);
        
    end
end