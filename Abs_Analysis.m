function [img,nd_image,fb] = Abs_Analysis(varargin)

atomType = 'Rb87';
imaging_system = 'high res';
% imaging_system = 'low res';
tof = 20e-3;
detuning = 12;
dispOD = [0,1];
plotOpt = 0;
plotROI = 0;
useFilt = 1;
filtWidth = 5e-6;
useJointFit = 0;
%% Set imaging region-of-interest (ROI)
% pixel_marker = [1053,210]; %mag zero for horizontal
pixel_marker = [1053,265]; %Target optical trap location
% pixel_marker = [900,611];
% pixel_marker = [];

roiRow = 1300 + 150*[-1,1];
% roiRow = [1370;1700] + 150*[-1,1];
roiCol = 1053 + 150*[-1,1];
% roiCol = [10,2000];
% roiCol = repmat(1053 + 200*[-1,1],size(roiRow,1),1);
roiStep = 2*[1,1];
% fittype = 'gauss2d'; %for thermal cloud
fittype = '2comp2d'; %for BEC
% fittype = {'2comp2d','gauss2d'};

% % % for a full frame image
% roiRow = [10,1900];
% roiCol = [50,2000];
% roiStep = 10*[1,1];
% offset_region.row = [100,300];
% offset_region.col = [100,300];
offset_region.row = [];
offset_region.col = [];
%% Imaging parameters

imgconsts = AtomImageConstants(atomType,'tof',tof,'detuning',detuning,...
    'pixelsize',5.5e-6,'exposureTime',4*40e-6,'polarizationcorrection',1.58,'satOD',11);
imgconsts.freqs = 2*pi*get_trap_freq(0.75,1.2,100e-6,150e-6);

if strcmpi(imaging_system,'low res')
    imgconsts.magnification = 0.809;
    imgconsts.photonsPerCount = 0.1126;
    image_rotation = -90;
elseif strcmpi(imaging_system,'high res')
    imgconsts.magnification = 3.049;
    imgconsts.photonsPerCount = 0.64;
%     imgconsts.detuning_function = @double_resonance;
    image_rotation = 0;
elseif strcmpi(imaging_system,'vertical')
    imgconsts.magnification = 5.637; % bottom of cell
    %     imgconsts.magnification = 2.31; % middle of cell
    imgconsts.photonsPerCount = 0.4747;
    image_rotation = -90;
else
    warning('image system name wrong');
end

% directory = 'D:\raw-images';
directory = 'D:\labview-images';
% directory = 'C:\Users\admin\Downloads';

%% Load raw data
if nargin == 0 || (nargin == 1 && strcmpi(varargin{1},'last')) || (nargin == 2 && strcmpi(varargin{1},'last') && isnumeric(varargin{2}))
    %
    % If no input arguments are given, or the only argument is 'last', or
    % if the arguments are 'last' and a numeric array, then load the last
    % image(s).  In the case of 2 arguments, the second argument specifies
    % the counting backwards from the last image
    %
    if nargin < 2
        idx = 1;
    else
        idx = varargin{2};
    end
    args = {'files','last','index',idx};
else
    %
    % Otherwise, parse arguments as name/value pairs for input into
    % RawImageData
    %
    if mod(nargin,2) ~= 0
        error('Arguments must occur as name/value pairs!');
    end
    args = varargin;
end
%
% This loads the raw image sets
%
raw = BinaryImageData.loadImageSets('directory',directory,'rotation',image_rotation,args{:});
% figure(2);clf;
% imagesc(const.butterworth2D(raw.images{1}(:,:,1) - raw.images{1}(:,:,end),5),[-Inf,Inf]);
% axis equal;axis tight;
% ylim(430 + 100*[-1,1]);xlim(980 + 100*[-1,1]);

numImages = numel(raw);
plotOpt = plotOpt || numImages == 1;    %This always enables plotting if only one image is analyzed

img = AbsorptionImage.empty;
for nn = 1:numImages
    %     if raw(nn).is_multi_camera
    %         raw(nn).images = raw(nn).images{2};
    %     end
    img(nn,1) = AbsorptionImage(BinaryImageData);
end


for jj = 1:numImages
    %
    % Copy immutable properties
    %
    img(jj).constants.copy(imgconsts);
    img(jj).raw.copy(raw(jj));
    if img(jj).raw.is_multi_camera
        img(jj).raw.images = img(jj).raw.images{2};
    end
    img(jj).setClouds(size(roiRow,1));
    img(jj).offset_region = offset_region;
    for nn = 1:numel(img(jj).clouds)
        img(jj).clouds(nn).fitdata.set('roirow',roiRow(nn,:),'roiCol',roiCol(nn,:),...
            'roiStep',roiStep,'fittype',fittype,'method','x');
    end
    %
    % Create image
    %
    if size(img(jj).raw.images,3) == 2
        img(jj).makeImage;
    elseif size(img(jj).raw.images,3) == 3
        img(jj).makeImage([1,2,3]);
    elseif size(img(jj).raw.images,3) > 3
        img(jj).makeImage(size(img(jj).raw.images,3) - 1 - 3 + (1:3));
        %         img(jj).makeImage(size(img(jj).raw.images,3) + (-2:0));
    else
        error('Not sure what to do here');
    end
    if useFilt
        img(jj).butterworth2D(filtWidth);
    end
    %
    % Fit clouds
    %
    if useJointFit
        img(jj).jointFit([1,2]);
    else
        img(jj).fit;
    end

    %% Plotting
    if plotOpt
        if numImages == 1
            %
            % Plot absorption data and marginal distributions when there is only 1 image
            %
            figure(10);clf;
            img(jj).plotAllData(dispOD,plotROI);
            axs = get(gcf,'children');
            temp = [axs(end).Title.String]; %%%% uncomment to save image %%%%
            %             axs(end).Title.String = [axs(end).Title.String,', ',sprintf('N1 = %.2e N2 = %.2e, R = %.2f',img(jj).clouds(1).N,img(jj).clouds(2).N,img(jj).clouds(2).N/img(jj).clouds(1).N)];
            axs(end).Title.String = [axs(end).Title.String,', ',sprintf('N = %.2e',img(jj).clouds(1).N)];
%             colormap("gray");

%             %%%% uncomment to save image %%%%
%             saveas(gcf,['C:\Users\admin\Desktop\matlab-control\gravimeter-interface\saved images\',erase(temp,'Image: '),'.png']) % save image
%             %%%%

            if (pixel_marker)
                %show the marker
                plot(axs(end),axs(end).XLim,[pixel_marker(2),pixel_marker(2)],'g--')
                plot(axs(end),[pixel_marker(1),pixel_marker(1)],axs(end).YLim,'g--')
            end

            %if theare are three raw images show them
            %             if size(raw.images,3) == 3
            %                 %plot them
            %                 axes('position',[0.8,0.8,0.1,0.1]);
            %                 imagesc(raw.images(:,:,1),[-Inf,Inf]);axis equal;axis tight;
            %                 axes('position',[0.8,0.6,0.1,0.1]);
            %                 imagesc(raw.images(:,:,2),[-Inf,Inf]);axis equal;axis tight;
            %                 axes('position',[0.8,0.4,0.1,0.1]);
            %                 imagesc(raw.images(:,:,3),[-Inf,Inf]);axis equal;axis tight;
            %             end
        else
            %
            % Plot only the absorption data in a grid when there is more than one image
            %
            %             if jj == 1
            %                 figure(3);clf;
            %                 dimSubPlot=ceil(sqrt(numImages));
            %             end
            %
            %             figure(3);
            %             subplot(dimSubPlot,dimSubPlot,jj);
            %             img(jj).plotAbsData(dispOD,plotROI);
            figure(1);clf;
            img(jj).plotAllData(dispOD,plotROI);
            pause(0.01);
        end
    end

    %% Print summaries
    [labelStr,numStr] = img(jj).labelOneROI;
    if jj == 1
        disp(labelStr);
    end
    %     for nn = 1:numel(img(jj).clouds)
    %         [labelStr,numStr] = img(jj).labelOneROI;
    disp(numStr);

end


if numImages == 1 && raw.is_multi_camera && size(raw.images{1},3) > 1
    nd_image = raw.images{1}(:,:,2:end) - raw.images{1}(:,:,1);
    nd_image = pagetranspose(nd_image);
    row = 1:size(nd_image,1);
    col = 1:size(nd_image,2);
    nd_image = nd_image - mean(mean(nd_image(1:8,1:8,:),1),2);
    nd_image = -nd_image;

    figure(3);clf;
    axes('position',[0.3,0.3,0.6,0.65]);
    imagesc(nd_image(:,:,1),[-5,5]);axis equal;axis tight;
    colorbar;
    %
    % Plot the Y distribution
    %
    axes('position',[0.075,0.35,0.15,0.6]);
    plot(sum(nd_image(:,col,1),2),1:size(nd_image,1),'.-');
    grid on;
%     xlim([0,600]);
    %
    % Plot the X distribution
    %
    axes('position',[0.1,0.075,0.8,0.15]);
    plot(1:size(nd_image,2),sum(nd_image(row,:,1),1),'.-');
    grid on;
%     ylim([0,500]);
else
    nd_image = [];
end

try
    fb = FeedbackData.loadFeedbackData('directory',directory,'files',raw.getImageNumbers);
catch err
    fb = [];
end

if ~isempty(fb)
    figure(4);clf;
    subplot(2,1,1);
    plot(fb.t,fb.xpos,'.-');
    plot_format('Sample','X position [px]','',8);
    subplot(2,1,2);
    plot(fb.t,fb.xwidth,'.-');
    plot_format('Sample','X width [px]','',8);
end

idx = find(strcmpi(args,'export'));
if ~isempty(idx) && args{idx + 1}
    assignin("base","img",img);
    assignin("base","nd",nd_image);
    assignin("base","fb",fb);
end

