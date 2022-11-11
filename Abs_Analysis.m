function [img,nd_image] = Abs_Analysis(varargin)

atomType = 'Rb87';
imaging_system = 'high res';
tof = 20e-3;
detuning = 0;
dispOD = [0,1];
plotOpt = 1;
plotROI = 0;
useFilt = 0;
filtWidth = 5e-6;
useJointFit = 0;
%% Set imaging region-of-interest (ROI)
roiRow = 1650 + 200*[-1,1]; %580,800
roiCol = repmat(1070 + 200*[-1,1],size(roiRow,1),1); %850,1400
% roiRow = [1,128];
% roiCol = [1,128];
%roiRow = [1480 + 125*[-1,1];
%          1750 + 125*[-1,1]];
% roiRow = roiRow([2,3,4],:);
%roiCol = repmat(1100 + 125*[-1,1],size(roiRow,1),1);


roiRow = [400,1600];
roiCol = [400,1600];
roiStep = 3*[1,1];
% fittype = 'gauss2d'; 
fittype = 'twocomp1d';

% %for a full frame image
%  roiRow = [1,2048];
%  roiCol = [1,2048];

offset_region.row = [100,300];
offset_region.col = [100,300];
% offset_region.row = [];
% offset_region.col = [];
%% Imaging parameters

imgconsts = AtomImageConstants(atomType,'tof',tof,'detuning',detuning,...
    'pixelsize',5.5e-6,'exposureTime',100e-6,'polarizationcorrection',1,'satOD',11);
imgconsts.freqs = 2*pi*get_trap_freq(0.725,1.34);

if strcmpi(imaging_system,'low res')
    imgconsts.magnification = 0.6;
    imgconsts.photonsPerCount = 0.2644;
    image_rotation = 90;
elseif strcmpi(imaging_system,'high res')
    imgconsts.magnification = 3.5;
    imgconsts.photonsPerCount = 0.4747;
    image_rotation = -90;
elseif strcmpi(imaging_system,'vertical')
    imgconsts.magnification = 10.176;
    imgconsts.photonsPerCount = 0.4747;
    image_rotation = -90;
elseif strcmpi(imaging_system,'MOT')
    imgconsts.magnification = 0.75;
    imgconsts.photonsPerCount = 0.27;
    image_rotation = 90;
else
    warning('image system name wrong');
end

% directory = 'D:\raw-images';
directory = 'D:\labview-images';

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
            figure(1);clf;
            img(jj).plotAllData(dispOD,plotROI);
            axs = get(gcf,'children');
            %             axs(end).Title.String = [axs(end).Title.String,', ',sprintf('N1 = %.2e N2 = %.2e, R = %.2f',img(jj).clouds(1).N,img(jj).clouds(2).N,img(jj).clouds(2).N/img(jj).clouds(1).N)];
            axs(end).Title.String = [axs(end).Title.String,', ',sprintf('N = %.2e',img(jj).clouds(1).N)];
            
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


if numImages == 1 && raw.is_multi_camera
    nd_image = raw.images{1}(:,:,end) - raw.images{1}(:,:,1);
%     nd_image = const.butterworth2D(nd_image,3);
    sum_idx_y = 5:80;
    sum_idx_x = 100:170;

    figure(3);clf;
    axes('position',[0.3,0.3,0.6,0.65]);
    imagesc(nd_image,[-2^2,2^5]);axis equal;axis tight;
    colorbar;
    %
    % Plot the Y distribution
    %
    axes('position',[0.075,0.35,0.15,0.6]);
    plot(sum(nd_image(:,sum_idx_x),2),1:size(nd_image,1),'.-');
    grid on;
    xlim([0,600]);
    %
    % Plot the X distribution
    %
    axes('position',[0.1,0.075,0.8,0.15]);
    plot(1:size(nd_image,2),sum(nd_image(sum_idx_y,:),1),'.-');
    grid on;
    ylim([0,500]);
elseif numImages == 1 && size(raw.images,3) > 5
%     sum_idx_y = 500:580;
%     sum_idx_x = 925:1025;
%     row = 510:570;
%     col = 1040:1120;
    row = 1:size(raw.images,1);
    col = 1:size(raw.images,2);
%     max_idx = size(raw.images,3) - 3;
%     nd_image = raw.images(:,:,3:max_idx)./raw.images(:,:,2) - 1;

    max_idx = size(raw.images,3) - 4 + 1;
    nd_image = raw.images(:,:,2:max_idx)./raw.images(:,:,end) - 1;
    nd_image(isinf(nd_image) | isnan(nd_image)) = 0;
    nd_image = nd_image(row,col,:);

%     nd_image = const.butterworth2D(nd_image,10);

%     for nn = 1:size(nd_image,3)
%         Y = fftshift(fft2(nd_image(:,:,nn)));
%         kx = 2*pi/(imgconsts.pixelSize/imgconsts.magnification)*linspace(-0.5,0.5,size(nd_image,2));
%         ky = 2*pi/(imgconsts.pixelSize/imgconsts.magnification)*linspace(-0.5,0.5,size(nd_image,1));
%         [KX,KY] = meshgrid(kx,ky);
%         F = 1./(KX.^2 + KY.^2);
%         F(isinf(F)) = eps;
%         nd_image(:,:,nn) = real(ifft2(ifftshift(F.*Y)));
%     end

%     nd_image = nd_image.^2;

%     sum_idx_y = 1:size(nd_image,1);
%     sum_idx_x = 1:size(nd_image,2);
%     row = 1:size(nd_image,1);
%     col = 1:size(nd_image,2);
%     sum_idx_y = 30:40;
%     sum_idx_x = 30:90;
%     nd_image = nd_image(row,col);

    figure(3);clf;
    axes('position',[0.3,0.3,0.6,0.65]);
%     imagesc(nd_image,[0,2e-2]);axis equal;axis tight;
%     imagesc(nd_image,[min(min(nd_image(row,col))),max(max(nd_image(row,col)))]);axis equal;axis tight;
    imagesc(nd_image(:,:,1),[-Inf,Inf]);axis equal;axis tight;
%     ylim([min(row),max(row)]);
%     xlim([min(col),max(col)]);
    colorbar;
    %
    % Plot the Y distribution
    %
    axes('position',[0.075,0.35,0.15,0.6]);
    plot(sum(nd_image(:,:,1),2),1:size(nd_image(:,:,1),1),'.-');
    grid on;
    xlim([-Inf,Inf]);
    %
    % Plot the X distribution
    %
    axes('position',[0.1,0.075,0.8,0.15]);
    plot(1:size(nd_image(:,:,1),2),sum(nd_image(:,:,1),1),'.-');
    grid on;
    ylim([-Inf,Inf]);

    axes('position',[0.8,0.075,0.1,0.1]);
    imagesc(raw.images,1,[-Inf,Inf]);axis equal;axis tight;

else
    nd_image = [];
end

end

