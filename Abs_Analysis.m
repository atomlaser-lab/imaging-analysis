function img = Abs_Analysis(varargin)

atomType = 'Rb87';
tof = 20e-3;
detuning = 0;
dispOD = [0,3];
plotOpt = 1;
plotROI = 1; 
useFilt = 0;
filtWidth = 75e-6;
useJointFit = 0;
%% Set imaging region-of-interest (ROI)
roiRow = [100,400]; %580,800
roiCol = [1000,1300]; %850,1400
% roiRow = [1500,1700];
%           1100,2000];
% roiCol = repmat([1000,1200],size(roiRow,1),1);
% roiRow = [1000,2000];
% roiCol = repmat([400,1800],size(roiRow,1),1);
% plotROI = {[min(roiRow(:)),max(roiRow(:))],[min(roiCol(:)),max(roiCol(:))]};
roiStep = [1,1]*1;
fittype = '2comp1d';
% mag = 0.758

offset_region.row = [100,300];
offset_region.col = [100,300];
%% Imaging parameters
imgconsts = AtomImageConstants(atomType,'tof',tof,'detuning',detuning,...
            'pixelsize',5.5e-6,'freqs',2*pi*[127,24,129],...
            'exposureTime',100e-6,'polarizationcorrection',15/8,'satOD',11);

imaging_system = 'high res';
if strcmpi(imaging_system,'low res')
    imgconsts.magnification = 0.6058;
    imgconsts.photonsPerCount = 0.2644;
    image_rotation = 90;
elseif strcmpi(imaging_system,'high res')
    imgconsts.magnification = 3.3851;
    imgconsts.photonsPerCount = 0.4747;
    image_rotation = -90;
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

numImages = numel(raw);
plotOpt = plotOpt || numImages == 1;    %This always enables plotting if only one image is analyzed

img = AbsorptionImage.empty;
for nn = 1:numImages
    if raw(nn).is_multi_camera
        raw(nn).images = raw(nn).images{2};
    end
    img(nn,1) = AbsorptionImage(BinaryImageData);
end


for jj = 1:numImages
    %
    % Copy immutable properties
    %
    img(jj).constants.copy(imgconsts);
    img(jj).raw.copy(raw(jj));
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

end


