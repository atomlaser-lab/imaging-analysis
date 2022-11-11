function img = Abs_Analysis_DualState(varargin)

atomType = 'Rb87';
imaging_system = 'high res';
tof = 20e-3;
detuning = 0;
dispOD = [0,1];
plotOpt = 1;
plotROI = 0;
useFilt = 0;
filtWidth = 5e-6;
%% Set imaging region-of-interest (ROI)

%fit types
%none
% sum
% gauss1d
% twocomp1d
% tf1d
% gauss2d
% tf2d
% twocomp2d

% roiRow{1} = [1,2048];
% roiCol{1} = 1200 + 600*[-1,1];
% roiRow{1} = [400,1750];
% roiCol{1} = [400,1750];
roiRow{1} = [1450,1850];
roiCol{1} = [900,1350];
roiStep{1} = 3*[1,1];
fittype{1} = 'twocomp2d';

roiRow{2} = [400,800]+50;
roiCol{2} = [900,1300];
% roiRow{2} = [1,2048];
% roiCol{2} = 1200 + 600*[-1,1];

roiStep{2} = 3*[1,1];
% fittype{2} = 'gauss1d';
fittype{2} = 'twocomp2d';

% offset_region.row = [100,300];
% offset_region.col = [100,300];
offset_region.row = [1,100];
offset_region.col = [1,100];
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
numImages = numel(raw);
plotOpt = plotOpt || numImages == 1;    %This always enables plotting if only one image is analyzed

img = AbsorptionImage.empty;
if size(raw(1).images,3) <= 3
    for nn = 1:numImages
        img(nn,1) = AbsorptionImage(BinaryImageData);
    end
else
    for nn = 1:numImages
        img(nn,1) = AbsorptionImage(BinaryImageData);
        img(nn,2) = AbsorptionImage(BinaryImageData);
    end
end

for kk = 1:2
    for jj = 1:numImages
        %
        % Copy immutable properties
        %
        img(jj,kk).constants.copy(imgconsts);
        img(jj,kk).raw.copy(raw(jj));
        if kk == 1
            img(jj,kk).raw.images = img(jj,kk).raw.images(:,:,[2,3,4]);
        else
            img(jj,kk).raw.images = img(jj,kk).raw.images(:,:,[1,3,4]);
        end
    
        
        img(jj,kk).setClouds(size(roiRow{kk},1));
        img(jj,kk).offset_region = offset_region;
        for nn = 1:numel(img(jj,kk).clouds)
            img(jj,kk).clouds(nn).fitdata.set('roirow',roiRow{kk}(nn,:),'roiCol',roiCol{kk}(nn,:),...
                'roiStep',roiStep{kk},'fittype',fittype{kk},'method','x');
        end
        %
        % Create image
        %
        img(jj,kk).makeImage([1,2,3]);
        if useFilt
            img(jj).butterworth2D(filtWidth);
        end
        %
        % Fit clouds
        %
        img(jj,kk).fit;
    
        %% Plotting
        if plotOpt
            if numImages == 1
                %
                % Plot absorption data and marginal distributions when there is only 1 image
                %
                figure(kk);clf;
                img(jj,kk).plotAllData(dispOD,plotROI);
                axs = get(gcf,'children');
                %             axs(end).Title.String = [axs(end).Title.String,', ',sprintf('N1 = %.2e N2 = %.2e, R = %.2f',img(jj).clouds(1).N,img(jj).clouds(2).N,img(jj).clouds(2).N/img(jj).clouds(1).N)];
                axs(end).Title.String = [axs(end).Title.String,', ',sprintf('F = %d, N = %.2e',kk,img(jj,kk).clouds(1).N)];
                
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
                figure(kk);clf;
                img(jj,kk).plotAllData(dispOD,plotROI);
                pause(0.01);
            end
        end
    
        %% Print summaries
        [labelStr,numStr] = img(jj,kk).labelOneROI;
        if jj == 1
            disp(labelStr);
        end
        disp(numStr);
    
    end
end


end

