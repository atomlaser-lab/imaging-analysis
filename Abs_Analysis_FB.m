function [img,nd_image,fb] = Abs_Analysis_FB(varargin)
pause(1);
atomType = 'Rb87';
imaging_system = 'low res';
% imaging_system = 'low res';
tof = evalin('base', 'opt.tof'); % get the 'tof' variable from the base workspace
detuning = evalin('base', 'opt.detuning'); % get the 'tof' variable from the base workspace
% detuning = 8;
dispOD = [0,1.5]; % [0.0, 1.5];
plotOpt = 1;
plotROI = 0;
useFilt = 0;
filtWidth = 50e-6;
%% Set imaging region-of-interest (ROI)
roiRow = [1,2048]; % FULL FRAME 
roiCol = [1,2048]; % FULL FRAME 
% roiRow =  1600 + 150*[-1,1]; % 150
% roiCol =  1090 + 150*[-1,1]; % 150
% roiCol =  880 + 150*[-1,1]; % 150
roiStep = 1*[1,1];
% fittype = 'gauss2d';
fittype = 'gauss1d';
% fittype = 'twocomp2d';

% offset_region.row = [1600,1800];
% offset_region.col = [200,400];
offset_region.row = [];
offset_region.col = [];
%% Imaging parameters

imgconsts = AtomImageConstants(atomType,'tof',tof,'detuning',detuning,...
    'pixelsize',5.5e-6,'exposureTime',40e-6,'polarizationcorrection',1,'satOD',11);
redpower = evalin('base', 'opt.redpower');
raycus = evalin('base', 'opt.raycus');
imgconsts.freqs = 2*pi*get_trap_freq(raycus, redpower); % double power of raykus
% imgconsts.freqs = 2*pi*get_trap_freq(2*1,2);
if strcmpi(imaging_system,'high res')
    imgconsts.magnification = 3.49;% 3.5;
elseif strcmpi(imaging_system,'low res')
    imgconsts.magnification = 0.955;
elseif strcmpi(imaging_system,'vertical')
    imgconsts.magnification = 5.42;
    imgconsts.polarizationCorrection = 3/2;    %This is only approximate for sigma and pi polarised light
end
imgconsts.photonsPerCount = 0.4747;
image_rotation = -90;
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
try
    fb = FeedbackData.loadFeedbackData('directory',directory,'files',raw.getImageNumbers);
catch err
    fb = [];
end
numImages = numel(raw);
% plotOpt = plotOpt || numImages == 1;    %This always enables plotting if only one image is analyzed

img = AbsorptionImage.empty;
for nn = 1:numImages
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
    elseif size(img(jj).raw.images,3) >= 3
        img(jj).makeImage([1,2,3]);
    elseif size(img(jj).raw.images,3) > 3
        img(jj).makeImage(size(img(jj).raw.images,3) + (-2:0));
    else
        error('Not sure what to do here');
    end
    if useFilt
        img(jj).butterworth2D(filtWidth);
    end
    %
    % Fit clouds
    %
%     img(jj).clouds.fitdata.ex = 6;
    img(jj).fit;

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
            axs(end).Title.String = [axs(end).Title.String,', ',sprintf('N = %.2e,x = %.0f,y = %.0f',img(jj).clouds(1).N,img(jj).clouds(1).pos(1)/(imgconsts.pixelSize/imgconsts.magnification),img(jj).clouds(1).pos(2)/(imgconsts.pixelSize/imgconsts.magnification))];
            %hold(axs(end),'on');
            %plot(x corners, y corners,'--');
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
    if plotOpt
        [labelStr,numStr] = img(jj).labelOneROI;
        if jj == 1
            disp(labelStr);
        end
        %     for nn = 1:numel(img(jj).clouds)
        %         [labelStr,numStr] = img(jj).labelOneROI;
        disp(numStr);
    end

end


if (numImages == 1 && raw.is_multi_camera && size(raw.images{1},3) > 1) || (numImages == 1 && size(raw.images,3) >= 5)
% Old things for debugging 
%     sum_idx_y = 500:580;
%     sum_idx_x = 925:1025;
%     row = 400:550;
%     col = 950:1150;
%     raw.images = pagetranspose(flipud(raw.images{1}));
    
    %%% Save the image and get rows and columns
    raw.images = raw.images{1};
    row = 1:size(raw.images,1);
    col = 1:size(raw.images,2);

% Old things for debugging: 
%     max_idx = size(raw.images,3) - 3;
%     nd_image = raw.images(:,:,3:max_idx)./raw.images(:,:,2) - 1;

    %%% Obtain the number of images, and calculate signal and remove bad images
    max_idx = size(raw.images,3);                                           % Count the number of ND images
    nd_image = raw.images(:,:,2:max_idx)./raw.images(:,:,1) - 1;            % Get the signal
    nd_image(isinf(nd_image) | isnan(nd_image)) = 0;                        % Remove bad images
    
    %%% Flip the images (equiv to rotation)
    nd_image = -nd_image(row,col,:);                                        % Flip image 
    nd_image = nd_image - sum(sum(nd_image(1:8,1:8,:),1),2)/128; % /64      % Take away minimum
    
    %%% Reindex images (accounts for dark image)
    if raw.images(:,:,1) > 10
        for nn = 1:size(nd_image,3)
            if raw.images(1,1,1 + nn) < 10
                nd_image(:,:,nn) = nd_image(:,:,nn - 1);
            end
        end
    end
    nd_image_orig = nd_image;
    

    % ~~~ ~~~ ~~~ ~~~ ~~~ 
    %%% From here to...
    if(0)
        %%% Apply a 4th order Butterworth filter to smooth the data:
        width = 15;                                                             % Width of the butterwork filter
        order = 4;                                                              % Order of the filter
        % nd_image = const.butterworth2D(nd_image,width,order);                   % Apply a 4th order BW filter... acts a LP filter
    
        %%% Original code by RT:
%         nd_image = const.butterworth2D(nd_image,15);                         % Apply a 4th order BW filter... acts a LP filter
    
        %%% Apply a low-pass filter to remove high freq noise
        order = 1;                                                              % Sets the filter order for dct filter
        % nd_image = dct_filter(nd_image,15,order);                               % Apply a discrte cosine transform filter (LP filter)    
    
        %%% Calculate the inverse Laplacian:
        for nn = 1:size(nd_image,3)
            %%% Calculate the Fourier transform:
    %         Y = fftshift(fft2(nd_image(:,:,nn)));
            Y = fft2(nd_image(:,:,nn));                                         % Apply a 2D Fourier transform    
            Y(1,1) = 0;                                                         % Regularisation parameter sets the origin to zero
            Y = fftshift(Y);                                                    % Shifts the data due to shift that is built into the Fourier transform
    
            %%% Generate the k vector that 
            kx = 2*pi/(imgconsts.pixelSize/imgconsts.magnification)*linspace(-0.5,0.5,size(nd_image,2));
            ky = 2*pi/(imgconsts.pixelSize/imgconsts.magnification)*linspace(-0.5,0.5,size(nd_image,1));
            
            %%% Make a mesh grid 
            [KX,KY] = meshgrid(kx,ky);
            K2 = KX.^2 + KY.^2;
            
            %%% Define the order 1 filter:
            F = 1./K2;
            F(isinf(F) | isnan(F)) = 1/eps;
    
            %%% Generate the image in real space:
            nd_image(:,:,nn) = real(ifft2(ifftshift(F.*Y)));                    % Apply my filter, and then shift back and then inverse transform:
        end
    
        %%% Apply a non-linear filter:
%         nd_image = nd_image.^6;                                                 % Apply a 6th order filter to supress low freq noise
        % nd_image = nd_image.^2;
    
        %%% Only include the ROI:
        % sum_idx_y = 1:size(nd_image,1);
        % sum_idx_x = 1:size(nd_image,2);
        row = 1:size(nd_image,1);                                               % Calculates the number of rows 
        col = 1:size(nd_image,2);                                               % Calculates the number of columns
        % sum_idx_y = 30:40;
        % sum_idx_x = 30:90;
        nd_image = nd_image(row,col);
    end 
    %%% ... here was disabled

    if plotOpt
        figure(3);clf;
        axes('position',[0.3,0.3,0.6,0.65]);
        imagesc(nd_image(:,:,1),[-Inf,Inf]);                                % Plots the first bright SGI:
        axis equal;
        axis tight; 
        colorbar;
        title(sprintf('%.5f',mean(max(max(nd_image(:,:,1),[],1),[],2))));
        
        %%% Plot the Y distribution
        axes('position',[0.075,0.35,0.15,0.6]);
        plot(sum(nd_image(:,:,1),2),size(nd_image(:,:,1),1):-1:1,'.-');
        grid on;
        xlim([-Inf,Inf]);
    %     xlim([0,2.5]);
        
        %%% Plot the X distribution
        axes('position',[0.1,0.075,0.8,0.15]);
        plot(1:size(nd_image(:,:,1),2),sum(nd_image(:,:,1),1),'.-');
        grid on;
        ylim([-Inf,Inf]);
    %     ylim([0,2.5]);
    end

    %%% Return the unfiltered image:
    if(1)
        nd_image = nd_image_orig;
    end 
else
    nd_image = [];
end

end

