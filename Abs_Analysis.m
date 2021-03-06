function img = Abs_Analysis(varargin)

atomType = 'Rb87';
tof = 30e-3;

dispOD = [0,1];
plotOpt = 0;
plotROI = 0;
useFilt = 0;
filtWidth = 50e-6;
useJointFit = 0;
%% Imaging first spot
roiCol = [300,500];
roiRow = repmat([500,750],size(roiCol,1),1);
roiStep = 20;
fittype = 'tf2d'; %gauss2d

%% Imaging Second spot
% roiRow = [150,550];
% roiCol = repmat([550,850],1,1);
% roiRow = [100,470;
%           470,800];
% roiCol = repmat([550,850],2,1);
% roiStep = 4;
% fittype = 'tf2d';

%% Fit parameters
% lb = CloudParameters('becwidth',[650e-6,650e-6]);
% ub = CloudParameters('becwidth',[850e-6,850e-6]);
% guess = CloudParameters('becwidth',[750e-6,750e-6]);

%% Imaging parameters
imgconsts = AtomImageConstants(atomType,'tof',tof,'detuning',0,...
            'pixelsize',6.45e-6,'magnification',0.99,...
            'freqs',2*pi*[53,53,25],'exposureTime',30e-6,...
            'polarizationcorrection',1.5,'satOD',5);

% directory = '../raw-images';
directory = '\\TANIT\2021';
% directory = 'E:\RawImages\2021';

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
raw = RawImageData.loadImageSets('directory',directory,args{:});

numImages = numel(raw);
plotOpt = plotOpt || numImages == 1;    %This always enables plotting if only one image is analyzed

img = AbsorptionImage.empty;
for nn = 1:numImages
    img(nn,1) = AbsorptionImage;
end


for jj = 1:numImages
    %
    % Copy immutable properties
    %
    img(jj).constants.copy(imgconsts);
    img(jj).raw.copy(raw(jj));
    img(jj).setClouds(size(roiRow,1));
    for nn = 1:numel(img(jj).clouds)
%         img(jj).clouds(nn).fitdata.set('roirow',roiRow(nn,:),'roiCol',roiCol(nn,:),...
%             'roiStep',roiStep,'fittype',fittype,'lb',lb,'ub',ub,'guess',guess);
        img(jj).clouds(nn).fitdata.set('roirow',roiRow(nn,:),'roiCol',roiCol(nn,:),...
            'roiStep',roiStep,'fittype',fittype,'method','y');
    end
    %
    % Create image
    %
    if size(img(jj).raw.images,3) == 2
        img(jj).makeImage;
    elseif size(img(jj).raw.images,3) == 3
        img(jj).makeImage([1,2,3]);
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
            figure(3);clf;
            img(jj).plotAllData(dispOD,plotROI);
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
            figure(3);clf;
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


