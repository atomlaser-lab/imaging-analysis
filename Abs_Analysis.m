function img = Abs_Analysis(varargin)

atomType = 'Rb87';
tof = evalin('base', 'opt.tof');
detuning = evalin('base', 'opt.detuning');

dispOD = [0,0.1];
plotOpt = 1;
plotROI = 0;
useFilt = 1;
filtWidth = 50e-6;
% detuning = 8;

%% Imaging first spot
% roiRow = [1600 + 250*[-1,1];1800 + 500*[-1,1]];
roiRow = [1200,1500;1500,1750];
roiCol = repmat(1050 + 150*[-1,1],size(roiRow,1),1);
% roiRow = [600,1000];
% roiCol = 1050 + 150*[-1,1];
roiStep = 2;
fittype = 'gauss2d'; %gauss2d
% fittype  = '2comp1d';
img_type = 'drop 1';

%% Imaging parameters
freqs = get_trap_freq(2,2);
imgconsts = AtomImageConstants(atomType,'tof',tof,'detuning',0,...
            'pixelsize',5.5e-6,'magnification',1.0285,...
            'freqs',2*pi*freqs,'exposureTime',5e-6,...
            'polarizationcorrection',1,'satOD',5);
 
if strcmpi(img_type,'drop 1')
    rot = 180;
    imgconsts.magnification = 1.0285;
    imgconsts.exposureTime = 5e-6;
elseif strcmpi(img_type,'drop 2')
    rot = 90;
    imgconsts.magnification = 1.0801;
    imgconsts.exposureTime = 50e-6;
end       
       

%% Load raw data
directory = 'E:\labview-images';
% directory = 'E:\SpatialFringes_data\tides-take2';
args = parse_arguments(varargin{:});
%
% This loads the raw image sets
%
raw = BinaryImageData.loadImageSets('directory',directory,'rotation',rot,args{:});

numImages = numel(raw);
plotOpt = plotOpt || numImages == 1;    %This always enables plotting if only one image is analyzed

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
    img(jj).setClouds(size(roiRow,1));
    for nn = 1:numel(img(jj).clouds)
%         img(jj).clouds(nn).fitdata.set('roirow',roiRow(nn,:),'roiCol',roiCol(nn,:),...
%             'roiStep',roiStep,'fittype',fittype,'lb',lb,'ub',ub,'guess',guess);
        imgsize = size(img(jj).raw.images);
        img(jj).clouds(nn).fitdata.set('imgsize',imgsize(1:2),'roirow',roiRow(nn,:),'roiCol',roiCol(nn,:),...
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
    img(jj).fit;
        
    %% Plotting
    if plotOpt
        if numImages == 1
            %
            % Plot absorption data and marginal distributions when there is only 1 image
            %
            figure(1);clf;
            img(jj).plotAllData(dispOD,plotROI);
            h = gcf;
            ax = h.Children(end);
%             title(ax,sprintf('%s, %.2f um',ax.Title.String,img.clouds.pos(2)*1e6 - 3.3048e3));
            title(ax,sprintf('%s, %.2e',ax.Title.String,img.clouds.N));
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

function args = parse_arguments(varargin)

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
    args = {'files','last','index',idx,'len',3};
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

end
