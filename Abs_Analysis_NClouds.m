function cloud = Abs_Analysis_NClouds(varargin)

atomType = 'Rb87';
tof = 35e-3;

col1 = 'b.-';
col2 = 'r--';
dispOD = [0,2.5];
plotOpt = 0;
plotROI = 0;

%% In-trap imaging ROI


% roiRow = repmat([575,725],3,1);
roiRow = [400,660;
          900,660];
roiCol = [60,1200;
          0,60];
      
      

% plotROI = {[575,675],[0,175]};

%% Second Spot Imaging ROI 
% roiRow = [175,450;
%           450,720];
% roiCol = repmat([550,850],2,1);

% plotROI = {[200,850],[550,850]};



%% Fit data
fitdata = AtomCloudFit('roiRow',roiRow(1,:),...
                       'roiCol',roiCol(1,:),...
                       'roiStep',1,...
                       'fittype','tf2d');    %Options: none, gauss1d, twocomp1d, bec1d, gauss2d, twocomp2d, bec2d

imgconsts = AtomImageConstants(atomType,'tof',tof,'detuning',0,...
            'pixelsize',6.45e-6,'magnification',0.99,...
            'freqs',2*pi*[53,53,25],'exposureTime',15e-6,...
            'polarizationcorrection',1.5,'satOD',5);

directory = 'E:\RawImages\2021';
% directory = 'Z:';

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
plotOpt = plotOpt || numImages==1;

numClouds = size(roiCol,1);
cloud(numImages,numClouds) = AbsorptionImage;

exRegion = cell(numClouds,2);
for mm = 1:numClouds
    exRegion(mm,:) = {roiRow(mm,1):roiRow(mm,2),roiCol(mm,1):roiCol(mm,2)};
end

for nn = 1:numImages
    for mm = 1:numClouds
        cloud(nn,mm).constants.copy(imgconsts);
        cloud(nn,mm).fitdata.copy(fitdata);
        cloud(nn,mm).fitdata.set('roirow',roiRow(mm,:),'roicol',roiCol(mm,:));
        cloud(nn,mm).raw.copy(raw(nn));
        cloud(nn,mm).makeImage;
        cloud(nn,mm).fit('method','x');
            
        %% Plotting
        if plotOpt
            if numImages == 1
                figure(3+mm-1);clf;
                cloud(nn,mm).plotAllData(dispOD,col1,col2,plotROI);
                pause(0.1);
            else
                if nn == 1
                    figure(3);clf;
                    dimSubPlot=ceil(sqrt(numImages));
                end

                figure(3);
                subplot(dimSubPlot,dimSubPlot,nn);
                cloud(nn,mm).plotAbsData(dispOD,plotROI);
            end
        end

        %% Print summaries
        [labelStr,numStr]=cloud(nn,mm).labelOneROI;
        if nn == 1
            disp(labelStr);
        end
        disp(numStr);
    end
    
end

end


