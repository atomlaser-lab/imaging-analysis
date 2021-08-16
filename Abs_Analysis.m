function cloud = Abs_Analysis(varargin)

atomType = 'Rb87';
% tof = varargin{2};
tof = 217e-3;%camera two

col1 = 'b.-';
col2 = 'r--';
dispOD = [0,.25];
plotOpt = 0;
plotROI = 0;
useFilt = 1;
filtWidth = 25e-6;
%% Imaging ROI first spot
fitdata = AtomCloudFit('roiRow',[10,1300],...
                       'roiCol',[50,1300],...
                       'roiStep',10,...
                       'fittype','tf2d');   


% fitdata = AtomCloudFit('roiRow',[400,800],...
%                        'roiCol',[50,370],...
%                        'roiStep',2,...
%                        'fittype','tf2d'); 

% fitdata = AtomCloudFit('roiRow',[500,720],...
%                        'roiCol',[0,250],...
%                        'roiStep',2,...
%                        'fittype','tf2d');    %Options: none, gauss1d, twocomp1d, bec1d, gauss2d, twocomp2d, bec2d

%% Imaging Second spot
% fitdata = AtomCloudFit('roiRow',[00,1000],...
%                        'roiCol',[550,850],...
%                        'roiStep',1,...
%                        'fittype','tf2d'); 

% fitdata = AtomCloudFit('roiRow',[180,700],...
%                        'roiCol',[550,850],...
%                        'roiStep',1,...
%                        'fittype','tf2d'); 

% fitdata = AtomCloudFit('roiRow',[0,1000],...
%                        'roiCol',[550,850],...
%                        'roiStep',2,...
%                        'fittype','none'); 

%% Imaging parameters
imgconsts = AtomImageConstants(atomType,'tof',tof,'detuning',0,...
            'pixelsize',6.45e-6,'magnification',0.99,...
            'freqs',2*pi*[53,53,25],'exposureTime',14e-6,...
            'polarizationcorrection',1.5,'satOD',5);

directory = 'E:\RawImages\2021';

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

cloud = AbsorptionImage.empty;
for nn = 1:numImages
    cloud(nn,1) = AbsorptionImage;
end


for jj = 1:numImages
    %
    % Copy immutable properties
    %
    cloud(jj).constants.copy(imgconsts);
    cloud(jj).fitdata.copy(fitdata);
    cloud(jj).raw.copy(raw(jj));
    %
    % Create image and fit
    %
    if size(cloud(jj).raw.images,3) == 2
        cloud(jj).makeImage;
    elseif size(cloud(jj).raw.images,3) == 3
        cloud(jj).makeImage([1,2,3]);
    else
        error('Not sure what to do here');
    end
    if useFilt
        cloud(jj).butterworth2D(filtWidth);
    end
    cloud(jj).fit('method','y');
        
    %% Plotting
    if plotOpt
        if numImages == 1
            %
            % Plot absorption data and marginal distributions when there is only 1 image
            %
            figure(3);clf;
            cloud(jj).plotAllData(dispOD,col1,col2,plotROI);
        else
            %
            % Plot only the absorption data in a grid when there is more than one image
            %
            if jj == 1
                figure(3);clf;
                dimSubPlot=ceil(sqrt(numImages));
            end
            
            figure(3);
            subplot(dimSubPlot,dimSubPlot,jj);
            cloud(jj).plotAbsData(dispOD,plotROI);
        end
    end
    
    %% Print summaries
    [labelStr,numStr]=cloud(jj).labelOneROI;
    if jj == 1
        disp(labelStr);
    end
    disp(numStr);
    
end

end


