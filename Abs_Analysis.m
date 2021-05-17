function cloud = Abs_Analysis(varargin)

atomType = 'Rb87';
tof = 216.5e-3;

col1 = 'b.-';
col2 = 'r--';
dispOD = [0,.25];
plotOpt = 0;
plotROI = 0;

% fitdata = AtomCloudFit('roiRow',[450,800],...
%                        'roiCol',[10,450],...
%                        'roiStep',2,...
%                        'fittype','tf2d');    %Options: none, gauss1d, twocomp1d, bec1d, gauss2d, twocomp2d, bec2d
fitdata = AtomCloudFit('roiRow',[150,550],...
                       'roiCol',[550,850],...
                       'roiStep',1,...
                       'fittype','tf1d'); 

imgconsts = AtomImageConstants(atomType,'tof',tof,'detuning',0,...
            'pixelsize',6.45e-6,'magnification',0.99,...
            'freqs',2*pi*[53,53,25],'exposureTime',15e-6,...
            'polarizationcorrection',1.5,'satOD',5);

directory = 'C:\Users\Ryan\MATLAB\spatial-fringes-analysis\images';

%% Load raw data
if nargin == 0 || (nargin == 1 && strcmpi(varargin{1},'last')) || (nargin == 2 && strcmpi(varargin{1},'last') && isnumeric(varargin{2}))
    %
    % If no input arguments are given, or the only argument is 'last', or
    % if the arguments are 'last' and a numeric array, then load the last
    % image(s).  In the case of 2 arguments, the second argument specifies
    % the counting backwards from the last image
    %
    args = {'files','last','index',varargin{2}};
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

cloud = AbsorptionImage;
if numImages > 1
    for nn = 2:numImages
        cloud(nn,1) = AbsorptionImage;
    end
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
    cloud(jj).makeImage;
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


