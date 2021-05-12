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

imgconsts = AtomImageConstants(atomType,'exposureTime',15e-6,...
            'pixelsize',6.45e-6,'magnification',0.99,...
            'freqs',2*pi*[53,53,25],'detuning',0,...
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
    if nargin == 2
        idx = varargin{2};
        raw(numel(idx),1) = RawImageData;
        for nn = 1:numel(idx)
            raw(nn).load('filenames','last','directory',directory,'index',idx(nn));
        end
    else
        raw = RawImageData('filenames','last','directory',directory);
    end
else
    %
    % Otherwise, parse arguments as name/value pairs for input into
    % RawImageData
    %
    if mod(nargin,2) ~= 0
        error('Arguments must occur as name/value pairs!');
    end
    
    filenames = {};
    index = 1;
    len = 2;
    for nn = 1:2:nargin
        v = varargin{nn+1};
        switch lower(varargin{nn})
            case {'filenames','files'}
                filenames = v;
            case 'directory'
                directory = v;
            case {'index','idx'}
                index = v;
            case {'length','len'}
                len = v;
        end
    end
    
    if ~isempty(filenames)
        numImages = numel(filenames);
        raw(numImages,1) = RawImageData;
        for mm = 1:numImages
            raw(mm).load('filenames',filenames{mm},'directory',directory);
        end
    else
        numImages = numel(index);
        raw(numImages,1) = RawImageData;
        for mm = 1:numImages
            raw(mm).load('filenames','last','directory',directory,'index',index(mm),'length',len);
        end
    end
    
end

numImages = numel(raw);
plotOpt = plotOpt || numImages==1;

cloud = AbsorptionImage;
if numImages > 1
    for nn = 2:numImages
        cloud(nn,1) = AbsorptionImage;
    end
end

for jj = 1:numImages

    cloud(jj).constants.copy(imgconsts);
    cloud(jj).fitdata.copy(fitdata);
    cloud(jj).raw.copy(raw(jj));
    cloud(jj).makeImage;
%     cloud(jj).fitdata.makeFitObjects(cloud(jj).x,cloud(jj).y,cloud(jj).image);
    cloud(jj).fit([],tof,'y');
%     cloud(jj).fit([],tof,'y',3);
        
    %% Plotting
    if plotOpt
        if numImages == 1
            figure(3);clf;
            cloud(jj).plotAllData(dispOD,col1,col2,plotROI);
        else
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


