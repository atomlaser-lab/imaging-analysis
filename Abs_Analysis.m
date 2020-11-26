function cloud = Abs_Analysis(varargin)

atomType = 'Rb87';
tof = 28e-3;

col1 = 'b.-';
col2 = 'r--';
maxOD = 2;
plotOpt = 0;

fitdata = AtomCloudFit('roiRow',[50,980],...
                       'roiCol',[50,1330],...
                       'roiStep',1,...
                       'fittype','gauss2d');    %Options: gauss1d, twocomp1d, bec1d, gauss2d, gauss2dangle

imgconsts = AtomImageConstants(atomType,'exposureTime',30e-6,...
            'pixelsize',6.45e-6,'magnificantion',0.97,...
            'freqs',2*pi*[220,230,100],'photonspercount',1/0.17);

directory = 'E:\RawImages\2019\01Jan\top';

%% Load raw data
if ischar(varargin{1}) && strcmpi(varargin{1},'last')
    if numel(varargin) == 1
        raw = RawImageData('filenames','last','directory',directory);
    elseif numel(varargin) == 2
        raw(numel(varargin{2}),1) = RawImageData;
        for nn = 1:numel(varargin{2})
            raw(nn).load('filenames','last','directory',directory,'index',varargin{2}(nn));
        end
    else
        error('Unsupported argument list');
    end
elseif mod(numel(varargin),2) ~= 0
    error('Must specify even numbers of file names');
else
    raw(round(numel(varargin)/2),1) = RawImageData;
    for nn = 1:2:numel(varargin)
        raw(nn).load('length',2,'filenames',varargin{nn:(nn+1)},'directory',directory);
    end
end

numImages = numel(raw);
plotOpt = plotOpt || N==1;

cloud(numImages,1) = AbsorptionImage;

for jj = 1:numImages

    cloud(jj).constants.copy(imgconsts);
    cloud(jj).fitdata.copy(fitdata);
    cloud(jj).makeImage;
    cloud(jj).fitdata.makeFitObjects(cloud(jj).x,cloud(jj).y,cloud(jj).image);
    cloud(jj).fit([],tof,'xy');
        
    %% Plotting
    if plotOpt
        if numImages == 1
            figure(3);clf;
            cloud(jj).plotAllData(maxOD,col1,col2);
        else
            if jj == 1
                figure(3);clf;
                dimSubPlot=ceil(sqrt(numImages));
            end
            
            figure(3);
            subplot(dimSubPlot,dimSubPlot,jj);
            cloud(jj).plotAbsData(maxOD);
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


