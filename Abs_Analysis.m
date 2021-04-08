function cloud = Abs_Analysis(varargin)

atomType = 'Rb87';
tof = 15e-3;

col1 = 'b.-';
col2 = 'r--';
dispOD = [0,0.5];
plotOpt = 0;
plotROI = 0;

fitdata = AtomCloudFit('roiRow',[301,601],...
                       'roiCol',[301,601],...
                       'roiStep',10,...
                       'fittype','gauss2d');    %Options: none, gauss1d, twocomp1d, bec1d, gauss2d, twocomp2d, bec2d, sum

imgconsts = AtomImageConstants(atomType,'exposureTime',100e-6,...
            'pixelsize',5.5e-6,'magnification',0.25,...
            'freqs',2*pi*[40,23,8],'detuning',8,...
            'polarizationcorrection',1.5,'satOD',3);

directory = 'D:\RawImages\2020\12December\';
% directory = 'Z:';

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
plotOpt = plotOpt || numImages==1;

cloud(numImages,1) = AbsorptionImage;

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


