function cloud = Abs_Analysis(varargin)
% function cloud = Abs_Analysis('last')
resetfolder=pwd;

atomType = 'Rb87';
tof = 13e-3;  %set expansion time here

col1 = 'b.-';
col2 = 'r--';
dispOD = [0,0.4];
plotOpt = 1;
plotROI = 0;

fitdata = AtomCloudFit('roiRow',[101,951],...
                       'roiCol',[11,1001],...
                       'roiStep',5,...
                       'fittype','gauss2d');    %Options: none, gauss1d, twocomp1d, bec1d, gauss2d, twocomp2d, bec2d, sum

imgconsts = AtomImageConstants(atomType,'exposureTime',100e-6,...
            'pixelsize',5.5e-6,'magnification',0.25,...
            'freqs',2*pi*[40,23,8],'detuning',0,... %set detuning here
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
    
    %% save output to text file
    %make directory
    fpathfull = [mfilename('fullpath'),'.m'];
    [fpath,fname,fext] = fileparts(fpathfull);
    dstr = datestr(datetime,'YYYY\\mm\\dd\\hh_MM_ss');
    subfolder = 'data-archive';
    dirname = sprintf('%s\\%s\\%s',fpath,subfolder,datestr(datetime,'YYYY\\mm\\dd'));
    if ~isfolder(dirname)
        mkdir(dirname);
    end
    %open directory
    cd(dirname)
    %make text file with name=time
    textname=datestr(datetime,'hh_MM_ss');
    fileID=fopen(textname,'w');
    blah='Image  |  x width/um  |  y width/um  |  Natoms  |  BEC  |  PeakOD  |  T/nk  |  PSD\n';
    
    %write on file so that it contains data
    fprintf(fileID,blah);
    fprintf(fileID,numStr);
    cd(resetfolder);
end

end


