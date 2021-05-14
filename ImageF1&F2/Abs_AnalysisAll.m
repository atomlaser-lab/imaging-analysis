function cloud = Abs_AnalysisAll(varargin)
% function cloud = Abs_Analysis('last')
resetfolder=pwd;

atomType = 'Rb87';
tof = 13e-3;  %set expansion time here

col1 = 'b.-';
col2 = 'r--';
dispOD = [0,0.4];
plotOpt = 0;
plotROI = 0;

fitdata = AtomCloudFit('roiRow',[1,512],...   %ensure this matches the size in RawImageData
                       'roiCol',[201,800],...%1216
                       'roiStep',5,...
                       'fittype','gauss2d');    %Options: none, gauss1d, twocomp1d, bec1d, gauss2d, twocomp2d, bec2d, sum

imgconsts = AtomImageConstants(atomType,'exposureTime',100e-6,'tof',tof,...
            'pixelsize',5.5e-6,'magnification',0.25,...
            'freqs',2*pi*[40,23,8],'detuning',0,... %set detuning here
            'polarizationcorrection',1.5,'satOD',3);

directory = 'D:\RawImages\2020\12December\';
% directory = 'Z:';

%% Load raw data
%
% Check that arguments appear as name/value pairs
%
if mod(nargin,2) ~= 0
    error('Arguments must occur as name/value pairs!');
end
%
% This loads the raw image sets
%
raw = RawImageData.loadImageSets('directory',directory,varargin{:});


numImages = numel(raw);
plotOpt = plotOpt || numImages==1;
%
% This creates an array of AbsorptionImage objects, and it has to be done
% this way because of some strange MATLAB problem where it doesn't create
% clean copies of the objects.
%
for nn = 1:numImages
    cloud(nn,1) = AbsorptionImage;
    cloud(nn,2) = AbsorptionImage;
end

for jj = 1:numImages
    for mm = 1:size(cloud,2)

        cloud(jj,mm).constants.copy(imgconsts);
        cloud(jj,mm).fitdata.copy(fitdata);
        cloud(jj,mm).raw.copy(raw(jj));
        cloud(jj,mm).makeImage([mm,3]); %3 is always the background image
        cloud(jj,mm).fit('method','y');

        %% Plotting
        if plotOpt
            %
            % If plotting is enabled, plot images
            %
            if numImages == 1
                %
                % If only one image set is being analyzed, plot OD with
                % marginal distributions and image information
                %
                if mm == 1
                    %
                    % Plot first image in set
                    %
                    figure(3);clf;
                    cloud(jj,mm).plotAllData(dispOD,col1,col2,plotROI);   
                else
                    %
                    % Plot second image in set
                    %
                    figure(4);clf;
                    cloud(jj,mm).plotAllData(dispOD,col1,col2,plotROI);  
                end
            else
                %
                % If there is more than 1 image set being analyzed, plot
                % only the OD as subplots
                %
                if mm == 1
                    figure(3);clf;
                else
                    figure(4);clf;
                end
                if jj == 1
                    dimSubPlot=ceil(sqrt(numImages));
                end
                figure(3);
                subplot(dimSubPlot,dimSubPlot,jj);
                cloud(jj,mm).plotAbsData(dispOD,plotROI);
            end
        end

        %% Print summaries
        [labelStr,numStr]=cloud(jj,mm).labelOneROI;
        if jj == 1 && mm == 1
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

end


