function img = Abs_Analysis_DualState_RT(varargin)

atomType = 'Rb87';
tof = evalin('base', 'opt.tof');
detuning = evalin('base', 'opt.detuning');
img_type = 'drop 1';
fittype = 'gauss2d';

dispOD = [0,.25];
plotOpt = 1;
plotROI = {[100,2048],[1,512]};
useFilt = 1;
filtWidth = 50e-6;

%% ROI for the F = 2 image
roiRow{1} = 1200 + 100*sort([-5,-3;-3,-1;-1,1;1,3;3,5]);
roiCol{1} = repmat([150,350],size(roiRow{1},1),1);
roiStep{1} = 3;
%% ROI for F = 1 image
roiRow{2} = 1200 + 200*sort([-1,1.2;1.2,3;3,5]);
% roiRow{2} = [1100,1500];
roiCol{2} = repmat([125,375],size(roiRow{2},1),1);
roiStep{2} = 3;

%% Imaging parameters
freqs = get_trap_freq(2,2);
imgconsts = AtomImageConstants(atomType,'tof',tof,'detuning',detuning,...
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
args = parse_arguments(varargin{:});
%
% This loads the raw image sets
%
raw = BinaryImageData.loadImageSets('directory',directory,'rotation',rot,args{:});

numImages = numel(raw);
plotOpt = plotOpt || numImages == 1;    %This always enables plotting if only one image is analyzed

img = AbsorptionImage.empty;
% for nn = 1:numImages
%     img(nn,1) = AbsorptionImage(BinaryImageData);
% end

for kk = 1:2
    for jj = 1:numImages
        %
        % Copy immutable properties
        %
        img(jj,kk) = AbsorptionImage(BinaryImageData);
        img(jj,kk).constants.copy(imgconsts);
        img(jj,kk).raw.copy(raw(jj));
        if kk == 1
            img(jj,kk).raw.images = img(jj,kk).raw.images(:,:,[1,3,4]); % atoms in F = 2 manifold
        else
            img(jj,kk).raw.images = img(jj,kk).raw.images(:,:,[2,3,4]); % atoms in F = 1 manifold
        end
        img(jj,kk).setClouds(size(roiRow{kk},1));
        for nn = 1:numel(img(jj,kk).clouds)
            imgsize = size(img(jj,kk).raw.images);
            img(jj,kk).clouds(nn).fitdata.set('imgsize',imgsize(1:2),'roirow',roiRow{kk}(nn,:),'roiCol',roiCol{kk}(nn,:),...
                'roiStep',roiStep{kk},'fittype',fittype,'method','y');
        end
        %
        % Create image
        %
        img(jj,kk).makeImage([1,2,3]);
        if useFilt
            img(jj,kk).butterworth2D(filtWidth);
        end
        %
        % Fit clouds
        %
        img(jj,kk).fit;
            
        %% Total Atom Number in ROIs
        
        for kk = 1:size(img,2)
            NsumArray = {img(kk).clouds.Nsum};
            NArray = {img(kk).clouds.N};
            tempData(kk).Ntotal = sum(cell2mat(NArray));
            tempData(kk).Nsumtotal = sum(cell2mat(NsumArray));

            tempData(kk).mfFrac_sum = cell2mat(NArray)/tempData(kk).Ntotal;
            tempData(kk).mfFrac_fit = cell2mat(NsumArray)/tempData(kk).Nsumtotal;
            if kk == size(img,2)
                ManifoldFrac_sum = tempData(kk).Nsumtotal/sum(tempData(kk).Nsumtotal);
                ManifoldFrac_fit = tempData(kk).Ntotal/sum(tempData(kk).Ntotal);
            end
        end

        %% Plotting
        if plotOpt
            figure(kk);clf;
            img(jj,kk).plotAllData(dispOD,plotROI);
            h = gcf;
            ax = h.Children(end);             
            title(ax,sprintf('%s, F = %d',ax.Title.String,3 - kk));            
            children = get(h, 'Children');
            subplot_handle = children(4);
            axes(subplot_handle);
            for ii = 1:size(tempData(kk).mfFrac_sum,2)
                text(plotROI{2}(1), roiRow{kk}(ii)-30,sprintf('P=%g',round(tempData(kk).mfFrac_sum(ii),3)),'Color','w')                
            end
        end
        
        %% Print summaries
%         [labelStr,numStr] = img(jj,kk).labelOneROI;
%         if jj == 1 && kk == 1
%             disp(labelStr);
%         end
%         disp(numStr);


    end
end

Nfit = [img(1).get('N'),img(2).get('N')];
Pfit = Nfit./sum(Nfit);
Nsum = [img(1).get('Nsum'),img(2).get('Nsum')];
Psum = Nsum./sum(Nsum);
fprintf('% 8s|% 10s|% 10s|% 10s|% 10s\n','Image','Nfit','Nsum','Pfit','Psum');
fprintf('% 6.1f|% 10.0f|% 10.0f|% 10.3f|% 10.3f\n',img(1).raw.getImageNumbers,sum(Nfit),sum(Nsum),0,0);
for nn = 1:numel(Nfit)
    if nn <= 5
        s = sprintf('2,%.0f',nn - 3);
    else
        s = sprintf('1,%.0f',nn - 7);
    end
    fprintf('% 8s|% 10.0f|% 10.0f|% 10.3f|% 10.3f\n',s,Nfit(nn),Nsum(nn),Pfit(nn),Psum(nn));
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
