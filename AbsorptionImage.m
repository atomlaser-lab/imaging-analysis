classdef AbsorptionImage < handle
    %ABSORPTIONIMAGE Defines a class that describes an absorption image.
    %
    %   Properties of the absorption image are things like the image
    %   itself, the x and y coordinates, the number of atoms, widths,
    %   temperatures, etc.
    %
    %   Three "constant" handle properties are used - RAW, CONSTANTS, and
    %   FITDATA of type RAWIMAGEDATA, ATOMIMAGECONSTANTS, and ATOMCLOUDFIT,
    %   respectively.  These describe the raw image files, the constants
    %   used for extracted useful information about the sample (like image
    %   detuning, pixel size, etc), and the fit function, data, and
    %   methods, respectively
    
    properties
        %
        % Image properties
        %
        x               %x position
        y               %y position
        ODraw           %raw absorption image: -log(imgWithAtoms/imgWithoutAtoms)
        ODcorr          %corrected absorption image
        imgidxs         %Indicies in raw image data corresponding to data used for this instance
    end
    
    properties(SetAccess = protected)
        clouds          %Array of AtomCloud objects describing all atom clouds in image
    end

    properties(SetAccess = immutable)
        raw             %RAWIMAGEDATA object describing the raw image data
        constants       %ATOMIMAGECONSTANTS object describing the conditions under which the measurement was done
    end

    methods
        function self = AbsorptionImage(varargin)
            %ABSORPTIONIMAGE Creates and ABSORPTIONIMAGE object
            %
            %   C = ABSORPTIONIMAGE() Creates object C with constant
            %   properties set to default values
            %
            %   C = ABSORPTIONIMAGE(RAW) creates object C with raw image
            %   data RAW
            %
            %   C = ABSORPTIONIMAGE(RAW,CONSTANTS) creates object C with
            %   raw image data RAW and constants CONSTANTS
            %
            %   C = ABSORPTIONIMAGE(RAW,CONSTANTS,FITDATA) creates object C
            %   with raw image data RAW, constants CONSTANTS, and fit data
            %   FITDATA
            if nargin < 1
                self.raw = RawImageData;
            else
                self.raw = varargin{1};
            end
            if nargin < 2
                self.constants = AtomImageConstants;
            else
                self.constants = varargin{2};
            end
        end
        
        function self = setClouds(self,v)
            %SETCLOUDS Sets the number of clouds or the CLOUDS property
            %
            %   SELF = SELF.SETCLOUDS(V) where V is an integer creates V
            %   blank ATOMCLOUD objects
            %
            %   SELF = SELF.SETCLOUDS(V) where V is an array of ATOMCLOUD
            %   objects sets the internal CLOUDS property to V
            if isnumeric(v)
                self.clouds = AtomCloud.empty;
                for nn = 1:v
                    self.clouds(nn,1) = AtomCloud(self.constants);
                end
            elseif isa(v,'AtomCloud')
                self.clouds = v;
            else
                error('Input must be either an integer or an array of AtomCloud objects');
            end
        end
                

        function self = makeImage(self,imgIndexes)
            %MAKEIMAGE Creates the uncorrected and corrected absorption
            %images
            %
            %   C = C.MAKEIMAGE() creates the images for object C based on
            %   constant properties
            %
            %   C = C.MAKEIMAGE(IMGINDEXES) creates the images for object C
            %   assuming that the image with atoms is the raw image with
            %   index IMGINDEXES(1) and that the image without atoms is the
            %   raw image with index IMGINDEXES(2)
            
            % Shorten names of constants and raw data
            c = self.constants;
            r = self.raw;
            % Compute the saturation intensity in camera counts
            Nsat = c.satN;
            %
            % Select "image with atoms" and "image without atoms". Either
            % use the default (first and second raw images) or use the
            % user-supplied indices.
            %
            if nargin < 2 && isempty(self.imgidxs)
                self.imgidxs = [1,2];
            elseif nargin >= 2
                self.imgidxs = imgIndexes;
            end
            if numel(self.imgidxs) == 2
                imgWithAtoms = r.images(:,:,self.imgidxs(1));
                imgWithoutAtoms = r.images(:,:,self.imgidxs(2));
            elseif numel(self.imgidxs) == 3
                imgWithAtoms = r.images(:,:,self.imgidxs(1)) - r.images(:,:,self.imgidxs(3));
                imgWithoutAtoms = r.images(:,:,self.imgidxs(2)) - r.images(:,:,self.imgidxs(3));
            else
                error('Number of images to use is set to %d. Not sure what to do here...',numel(self.imgidxs));
            end
            %
            % Create the uncorrected optical depth map and get peak OD.
            % Set NaNs and Infs to zero
            %
            tmp = real(-log(imgWithAtoms./imgWithoutAtoms));
            tmp(isnan(tmp) | isinf(tmp)) = 0;
            self.ODraw = tmp;
%             self.peakOD = max(max(self.ODraw));
            %
            % Correct for finite saturation OD
            %
            if ~isinf(c.satOD)
                ODmod = real(log((1-exp(-c.satOD))./(exp(-self.ODraw)-exp(-c.satOD))));
            else
                ODmod = self.ODraw;
            end
            %
            % Correct for polarization and intensity saturation - this is
            % the corrected optical depth map
            %
            self.ODcorr = c.polarizationCorrection*ODmod + (1 - exp(-ODmod)).*imgWithoutAtoms./Nsat;
            %
            % Create x and y vectors based on pixel size and magnification
            %
            self.x = (c.pixelSize/c.magnification)*(1:size(self.ODraw,2));
            self.y = (c.pixelSize/c.magnification)*(1:size(self.ODraw,1));
        end
        
        function self = butterworth2D(self,spatialWidth,order,filterType)
            %BUTTERWORTH2D Applies a 2D Butterworth filter to the corrected
            %   and raw OD images
            %
            %   C = C.BUTTERWORTH2D(WIDTH) Applies a 4th order Butterworth
            %   filter with spatial widths defined by WIDTH.  If WIDTH has
            %   2 elements, the x and y filter widths are defined by those
            %   elements.
            %
            %   C = C.BUTTERWORTH2D(__,ORDER) Applies a Butterworth filter
            %   of order ORDER.
            %
            %   C = C.BUTTERWORTH2D(__,TYPE) Applies a filter that is
            %   either low pass (TYPE = 'low') or high pass (TYPE = 'high')
            imgfft = fftshift(fft2(self.ODcorr));
            imgrawfft = fftshift(fft2(self.ODraw));
            kx = 1/(2*diff(self.x(1:2)))*linspace(-1,1,numel(self.x));
            ky = 1/(2*diff(self.y(1:2)))*linspace(-1,1,numel(self.y));
            [KX,KY] = meshgrid(kx,ky);
            if numel(spatialWidth) == 1
                spatialWidth(2) = spatialWidth;
            end
            K2 = (spatialWidth(1)*KX).^2 + (spatialWidth(2)*KY).^2;
            if nargin < 3 || isempty(order)
                order = 4;
            end
            if nargin < 4 || strcmpi(filterType,'low')
                F = (1+K2.^order).^(-1);
            elseif strcmpi(filterType,'high')
                F = K2.^order./(1+K2.^order);
            else
                error('Filter type %s unknown!',filterType);
            end
            self.ODcorr = real(ifft2(ifftshift(F.*imgfft)));
            self.ODraw = real(ifft2(ifftshift(F.*imgrawfft)));
        end

        function self = fit(self)
            %FIT Performs fits to all ATOMCLOUD objects CLOUDS
            %
            for nn = 1:numel(self.clouds)
                self.clouds(nn).fitdata.makeFitObjects(self.x,self.y,self.ODcorr);
                self.clouds(nn).fit;
            end
        end
        
        function self = jointFit(self,whichClouds,lbi,ubi,gi)
            %JOINTFIT Performs a joint 2D TF fit to specified clouds.
            %
            %   SELF = SELF.JOINTFIT(WHICHCLOUDS) Performs a joint TF fit
            %   using the clouds specified by the vector WHICHCLOUDS, which
            %   indicates the indices of the CLOUDS property to use.  The
            %   joint fit estimates parameter bounds and values, but
            %   each value is overwritten by the corresponding value in the
            %   CLOUDS.FITDATA obejct.
            %
            %   SELF = SELF.JOINTFIT(__,LB,UB,G) uses global lower and
            %   upper bounds LB and UB and global guess G (all are
            %   CLOUDPARAMETERS objects) to perform fit.  If not using,
            %   set to [].  Any parameter in these objects set to [] will
            %   not be used for the fit constraints.
            
            %
            % Select clouds and initialize bounds/guesses
            %
            c = self.clouds(whichClouds);
            lb = CloudParameters.empty;
            ub = CloudParameters.empty;
            g = CloudParameters.empty;
            lbarray = [];
            ubarray = [];
            garray = [];
            order = {'becamp','posx','becwidthx','posy','becwidthy'};
            xmin = Inf;ymin = Inf;
            xmax = 1;ymax = 1;
            xstep = 1;ystep = 1;
            for nn = 1:numel(c)
                %
                % Make fit objects, and overwrite the upper and lower
                % bounds on the BEC width
                %
                c(nn).fitdata.makeFitObjects(self.x,self.y,self.ODcorr);
                %
                % Estimate bounds and initial guesses, then replace with
                % user-supplied guesses as appropriate
                %
                [lb(nn),ub(nn)] = AtomCloudFit.guessBounds2D(c(nn).fitdata.x,c(nn).fitdata.y,c(nn).fitdata.image);
                g(nn) = AtomCloudFit.guessTF2DParams(c(nn).fitdata.x,c(nn).fitdata.y,c(nn).fitdata.image);
                if nargin > 2 && ~isempty(lbi)
                    lb(nn).compare(lbi);
                else
                    lb(nn).compare(c(nn).fitdata.lb);
                end
                if nargin > 3 && ~isempty(ubi)
                    ub(nn).compare(ubi);
                else
                    ub(nn).compare(c(nn).fitdata.ub);
                end
                if nargin > 4 && ~isempty(gi)
                    g(nn).compare(gi);
                else
                    g(nn).compare(c(nn).fitdata.guess);
                end
                lbarray = [lbarray,lb(nn).convert2array(order{:})];
                ubarray = [ubarray,ub(nn).convert2array(order{:})];
                garray = [garray,g(nn).convert2array(order{:})];
                % 
                % Find pixel bounds
                %
                xmin = min(xmin,c(nn).fitdata.roiCol(1));
                xmax = max(xmax,c(nn).fitdata.roiCol(2));
                ymin = min(ymin,c(nn).fitdata.roiRow(1));
                ymax = max(ymax,c(nn).fitdata.roiRow(2));
                xstep = max(xstep,c(nn).fitdata.roiStep(1));
                ystep = max(ystep,c(nn).fitdata.roiStep(2));
            end
            %
            % Create position grid and data
            %
            row = ymin:ystep:ymax;
            col = xmin:xstep:xmax;
            xx = self.x(col);
            yy = self.y(row);
            [X,Y] = meshgrid(xx,yy);
            Z = cat(3,X,Y);
            data = self.ODcorr(row,col);
            %
            % Append background terms
            %
            lbarray = [lbarray,-0.05,-0.02/range(xx),-0.02/range(yy)];
            ubarray = [ubarray,+0.05,+0.02/range(xx),+0.02/range(yy)];
            garray = [garray,0,0,0];
            %
            % Perform fit
            %
            options = AtomCloudFit.getoptions;
            func = @(c,pos) AtomCloudFit.multiFitBEC2D(c,pos);
            params = AtomCloudFit.attemptFit(func,garray,Z,data,lbarray,ubarray,options);
            p0 = [Z(1,1,1),Z(1,1,2)];
            [X,Y] = meshgrid(self.x,self.y);
            Znew = cat(3,X,Y);
            f = func(params,Znew);
            bgWrong = AtomCloudFit.bg2D(params(end-2:end),Znew);
            bg = AtomCloudFit.bg2D(params(end-2:end),Znew,p0);
            f = f - bgWrong + bg;
            res = self.ODcorr - f;
            %
            % Extract parameters
            %
            for nn = 1:numel(c)
                idx = ((nn-1)*5 + 1):(nn*5);
                c(nn).fitdata.params = CloudParameters(0).setFromArray(order,params(idx));
                [row,col] = c(nn).fitdata.makeROIVectors;
                c(nn).fitdata.residuals = res(row,col);
                c(nn).fitdata.bg = bg(row,col);
                c(nn).fitdata.xfit = sum(f(row,col),1);
                c(nn).fitdata.yfit = sum(f(row,col),2);
                %
                % This extracts parameters without doing the single-cloud
                % fit
                %
                c(nn).fit('dofit',false);
            end
            
        end

        function varargout = get(self,varargin)
            %GET Takes a vector of AbsorptionImage objects and returns the
            %properties as vectors.
            %
            %   VARARGOUT = C.GET('NAME1','NAME2',...) takes a vector of
            %   AbsorptionImage objects C and returns in VARARGOUT the vectors
            %   of properties 'NAME1', 'NAME2', etc.  Valid values for 'NAME1',
            %   etc correspond to property names in the AbsorptionImage class.
            if numel(varargin) == 1
                tmp = zeros(numel(self),numel(self(1).clouds));
                for nn = 1:numel(self)
                    for mm = 1:numel(self(1).clouds)
                        tmp(nn,mm) = self(nn).clouds(mm).get(varargin{1});
                    end
                end
                varargout{1} = tmp;
            else
                for nn = 1:numel(self)
                    for mm = 1:numel(varargin)
                        varargout{mm}(nn,:) = self(nn).clouds(1).(varargin{mm});
                    end
                end
            end
        end


        %% Plotting functions
        function self = plotROI(self,row,col)
            %PLOTROI Plots the region-of-interest (ROI) on the current axes
            %
            %   C = C.PLOTROI() plots the ROI stored in ABSORPTIONIMAGE
            %   object C
            %
            %   C = C.PLOTROI(ROW,COL) uses the ROI given by ROW and COL to
            %   plot the ROI
            %
            if nargin > 1
                plot([col(1),col(end),col(end),col(1),col(1)],[row(1),row(1),row(end),row(end),row(1)],'r--');
            else
                for nn = 1:numel(self.clouds)
                    row = self.clouds(nn).fitdata.roiRow;
                    col = self.clouds(nn).fitdata.roiCol;
                    plot([col(1),col(end),col(end),col(1),col(1)],[row(1),row(1),row(end),row(end),row(1)],'r--');
                end
            end
        end

        function str = makeImageSummary(self)
            %MAKEIMAGESUMMARY Creates a string used for summarizing a given
            %absorption image
            %
            %   STR = C.MAKEIMAGESUMMARY() Creates a string STR with the
            %   number of atoms from the fit, the peak OD, the temperature,
            %   and the number of atoms from summing over the ROI
            %
            c = self.clouds(1);
            str{1} = sprintf('N = %1.3g (%d%%)    OD_{peak} = %1.3g    T_{y} = %3.2f uK',...
                c.N,round(c.becFrac*100),c.peakOD,sqrt(prod(c.T))*1e6);
            str{2} = sprintf('Nsum = %1.3g',c.Nsum);
        end

        function self = plotYData(self)
            %PLOTYDATA Plots the marginal y distribution as a vertically
            %oriented plot
            %
            %   C = C.PLOTYDATA plots the data for all atom clouds in
            %   ABSORPTIONIMAGE object C
            %
            for nn = 1:numel(self.clouds)
                %
                % Plots the marginal Y distribution as a
                % vertically-oriented plot
                %
                f = self.clouds(nn).fitdata;
                plot(f.ydata,-f.y,'.-');
                hold on;
                plot(f.yfit,-f.y,'-','linewidth',1);
            end
            if numel(self.clouds) == 1
                %
                % Creates a summary of the Y distribution
                %
                str{1} = sprintf('Gauss_{y} = %3.1f um',self.clouds(1).gaussWidth(2)*1e6);
                str{2} = sprintf('TF_{y} = %3.1f um',self.clouds(1).becWidth(2)*1e6);
                hold off;
                xlabel(str,'fontsize',8);
            end
        end

        function self = plotXData(self)
            %PLOTXDATA Plots the marginal X distribution as a horizontally
            %oriented plot
            %
            %   C = C.PLOTXDATA plots the X distributions for each atomic
            %   cloud in the ABSORPTIONIMAGE object C
            %
            for nn = 1:numel(self.clouds)
                %
                % Plots the marginal X distribution as a
                % horizontally-oriented plot
                %
                f = self.clouds(nn).fitdata;
                h = plot(f.x,f.xdata,'.-');
                hold on;
                plot(f.x,f.xfit,'-','linewidth',1);
            end
            if numel(self.clouds) == 1
                %
                % Creates a summary of the X distribution
                %
                str{1} = sprintf('Gauss_{x} = %3.1f um',self.clouds(1).gaussWidth(1)*1e6);
                str{2} = sprintf('TF_{x} = %3.1f um',self.clouds(1).becWidth(1)*1e6);
                hold off;
                xlabel(str,'fontsize',8);
            end
        end

        function self = plotAllData(self,dispOD,plotROI)
            %PLOTALLDATA Plots the absorption image and the X and Y
            %marginal distributions on the current figure
            %
            %   C = C.PLOTALLDATA(DISPOD) Plots, for object C,
            %   the absorption image and the X and Y distributions on the
            %   current figure. DISPOD is the OD range to display on the
            %   image (using caxis).
            %
            %   C = C.PLOTALLDATA(__,PLOTROI) restricts the image to
            %   PLOTROI.  If PLOTROI is not a cell array, it is interpreted
            %   as a boolean value and restricts the image to the ROI
            %   defined in FITDATA if true.  If it is a cell array, then
            %   this is interpreted as the {ROW,COL} to use for limiting
            %   the image to a given ROI
            %
            if nargin < 3
                plotROI = false;
            end
            %
            % Create the main axes for the image and plot the absorption
            % data
            %
%             subplot(6,6,[2:5 8:11 14:17 20:23]);
            axes('position',[0.3,0.3,0.6,0.65]);
            self.plotAbsData(dispOD,plotROI);
            xlabel(self.makeImageSummary,'fontsize',8);
            %
            % Draw the ROI
            %
            hold on;
            self.plotROI;
            %
            % Plot the Y distribution
            %
%             subplot(6,6,[1 7 13 19]);
            axes('position',[0.075,0.35,0.15,0.6]);
            self.plotYData();
            %
            % Plot the X distribution
            %
%             subplot(6,6, 31:36)
            axes('position',[0.1,0.075,0.8,0.15]);
            self.plotXData();
        end

        function self = plotAbsData(self,dispOD,plotROI)
            %PLOTABSDATA Plots the absorption image data on the current
            %axes
            %
            %   C = C.PLOTABSDATA(DISPOD) Plots the absorption image data
            %   using DISPOD as the caxis limits
            %
            %   C = C.PLOTABSDATA(__,PLOTROI) restricts the image to
            %   PLOTROI.  If PLOTROI is not a cell array, it is interpreted
            %   as a boolean value and restricts the image to the ROI
            %   defined in FITDATA if true.  If it is a cell array, then
            %   this is interpreted as the {ROW,COL} to use for limiting
            %   the image to a given ROI
            %
            if nargin < 3
                plotROI = false;
            end
            %
            % Create the image make it look nice
            %
            imagesc(self.ODraw,dispOD);
            axis equal;
            axis tight;
            colorbar;
            colormap(jet);
            %
            % Set x and y limits to the stored ROI or the user-specified
            % ROI
            %
            if ~iscell(plotROI) && plotROI
                xlim(self.clouds(1).fitdata.roiCol);
                ylim(self.clouds(1).fitdata.roiRow);
            elseif iscell(plotROI)
                xlim(plotROI{2});
                ylim(plotROI{1});
            end
            %
            % Make the title the image number
            %
            imgNums = self.raw.getImageNumbers;
            strTitle = sprintf('Image: %d',imgNums(1));
            title(strTitle,'fontsize',14);
        end

        %% Labelling functions
        function [labelStr,numberStrTotal] = labelOneROI(self)
            %LABELONEROI Creates labels that can be printed on the command
            %line so that the user can see image properties
            %
            %   [LABEL,NUMBER] = C.LABELONEROI() creates a label string
            %   LABEL and a number string with parameters NUMBER
            %
            imgNum = self.raw.getImageNumbers;
            labelCell = {'Image','x width/um','y width/um','Nsum' ,'Nfit' ,'BEC %','PeakOD','T/nk' ,'PSD'};
            formatCell = {'% 5.1f','%0.3e'     ,'%0.3e'     ,'%0.2e','%0.2e','%0.2f','%0.2e' ,'%0.2e','%0.2e'};
            numberStrTotal = {};
            for nn = 1:numel(self.clouds)
                c = self.clouds(nn);
                numberCell = {imgNum(1)+nn/10,c.gaussWidth(1)*1e6,c.gaussWidth(2)*1e6,c.Nsum,c.N,c.becFrac*1e2,c.peakOD,sqrt(prod(c.T))*1e9,c.PSD};
                if nn == 1
                    [labelStr,numberStrTotal{nn}] = self.formatLabel(labelCell,formatCell,numberCell);
                else
                    [~,numberStrTotal{nn}] = self.formatLabel(labelCell,formatCell,numberCell);
                end
            end
            numberStrTotal = strjoin(numberStrTotal,'\n');
        end
        
        
        %% Saving and loading functions
        function s = struct(self)
            %STRUCT Creates a struct from the AbsorptionImage object
            %
            %   S = C.STRUCT() Creates a struct from the AbsorptionImage
            %   object C
            if numel(self) > 1
                for nn = 1:numel(self)
                    s(nn,1) = struct(self(nn)); %#ok<*AGROW>
                end
            else
                s.x = self.x;
                s.y = self.y;
                s.ODraw = self.ODraw;
                s.ODcorr = self.ODcorr;
                s.imgidxs = self.imgidxs;
                s.clouds = self.clouds.struct;
                s.raw = struct(self.raw);
                s.constants = struct(self.constants);
            end
        end
        
        function s = saveobj(self)
            if numel(self) > 1
                for nn = 1:numel(self)
                    s(nn,1) = saveobj(self(nn));
                end
            else
                s = self.struct;
                %
                % Remove unnecessary data
                %
                s.x = [];
                s.y = [];
                s.ODraw = [];
                s.ODcorr = [];
                s.raw.images = uint16(s.raw.images);
                for nn = 1:numel(s.clouds)
                    s.clouds(nn).fitdata.image = [];
                    s.clouds(nn).fitdata.bg = [];
                    s.clouds(nn).fitdata.residuals = [];
                end
            end
        end

    end

    methods(Static)
        function self = loadobj(a)
            %LOADOBJ Converts simpler object into ABSORPTIONIMAGE
            %
            %   C = LOADOBJ(A) converts simpler object A into instance C
            if numel(a) > 1
                for nn = 1:numel(a)
                    self(nn) = AbsorptionImage.loadobj(a(nn));
                end
            else
                self = AbsorptionImage;
                raw = RawImageData.loadobj(a.raw);
                c = AtomImageConstants.loadobj(a.constants);
                self.raw.copy(raw);
                self.constants.copy(c);

                self.makeImage;
                self.clouds = AtomCloud.empty;
                for nn = 1:numel(a.clouds)
                    self.clouds(nn,1) = AtomCloud.loadobj(a.clouds(nn));
                    self.clouds(nn,1).setConstants(self.constants);
                    self.clouds(nn,1).fitdata.makeFitObjects(self.x,self.y,self.ODcorr);
                end
            end
        end
        
        function [LabelStr,NumberStrTotal] = formatLabel(LabelCell,FormatCell,NumberCell)
            %FORMATLABEL Uses input labels, formats, and values to create a
            %nice label that can be printed to the command line
            %
            %   [LABELSTR,NUMBERSTR] =
            %   ABSORPTIONIMAGE.FORMATLABEL(LABELCELL,FORMATCELL,NUMBERCELL)
            %   creates a label string LABELSTR and number string NUMBERSTR
            %   using cell arrays LABELCELL, FORMATCELL, and NUMBERCELL
            N=numel(LabelCell);

            LabelStr='';
            for nn=1:N
                LabelMidPoint(nn)=numel(LabelStr)+round(numel(LabelCell{nn})/2);    %#ok
                if nn==N
                    LabelStr=[LabelStr,LabelCell{nn}];  %#ok
                else
                    LabelStr=[LabelStr,LabelCell{nn},'  |  '];  %#ok
                end

            end

            NumberStrTotal=repmat(' ',1,numel(LabelStr));
            for nn=1:N
                NumberStr=num2str(NumberCell{nn},FormatCell{nn});
                NumberLength=numel(NumberStr);
                NumberMidPoint=floor(NumberLength/2);
                NumberStrTotal((LabelMidPoint(nn)-NumberMidPoint+1):(LabelMidPoint(nn)+(NumberLength-NumberMidPoint)))=NumberStr;
            end
            NumberStrTotal(LabelStr=='|')='|';
        end
        
        function imgOut = match(img1,img2,exRegion)
            %MATCH Attempts to rescale a background image to more closely
            %match an image with atoms while excluding the region where the
            %atoms are
            %
            %   IMGOUT = MATCH(IMG1,IMG2,EXREGION) re-scales IMG2 to match
            %   IMG1 while excluding region EXREGION
            %
            Filt = ones(size(img1));
            for nn = 1:size(exRegion,1)
                row = exRegion{nn,1};
                col = exRegion{nn,2};
                Filt(row,col) = 0;
            end
            s = fminbnd(@(x) sum(sum(((img1-x*img2).*Filt)).^2),0.1,10);
            imgOut = img2*s;
        end
        
    end

end