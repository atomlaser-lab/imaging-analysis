classdef AtomCloudFit < handle
    %ATOMCLOUDFIT Defines a class for fitting data from absorption images
    %to particular distributions and extracting useful parameters from
    %those fits
    properties
        %
        % Settings for constructing the fit data
        %
        roiRow      %2-element vector specifying the start and end row values for the region of interest (ROI)
        roiCol      %2-element vector specifying the start and end column values for the ROI
        roiStep     %1 or 2 element vector specifying the step size for the ROI to reduce the number of points to fit
        imgSize     %The size of the image in question
        %
        % Properties describing the restricted data
        %
        image       %The actual image data, restricted to the ROI
        x           %The x coordinate for the fit
        xdata       %The marginal X distribution
        xfit        %The fit to xdata
        y           %The y coordinate for the fit
        ydata       %The marginal Y distribution
        yfit        %the fit to ydata
        %
        % Fitting parameters
        fittype     %The type of fit to use
        ex          %ODs above ex are excluded
        lb          %Lower bound for fit
        ub          %Upper bound for fit
        guess       %Guess for fit
        %
        % Fit results
        %
        residuals   %The fit residuals
        bg          %The fit background
        params      %The parameters as a CloudParameters object
    end

    methods
        function self = AtomCloudFit(varargin)
            %ATOMCLOUDFIT Creates an ATOMCLOUDFIT object
            %
            %   C = ATOMCLOUDFIT(NAME,VALUE,...) creates an ATOMCLOUDFIT
            %   object using parameters specified as NAME/VALUE pairs.  See
            %   SETUP() for more information
            self.imgSize = RawImageData.DEFAULT_SIZE;
            self.roiRow = [1,RawImageData.DEFAULT_SIZE(1)];
            self.roiCol = [1,RawImageData.DEFAULT_SIZE(2)];
            self.roiStep = 1;
            self.fittype = 'none';
            self.lb = CloudParameters([]);
            self.ub = CloudParameters([]);
            self.guess = CloudParameters([]);
            if nargin > 0
                self.set(varargin{:});
            end
        end

        function self = copy(self,obj)
            %COPY Copies object properties from one object to this one
            %
            %   C = C.COPY(OBJ) copies properties from OBJ to the current
            %   object C
            p = properties(self);
            for nn = 1:numel(p)
                self.(p{nn}) = obj.(p{nn});
            end
        end

        function self = set(self,varargin)
            %SET Sets parameters for the object
            %
            %   C = C.SET(NAME,VALUE,...) sets parameters given by NAME
            %   to VALUE for each pair of NAME and VALUE.  Valid options
            %   are 'roirow', 'roicol', 'roistep', and 'fittype'
            if mod(numel(varargin),2) ~= 0
                error('Arguments must be as name/value pairs');
            else
                for nn = 1:2:numel(varargin)
                    p = lower(varargin{nn});
                    v = varargin{nn+1};
                    switch p
                        case 'roirow'
                            self.roiRow = v;
                        case 'roicol'
                            self.roiCol = v;
                        case 'roistep'
                            self.roiStep = v;
                        case 'fittype'
                            self.fittype = v;
                        case 'imgsize'
                            self.imgSize = v;
                        case 'lb'
                            self.lb = v;
                        case 'ub'
                            self.ub = v;
                        case 'guess'
                            self.guess = v;
                        otherwise
                            error('Unsupported parameter %s',p);
                    end
                end
            end
        end

        function set.fittype(self,v)
            %SET.FITTYPE Sets the fit type
            %
            %   C.FITTYPE = V sets the fit type to V
            %
            v = lower(v);
            switch v
                case 'none'
                    self.fittype = 'none';
                case 'sum'
                    self.fittype = 'sum';
                case 'gauss1d'
                    self.fittype = 'gauss1d';
                case {'twocomp1d','2comp1d'}
                    self.fittype = 'twocomp1d';
                case {'tf1d','bec1d'}
                    self.fittype = 'tf1d';
                case 'gauss2d'
                    self.fittype = 'gauss2d';
                case {'bec2d','tf2d'}
                    self.fittype = 'tf2d';
                case {'twocomp2d','2comp2d'}
                    self.fittype = 'twocomp2d';
                otherwise
                    error('Fit type ''%s'' not supported',v);
            end
        end

        function set.roiRow(self,v)
            %SET.ROIROW Sets the row for the ROI
            %
            %   C.ROIROW = V sets the ROI row to V
            %
            if numel(v) ~= 2
                error('Must supply a 2-element vector');
            else
                v(1) = max(v(1),1);
                v(2) = min(v(2),self.imgSize(1));
                self.roiRow = v;
            end
        end

        function set.roiCol(self,v)
            %SET.ROICOL Sets the column for the ROI
            %
            %   C.ROICOL = V sets the ROI column to V
            %
            if numel(v) ~= 2
                error('Must supply a 2-element vector');
            else
                v(1) = max(v(1),1);
                v(2) = min(v(2),self.imgSize(2));
                self.roiCol = v;
            end
        end

        function set.roiStep(self,v)
            %SET.ROISTEP Sets the step size for the ROI
            %
            %   C.ROISTEP = V sets the ROI step to V. V can be either a 1
            %   or 2 element vector, where the elements of the latter
            %   correspond to the step size for the row and column
            %   directions, respectively
            %
            if numel(v) == 1
                self.roiStep = [v,v];
            elseif numel(v) == 2
                self.roiStep = v;
            else
                error('Must supply a 1 or 2 element vector');
            end
        end

        function self = checkROI(self)
            %CHECKROI Checks the ROI to make sure it is within the image bounds
            %and throws and error if it is not
            if self.roiRow(1) < 1
                error('ROI starts at a row (%d) less than 1',self.roiRow(1));
            elseif self.roiRow(2) > self.imgSize(1)
                error('ROI ends at a row (%d) larger than the image size of %d. Check the ROI or the value of ''imgSize''',self.roiRow(2),self.imgSize(1));
            end
            if self.roiCol(1) < 1
                error('ROI starts at a column (%d) less than 1',self.roiCol(1));
            elseif self.roiCol(2) > self.imgSize(2)
                error('ROI ends at a column (%d) larger than the image size of %d. Check the ROI or the value of ''imgSize''',self.roiCol(2),self.imgSize(2));
            end
        end

        function self = coerceROI(self)
            %COERCEROI Coerces ROI to be within the image bounds
            self.roiRow = [max(self.roiRow(1),1),min(self.roiRow(2),self.imgSize(1))];
            self.roiCol = [max(self.roiCol(1),1),min(self.roiCol(2),self.imgSize(2))];
        end

        function r = is1D(self)
            %IS1D Indicates if the fit type is a 1D fit type
            %
            %   R = C.IS1D() returns true if the fit type is a 1D type and
            %   false if it is not
            %
            v = {'gauss2d','gauss2dangle','twocomp2d','tf2d','none'};
            r = true;
            for nn = 1:numel(v)
                if strcmpi(v{nn},self.fittype)
                    r = false;
                    return
                end
            end
        end

        function [row,col] = makeROIVectors(self)
            %MAKEROIVECTORS Creates vectors corresponding to the ROI
            %
            %   [ROW,COL] = MAKEROIVECTORS() returns vectors ROW and COL
            %   corresponding to the indices of the ROI row and col.
            %
            if numel(self.roiStep) == 1
                row = self.roiRow(1):self.roiStep:self.roiRow(2);
                col = self.roiCol(1):self.roiStep:self.roiCol(2);
            elseif numel(self.roiStep) == 2
                row = self.roiRow(1):self.roiStep(1):self.roiRow(2);
                col = self.roiCol(1):self.roiStep(2):self.roiCol(2);
            else
                error('ROI step should be either 1 or 2 elements!');
            end
        end

        function self = makeFitObjects(self,x,y,image)
            %MAKEFITOBJECTS Creates internal input fit objects based on the
            %x and y positions and the image
            %
            %   C = C.MAKEFITOBJECTS(X,Y,IMAGE) Sets internal fit objects
            %   based on the input x and y positions and the image data.
            %   This function sets the properties x, xdata, y, ydata, and
            %   image
            %
            
            %
            %Coerce ROI to lie within image boudns
            %
            self.imgSize = size(image);
            self.coerceROI();
            %
            % Restrict input data to the ROI
            %
            [row,col] = self.makeROIVectors;
            self.image = image(row,col);
            if numel(self.image) == 0
                error('ROI contains no elements!');
            end
            self.x = x(col);
            self.xdata = sum(image(row,col),1);
            self.y = y(row);
            self.ydata = sum(image(row,col),2);
            %
            % Sanitise data so that there are no NaNs or Infs or complex
            % numbers
            %
            idx = ~(isnan(self.xdata) | (imag(self.xdata) ~= 0) | isinf(self.xdata));
            self.xdata = self.xdata(idx);
            self.x = self.x(idx);
            self.xfit = zeros(size(self.x));
            idx = ~(isnan(self.ydata) | (imag(self.ydata) ~= 0) | isinf(self.ydata));
            self.ydata = self.ydata(idx);
            self.y = self.y(idx);
            self.yfit = zeros(size(self.y));
        end

        function self = fit(self,fittype,ex)
            %FIT Fits atomic cloud and sets internal params property
            %
            %   C = C.FIT() fits the image data according to the internal
            %   fit type
            %
            %   C = C.FIT(FITTYPE) fits the image data according to FITTYPE
            %
            %   C = C.FIT(__,EX) fits the image data while excluding ODs
            %   above EX
            if nargin > 1 && ~isempty(fittype)
                self.fittype = fittype;
            end
            if nargin < 3
                ex = [];
            end
            %
            % Perform fit according to fit type
            %
            self.bg = [];
            self.residuals = [];
            self.params = CloudParameters([]);
            switch self.fittype
                case {'none','sum'}
                    self.params = CloudParameters();
                case {'gauss1d','twocomp1d','tf1d'}
                    [px,self.xfit,self.bg.x] = self.fit1D(self.fittype,self.x,self.xdata,self.lb,self.ub,self.guess,ex);
                    [py,self.yfit,self.bg.y] = self.fit1D(self.fittype,self.y,self.ydata,self.lb,self.ub,self.guess,ex);
                    self.residuals.x = self.xdata(:) - self.xfit(:);
                    self.residuals.y = self.ydata(:) - self.yfit(:);
                    self.params = CloudParameters(px,py);
                case {'gauss2d','gauss2dangle','twocomp2d','tf2d'}
                    [self.params,f,self.bg] = self.fit2D(self.fittype,self.x,self.y,self.image,self.lb,self.ub,self.guess,ex);
                    self.xfit = sum(f,1);self.yfit = sum(f,2);
                    self.residuals = self.image - f;
                otherwise
                    error('Fit type %s not supported',self.fittype);
            end
        end
        
        %% Saving and loading functions
        function s = struct(self)
            %STRUCT Creates a structure from the object instance
            %
            %   S = F.STRUCT() creates a structure S from the object
            %   instance F
            
            %
            % Settings for constructing fit data
            %
            s.roiRow = self.roiRow;
            s.roiCol = self.roiCol;
            s.roiStep = self.roiStep;
            s.imgSize = self.imgSize;
            %
            % Properties describing restricted data
            %
            s.image = self.image;
            s.x = self.x;
            s.xdata = self.xdata;
            s.xfit = self.xfit;
            s.y = self.y;
            s.ydata = self.ydata;
            s.yfit = self.yfit;
            %
            % Fitting parameters
            %
            s.fittype = self.fittype;
            s.ex = self.ex;
            s.lb = self.lb.struct;
            s.ub = self.ub.struct;
            s.guess = self.guess.struct;
            %
            % Fit results
            %
            s.residuals = self.residuals;
            s.bg = self.bg;
            s.params = self.params.struct;
        end
        
        function s = saveobj(self)
            %SAVEOBJ Saves an instance of the ATOMCLOUDFIT class as a
            %simpler structure
            %
            %   S = F.SAVEOBJ() saves instance F as simpler structure S
            s = self.struct;
        end

    end

    methods(Static)
        function b = loadobj(a)
            %LOADOBJ Converts saved structure into class instance
            %
            %   F = LOADOBJ(A) converts simpler structure A into proper
            %   instance of ATOMCLOUDFIT F
            p = fields(a);
            b = AtomCloudFit;
            for nn = 1:numel(p)
                b.(p{nn}) = a.(p{nn});
            end
            b.params = CloudParameters.loadobj(b.params);
            b.lb = CloudParameters.loadobj(b.lb);
            b.ub = CloudParameters.loadobj(b.ub);
            b.guess = CloudParameters.loadobj(b.guess);
        end
        
        function opt = getoptions
            %GETOPTIONS creates options object with default parameters
            %
            %   OPT = GETOPTIONS returns OPT as the options object
            opt = optimset('Display','off', 'MaxFunEvals',1000, 'TolFun', 1e-7, 'TolX', 1e-6);
        end
        
        function y = bg1D(c,pos,p0)
            %BG1D defines a 1D background function comprising an offset and
            %a linearly varying term
            %
            %   F = BG1D(C,POS) returns a 1D background term where C(1) is
            %   a constant offset and C(2) is the first derivative. POS is
            %   the coordinate variable
            %
            %   F = BG1D(__,P0) shifts the linear term to be zero at P0
            if nargin < 3
                p0 = 0;
            end
            y = c(1) + c(2)*(pos - p0);
        end
        
        function z = bg2D(c,pos,p0)
            %BG2D defines a 2D background function comprising an offset and
            %a linearly varying term in each of the x and y directions
            %
            %   F = BG2D(C,POS) returns a 2D background term where C(1) is
            %   a constant offset, C(2) is the first derivative along X,
            %   and C(3) is the first derivative along Y. POS is an NxNx2
            %   array where POS(:,:,1) is the gridded X values and
            %   POS(:,:,2) is the gridded Y values
            %
            %   F = BG1D(__,P0) shifts the linear terms to be zero at P0
            X = pos(:,:,1);
            Y = pos(:,:,2);
            if nargin < 3
                p0 = [0,0];
            end
            z = c(1) + c(2)*(X - p0(1)) + c(3)*(Y - p0(2));
        end
        
        function y = gauss1D(c,pos)
            %GAUSS1D defines a 1D gaussian function with offset and linear
            %offset term
            %
            %   Y = GAUSS1D(C,X) returns a 1D Gaussian shape where C(1) is
            %   the amplitude, C(2) is the offset, C(3) is the standard
            %   deviation, C(4) is a constant offset, and C(5) is a
            %   linearly varying offset term (see BG1D())
            %
            A = c(1);
            x0 = c(2);
            xw = c(3);
            y = A*exp(-(pos-x0).^2./(2*xw^2)) + AtomCloudFit.bg1D(c(end-1:end),pos,x0);
        end

        function y = bec1D(c,pos)
            %BEC1D defines a 1D Thomas-Fermi function with offset and
            %linear offset term
            %
            %   Y = BEC1D(C,X) returns a 1D Thomas-Fermi distribution with
            %   amplitude C(1), offset C(2), TF radius C(3), constant
            %   offset C(4), and linearly varying offset C(5)
            %
            A = c(1);
            x0 = c(2);
            xw = c(3);
            s2 = ((pos-x0)/xw).^2;
            y = A*((1 - s2).*(s2 <= 1)).^2 + AtomCloudFit.bg1D(c(end-1:end),pos,x0);
        end

        function y = twoComp1D(c,pos)
            %TWOCOMP1D defines a 1D two-component fit function with offset
            %and linear offset term
            %
            %   Y = TWOCOMP1D(C,POS) returns a two-component distribution
            %   with Gaussian amplitude C(1), position offset C(2),
            %   Gaussian width C(3), TF amplitude C(4), TF radius C(5),
            %   constant offset C(6), and linearly varying offset C(7)
            %
            Ag = c(1);
            x0 = c(2);
            xwg = c(3);
            Ab = c(4);
            xwb = c(5);
            s2 = ((pos-x0)/xwb).^2;
            y = Ag*exp(-(pos-x0).^2./(2*xwg^2)) + Ab.*((1-s2).*(s2<=1)).^2 ...
                + AtomCloudFit.bg1D(c(end-1:end),pos,x0);
        end

        function F = gauss2D(c,Z)
            %GAUSS2D defines a 2D Gaussian distribution
            %
            %   F = GAUSS2D(C,Z) where Z is a 3D array where Z(:,:,1)
            %   the X array and Z(:,:,2) the Y array.  C is s.t. C(1) is
            %   the amplitude, C(2) is the x offset, C(3) is the x standard
            %   deviation, C(4) is the y offset, C(5) is the y standard
            %   deviation, C(6) is the linear term in X, C(7) is the linear
            %   term in Y, c(8) is the offset
            %
            x = Z(:,:,1);
            y = Z(:,:,2);
            
            A = c(1);
            x0 = c(2);
            xw = c(3);
            y0 = c(4);
            yw = c(5);
            
            F = A*exp(-(x-x0).^2./(2*xw.^2)-(y-y0).^2./(2*yw.^2))...
                + AtomCloudFit.bg2D(c(end-2:end),Z,[x0,y0]);
        end
        
        function F = gauss2DAngle(c,Z)
            %GAUSS2DANGLE defines a 2D Gaussian distribution that is
            %rotated by an angle
            %
            %   F = GAUSS2DANGLE(C,Z) where Z is a 3D array where Z(:,:,1)
            %   the X array and Z(:,:,2) the Y array.  C is s.t. C(1) is
            %   the amplitude, C(2) is the x offset, C(3) is the x standard
            %   deviation, C(4) is the y offset, C(5) is the y standard
            %   deviation, C(6) is the linear term in X, C(7) is the linear
            %   term in Y, c(8) is the offset, and C(9) is the rotation
            %   angle
            %
            x = Z(:,:,1);
            y = Z(:,:,2);
            
            A = c(1);
            x0 = c(2);
            xw = c(3);
            y0 = c(4);
            yw = c(5);
            th = c(6);
            
            xr = x - x0;
            yr = y - y0;
            
            F = A*exp(-(xr*cos(th)+yr*sin(th)).^2./(2*xw.^2)-(-xr*sin(th)+yr*cos(th)).^2./(2*yw.^2))...
                + AtomCloudFit.bg2D(c(end-2:end),Z,[x0,y0]);
        end
        
        function F = bec2D(c,Z)
            %BEC2D defines a 2D Thomas-Fermi distribution
            %
            %   F = BEC2D(C,Z) where Z is a 3D array where Z(:,:,1) the X
            %   array and Z(:,:,2) the Y array.  C is s.t. C(1) is the
            %   amplitude, C(2) is the x offset, C(3) is the x width, C(4)
            %   is the y offset, C(5) is the y width, C(6) is the linear
            %   term in X, C(7) is the linear term in Y, c(8) is the offset
            %
            n0 = c(1);
            x0 = c(2);
            xw = c(3);
            y0 = c(4);
            yw = c(5);
%             linx = c(6);
%             liny = c(7);
%             offset = c(8);
            
            x = (Z(:,:,1)-x0)/xw;
            y = (Z(:,:,2)-y0)/yw;
            s2 = x.^2+y.^2;
            F = n0*((1-s2).*(s2<=1)).^1.5...
                + AtomCloudFit.bg2D(c(end-2:end),Z,[x0,y0]);
        end
        
        function F = twoComp2D(c,Z)
            %TWOCOMP2D defines a 2D two-component distribution
            %
            %   F = TWOCOMP2D(C,Z) where Z is a 3D array where Z(:,:,1) the X
            %   array and Z(:,:,2) the Y array.  C is s.t. C(1) is the
            %   BEC amplitude, C(2) is the x offset, C(3) is the x BEC width, C(4)
            %   is the y offset, C(5) is the y BEC width, C(6) is the linear
            %   term in X, C(7) is the linear term in Y, c(8) is the
            %   offset, C(9) is the Gaussian amplitude, C(10) is the x
            %   Gaussian standard deviation, and C(11) is the y Gaussian
            %   standard deviation
            %
            Ag = c(1);
            x0 = c(2);
            xwg = c(3);
            y0 = c(4);
            ywg = c(5);
            Ab = c(6);
            xwb = c(7);
            ywb = c(8);
            x = (Z(:,:,1)-x0);
            y = (Z(:,:,2)-y0);
            s2 = (x/xwb).^2 + (y/ywb).^2;
            F = Ag*exp(-x.^2./(2*xwg.^2)-y.^2./(2*ywg.^2))...
                + Ab*((1 - s2).*(s2 <= 1)).^1.5...
                + AtomCloudFit.bg2D(c(end-2:end),Z,[x0,y0]);
            
        end

        function guess = guessGauss1DParams(x,y)
            %GUESSGAUSS1DPARAMS Guesses 1D Gaussian parameters
            %
            %   GUESS = GUESSGAUSS1DPARAMS(X,Y) Guesses 1D Gaussian
            %   parameters using position data X and OD data Y
            %
            
            % Estimate of fitting params...
            Az = min(y);
            [Max,I] = max(smooth(y,7));
            Bz = Max - Az;
            Cz = x(I);
            tmp=y-Az;
            tmp=tmp/max(tmp);
            width_point = exp(-0.5);

            I1 = find(smooth(tmp,3) >= width_point, 2, 'first');
            I2 = find(smooth(tmp,3) >= width_point, 2, 'last');
            Dz1 = mean(x(I1));
            Dz2 = mean(x(I2));
            Dztrial = Dz2 - Dz1;

            if length(I1) < 1 || Dztrial < 10e-6
                Dz = 10e-6;
            else
                Dz = Dztrial;
            end

            guess = CloudParameters([]);
            guess = guess.set('offset',Az,'pos',Cz,'gaussamp',Bz,...
                'gausswidth',Dz,'lin',0);
        end
        
        function g = guessGauss2DParams(x,y,z)
            %GUESSGAUSS2DPARAMS Guesses 2D Gaussian parameters
            %
            %   GUESS = GUESSGAUSS2DPARAMS(X,Y,Z) Guesses 1D Gaussian
            %   parameters using position data X and Y, and OD data Z
            %
            gx = AtomCloudFit.guessGauss1DParams(x,sum(z,1));
            gy = AtomCloudFit.guessGauss1DParams(y,sum(z,2));
            g = CloudParameters([]);
            g.gaussAmp = max(z(:));
            g.offset = min(z(:));
            g.pos = [gx.pos,gy.pos];
            g.gaussWidth = [gx.gaussWidth,gy.gaussWidth];
            g.lin = [0,0];
        end
        
        function g = guessTF2DParams(x,y,z)
            %GUESSTF2DPARAMS Guesses 2D Thomas-Fermi parameters
            %
            %   GUESS = GUESSTF2DPARAMS(X,Y,Z) Guesses 2D TF
            %   parameters using position data X and Y, and OD data Z
            %
            gx = AtomCloudFit.guessTF1DParams(x,sum(z,1));
            gy = AtomCloudFit.guessTF1DParams(y,sum(z,2));
            g = CloudParameters([]);
            g.becAmp = max(z(:));
            g.offset = min(z(:));
            g.pos = [gx.pos,gy.pos];
            g.becWidth = [gx.becWidth,gy.becWidth];
            g.lin = [0,0];
        end

        function guess = guessTF1DParams(x,y)
            %GUESTF1DPARAMS Guesses 1D TF parameters
            %
            %   GUESS = GUESSTF1DPARAMS Guesses 1D TF parameters
            %   using position data X and OD data Y
            %
            Az = min(y);
            [B,I] = max(y);
            Bz = 2/pi*(B-Az);
            Cz = x(I);

            tmp = y - Az;
            tmp = tmp./max(tmp);
            wp = 0.5;

            I1 = find(smooth(tmp,3) >= wp, 2, 'first');
            I2 = find(smooth(tmp,3) >= wp, 2, 'last');
            Dz1 = mean(x(I1));
            Dz2 = mean(x(I2));
            Dztrial = Dz2 - Dz1;

            if length(I1) < 1 || Dztrial < 5e-6
                Dz = 5e-6;
            else
                Dz = Dztrial;
            end

            guess = CloudParameters([]);
            guess = guess.set('offset',Az,'pos',Cz,'becamp',Bz,...
                'becwidth',Dz,'lin',0);
        end

        function [lb,ub] = guessBounds1D(x,y)
            %GUESSBOUNDS1D Guesses parameter bounds for 1D fits
            %
            %   [LB,UB] = GUESSBOUNDS1D(X,Y) Guesses upper and lower bounds
            %   UB and LB given coordinate data X and ordinate data Y.  UB
            %   and LB are returned as CLOUDPARAMETER objects
            lb = CloudParameters('offset',-100,'gaussamp',0,'becamp',0,...
                'pos',min(x),'gausswidth',5*abs(diff(x(1:2))),...
                'becwidth',5*abs(diff(x(1:2))),'lin',-10/range(x));
            ub = CloudParameters('offset',100,'gaussamp',1e3,'becamp',1e3,...
                'pos',max(x),'gausswidth',range(x),...
                'becwidth',range(x),'lin',+10/range(x));
        end
        
        function [lb,ub] = guessBounds2D(x,y,z)
            %GUESSBOUNDS2D Guesses parameter bounds for 2D fits
            %
            %   [LB,UB] = GUESSBOUNDS2D(X,Y,Z) Guesses upper and lower
            %   bounds UB and LB given coordinate data X and Y and ordinate
            %   data Y.  UB and LB are returned as CLOUDPARAMETER objects
            lb = CloudParameters('offset',-0.1,'gaussamp',0,'becamp',0,...
                'pos',[min(x),min(y)],'gausswidth',5*abs([diff(x(1:2)),diff(y(1:2))]),...
                'becwidth',5*abs([diff(x(1:2)),diff(y(1:2))]),'lin',-0.1./[range(x),range(y)],...
                'cloudangle',-pi/6);
            ub = CloudParameters('offset',0.1,'gaussamp',10,'becamp',10,...
                'pos',[max(x),max(y)],'gausswidth',[range(x),range(y)],...
                'becwidth',[range(x),range(y)],'lin',+10./[range(x),range(y)],...
                'cloudangle',pi/6);
        end
        
        function params = attemptFit(f,guess,x,y,lb,ub,options)
            %ATTEMPTFIT Attempts to fit data to a given fit function
            %
            %   PARAMS = ATTEMPTFIT(FUNC,GUESS,X,Y,LB,UB,OPTIONS) attempts
            %   to fit position data X and OD data Y to function FUNC using
            %   an initial guess GUESS and upper and lower bounds UB and
            %   LB with options OPTIONS.
            try
                params = lsqcurvefit(f, guess, x, y, lb, ub, options);
            catch err
                if (strcmp(err.identifier,'optimlib:snls:UsrObjUndefAtX0'))
                    disp('No atoms found in ROI for fit');
                    params = zeros(size(guess));
                else
                    rethrow(err);
                end
            end
        end

        function [p,f,bg] = fit1D(fittype,x,y,lbi,ubi,gi,ex)
            %FIT1D Fits a 1D distribution to 1D data
            %
            %   [P,F,BG] = FIT1D(FITTYPE,X,Y) fits 1D data described by
            %   coordainte X and ordinate Y using the fit type given by
            %   FITTYPE. Returns parameters as CLOUDPARAMETERS object P and
            %   fit data F. Background field is given by BG.  FITTYPE is
            %   any 1D fit allowed by the ATOMCLOUDFIT classes 'fittype'
            %   property
            %
            %   [P,F,BG] = FIT1D(__,LBI) uses user-supplied lower bounds
            %   LBI as a CLOUDPARAMETERS object.  Set values to non-empty
            %   values to use as lower bound.  Set LBI to [] to disable
            %
            %   [P,F,BG] = FIT1D(__,UBI) uses user-supplied upper bounds
            %   UBI. Same constraints as for LBI
            %
            %   [P,F,BG] = FIT1D(__,GI) uses user-supplied guesses GI.
            %   Same constraints as for LBI and UBI
            %
            %   [P,F.BG] = FIT1D(__,EX) excludes ordinate data with values
            %   larger than EX. Set to [] to disable
            %
            x = x(:);y = y(:);  %Reshape data
            %
            % Exclude data if ex is set and not-empty
            %
            if nargin >= 7 && ~isempty(ex)
                [xex,yex] = AtomCloudFit.exclusion(ex,x,y);
            else
                xex = x;yex = y;
            end
            %
            % Set fit function and guesses based on fit type
            %
            switch lower(fittype)
                case {'gauss1d','gauss'}
                    g = AtomCloudFit.guessGauss1DParams(xex,yex);
                    func = @(c,x) AtomCloudFit.gauss1D(c,x);
                    order = {'gaussamp','pos','gausswidth','offset','lin'};
                case {'bec1d','tf1d','bec','tf'}
                    g = AtomCloudFit.guessTF1DParams(xex,yex);
                    func = @(c,x) AtomCloudFit.bec1D(c,x);
                    order = {'becamp','pos','becwidth','lin','offset'};
                case {'2comp1d','2comp','twocomp1d','twocomp'}
                    g = AtomCloudFit.guessGauss1DParams(xex,yex);
                    gtf = AtomCloudFit.guessTF1DParams(xex,yex);
                    g.compare(gtf);
                    g.gaussAmp = g.gaussAmp/2;
                    g.becAmp = g.becAmp/2;
                    func = @(c,x) AtomCloudFit.twoComp1D(c,x);
                    order = {'gaussamp','pos','gausswidth','becamp','becwidth','offset','lin'};
            end
            %
            % Get default bounds and guesses
            %
            options = AtomCloudFit.getoptions;
            [lb,ub] = AtomCloudFit.guessBounds1D(x,y);
            %
            % Overwrite those bounds/guesses with user-supplied
            % bounds/gueses
            %
            if nargin > 3 && ~isempty(lbi)
                lb = lb.compare(lbi);
            end
            if nargin > 4 && ~isempty(ubi)
                ub = ub.compare(ubi);
            end
            if nargin > 5 && ~isempty(gi)
                g = g.compare(gi);
            end
            %
            % Convert to array values
            %
            lb = lb.convert2array(order{:});
            ub = ub.convert2array(order{:});
            guess = g.convert2array(order{:});
            %
            % Perform fit
            %
            params = AtomCloudFit.attemptFit(func,guess,xex,yex,lb,ub,options);
            %
            % Extract parameters
            %
            p = CloudParameters.convertFromArray(order,params);
            f = func(params,x);
            bg = AtomCloudFit.bg1D(params(end-1:end),x,p.pos);
        end
        
        function [p,f,bg] = fit2D(fittype,x,y,z,lbi,ubi,gi,ex)
            %FIT2D Fits a 2D distribution to 1D data
            %
            %   [P,F,BG] = FIT2D(FITTYPE,X,Y) fits 2D data described by
            %   coordaintes X and Y and ordinate Z using the fit type given by
            %   FITTYPE. Returns parameters as CLOUDPARAMETERS object P and
            %   fit data F. Background field is given by BG.  FITTYPE is
            %   any 2D fit allowed by the ATOMCLOUDFIT classes 'fittype'
            %   property
            %
            %   [P,F,BG] = FIT1D(__,LBI) uses user-supplied lower bounds
            %   LBI as a CLOUDPARAMETERS object.  Set values to non-empty
            %   values to use as lower bound.  Set LBI to [] to disable
            %
            %   [P,F,BG] = FIT1D(__,UBI) uses user-supplied upper bounds
            %   UBI. Same constraints as for LBI
            %
            %   [P,F,BG] = FIT1D(__,GI) uses user-supplied guesses GI.
            %   Same constraints as for LBI and UBI
            %
            %   [P,F.BG] = FIT1D(__,EX) excludes ordinate data with values
            %   larger than EX. Set to [] to disable
            %

            %
            % Exclude data if ex is set and not-empty
            %
            [X,Y] = meshgrid(x,y);
            if nargin >= 8 && ~isempty(ex)
                [Xex,Yex,zex] = AtomCloudFit.exclusion(ex,X,Y,z);
            else
                Xex = X;Yex = Y;zex = z;
            end
            Position(:,:,1) = Xex;
            Position(:,:,2) = Yex;
            %
            % Set fit function and guesses based on fit type
            %
            switch lower(fittype)
                case {'gauss2d','gauss'}
                    g = AtomCloudFit.guessGauss2DParams(x,y,zex);
                    func = @(c,x) AtomCloudFit.gauss2D(c,x);
                    order = {'gaussamp','posx','gausswidthx',...
                        'posy','gausswidthy','offset','linx','liny'};
                case {'gauss2dangle','gaussangle'}
                    g = AtomCloudFit.guessGauss2DParams(x,y,zex);
                    func = @(c,x) AtomCloudFit.gauss2D(c,x);
                    order = {'gaussamp','posx','gausswidthx',...
                        'posy','gausswidthy','cloudangle',...
                        'offset','linx','liny'};
                case {'bec2d','tf2d','bec','tf'}
                    g = AtomCloudFit.guessTF2DParams(x,y,zex);
                    func = @(c,x) AtomCloudFit.bec2D(c,x);
                    order = {'becamp','posx','becwidthx',...
                        'posy','becwidthy','offset','linx','liny'};
                case {'2comp2d','2comp','twocomp2d','twocomp'}
                    g = AtomCloudFit.guessGauss2DParams(x,y,zex);
                    gtf = AtomCloudFit.guessTF2DParams(x,y,zex);
                    
                    g.compare(gtf);
                    g.gaussAmp = g.gaussAmp/2;
                    g.becAmp = g.becAmp/2;
                    func = @(c,x) AtomCloudFit.twoComp2D(c,x);
                    order = {'gaussamp','posx','gausswidthx',...
                        'posy','gausswidthy','becamp','becwidthx',...
                        'becwidthy','offset','linx','liny'};
            end
            %
            % Get default bounds and guesses
            %
            options = AtomCloudFit.getoptions;
            [lb,ub] = AtomCloudFit.guessBounds2D(x,y,z);
            %
            % Overwrite those bounds/guesses with user-supplied
            % bounds/gueses
            %
            if nargin > 4 && ~isempty(lbi)
                lb.compare(lbi);
            end
            if nargin > 5 && ~isempty(ubi)
                ub.compare(ubi);
            end
            if nargin > 6 && ~isempty(gi)
                g.compare(gi);
            end
            %
            % Convert to array values
            %
            lb = lb.convert2array(order{:});
            ub = ub.convert2array(order{:});
            guess = g.convert2array(order{:});
            %
            % Perform fit
            %
            params = AtomCloudFit.attemptFit(func,guess,Position,zex,lb,ub,options);
            %
            % Extract parameters
            %
            p = CloudParameters.convertFromArray(order,params);
            f = func(params,Position);
            bg = AtomCloudFit.bg2D(params(end-2:end),Position,p.pos);
        end

    end

    


end