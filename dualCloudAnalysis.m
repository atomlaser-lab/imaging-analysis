function [p,f,bg] = dualCloudAnalysis(fitdata,varargin)

    xx = fitdata.x;
    yy = fitdata.y;
    %
    % Defaults
    %
    x = mean(fitdata.x)*[1,1];
    y = mean(fitdata.y) + 0.25*range(fitdata.y).*[-1,1];
    w = [750e-6,750e-6];
    a = [0.1,0.1];
    %
    % Process variable arguments
    %
    if mod(numel(varargin),2) ~= 0
        error('Arguments must appear as name/value pairs!');
    else
        for nn = 1:2:numel(varargin)
            v = varargin{nn+1};
            switch lower(varargin{nn})
                case 'x1'
                    x(1) = v;
                case 'x2'
                    x(2) = v;
                case 'x'
                    x = v;
                case 'y1'
                    y(1) = v;
                case 'y2'
                    y(2) = v;
                case 'y'
                    y = v;
                case 'w1'
                    w(1) = v;
                case 'w2'
                    w(2) = v;
                case 'w'
                    w = v;
            end
        end 
    end

    options = AtomCloudFit.getoptions;
    minWidth = 600e-6;
    maxWidth = 850e-6;
    lb1 = [0,min(xx),minWidth,min(yy),minWidth];
    lb2 = [0,min(xx),minWidth,mean(yy),minWidth];
    ub1 = [2,max(xx),maxWidth,mean(yy),maxWidth];
    ub2 = [2,max(xx),maxWidth,max(yy),maxWidth];
    g1 = [a(1),x(1),w(1),y(1),w(1)];
    g2 = [a(2),x(2),w(2),y(2),w(2)];

    lb = [lb1,lb2,-0.05,-0.02/range(xx),-0.02/range(yy)];
    ub = [ub1,ub2,+0.05,+0.02/range(xx),+0.02/range(yy)];
    guess = [g1,g2,0,0,0];

    z = fitdata.image;
    [X,Y] = meshgrid(fitdata.x  - mean(fitdata.x)*0,fitdata.y - mean(fitdata.y)*0);
    pos(:,:,1) = X;
    pos(:,:,2) = Y;
    
    params = lsqcurvefit(@dualCloudFitFunc,guess,pos,z,lb,ub,options);
    p = CloudParameters('offset',params(11),'lin',params(12:13),...
        'pos',params([2,4]),'becamp',params(1),'becwidth',params([3,5]));
    p(2) = CloudParameters('offset',params(11),'lin',params(12:13),...
        'pos',params([7,9]),'becamp',params(6),'becwidth',params([8,10]));
    f = dualCloudFitFunc(params,pos);
    bg = p(1).offset + p(1).lin(1).*(X - p(1).pos(1)) + p(1).lin(2).*(Y - p(1).pos(2));

end

function z = dualCloudFitFunc(c,pos)

    X = pos(:,:,1);
    Y = pos(:,:,2);

    a(1) = c(1);
    x(1) = c(2);
    wx(1) = c(3);
    y(1) = c(4);
    wy(1) = c(5);
    s1 = ((X - x(1))./wx(1)).^2 + ((Y - y(1))./wy(1)).^2;

    a(2) = c(6);
    x(2) = c(7);
    wx(2) = c(8);
    y(2) = c(9);
    wy(2) = c(10);
    s2 = ((X - x(2))./wx(2)).^2 + ((Y - y(2))./wy(2)).^2;

    offset = c(11);
    linx = c(12);
    liny = c(13);

    z = a(1).*((1 - s1).*(s1 <= 1)).^1.5...
        + a(2).*((1 - s2).*(s2 <= 1)).^1.5...
        + offset + linx.*(X - x(1)) + liny.*(Y - y(1));

    % z = a(1)*exp(-(X-x(1)).^2./(2*wx(1)^2)-(Y-y(1)).^2./(2*wy(1)^2))...
    %     + a(2)*exp(-(X-x(2)).^2./(2*wx(2)^2)-(Y-y(2)).^2./(2*wy(2)^2))...
    %     + offset + linx.*X + liny.*Y;

end