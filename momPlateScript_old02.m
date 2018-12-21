%behold the unorganized mind of sarah!
%universal constants
kk=2*pi;
netta = 377;
gamma = 1.781072418; %cme-petterson-pg39, (2.13)

phiInc = pi/4;
thetaInc = pi/2;
%plate info
NumCells = 5;
NumEdges = NumCells-1;
lenx = 2.23;%lambda
leny= 2.23;%lambda

delx = lenx/NumCells; %a in the book
dely = leny/NumCells; %b in the book
%-----------------------Generate points---------------------------------


%--------------------------generate matrix equation-----------------------
%need to first generate x, y points
%going to have two sets depending on Ex or Ey
ex = zeros(NumCells,1);
ey = zeros(NumCells,1);

TxmBxn = zeros(NumCells);
TxmBxnsp = zeros(NumCells);

TxmByn = zeros(NumCells);

TymByn = zeros(NumCells);
TymBynsp = zeros(NumCells);

for mm = 1:NumCells %marches with x
    
    ex(mm) = 2*sin(phiInc)*sin(thetaInc)...
                *exp(-1j*kk*sin(thetaInc)*mm*(delx*cos(phiInc)+dely*sin(phiInc)))...
                *sin(kk*sin(thetaInc)*cos(phiInc)*delx/2);
    ey(mm) = 2*cos(phiInc)*sin(thetaInc)...
                *exp(-1j*kk*sin(thetaInc)*mm*(delx*cos(phiInc)+dely*sin(phiInc)))...
                *sin(kk*sin(thetaInc)*sin(phiInc)*delx/2);
    for nn = 1:NumCells %marches with y
        delmn = mm-nn;


            %define function
            %green's function shifted to xm-xn and ym-yn
            ggfun = @(xp,yp) exp(-1j*kk*sqrt( (xp-delx*delmn).^2 + (yp-dely*delmn).^2))...
                                ./ (4*pi*sqrt( (xp-delx*delmn).^2 + (yp-dely*delmn).^2));
            %spleen
            ggpp1xx = @(xp,yp) (9/8+xp.*3/(2*delx) +xp.^2./(2*delx^2)) .* ...
                                exp(-1j*kk*sqrt( (xp-delx*delmn).^2 + (yp-dely*delmn).^2))...
                                ./ (4*pi*sqrt( (xp-delx*delmn).^2 + (yp-dely*delmn).^2));          
                            
                            
            ggpp2xx = @(xp,yp) (3/4-xp.^2./ delx^2).* ...
                                exp(-1j*kk*sqrt( (xp-delx*delmn).^2 + (yp-dely*delmn).^2))...
                                ./ (4*pi*sqrt( (xp-delx*delmn).^2 + (yp-dely*delmn).^2));
                            
            ggpp3xx = @(xp,yp) (9/8-3*xp./(2*delx) +xp.^2/(2*delx.^2)).* ...
                                exp(-1j*kk*sqrt( (xp-delx*delmn).^2 + (yp-dely*delmn).^2))...
                                ./ (4*pi*sqrt( (xp-delx*delmn).^2 + (yp-dely*delmn).^2));
                            
            
            ggpp1yy = @(xp,yp) (9/8+yp.*3/(2*dely) +yp.^2./(2*dely^2)) .* ...
                                exp(-1j*kk*sqrt( (xp-delx*delmn).^2 + (yp-dely*delmn).^2))...
                                ./ (4*pi*sqrt( (xp-delx*delmn).^2 + (yp-dely*delmn).^2));          
                            
                            
            ggpp2yy = @(xp,yp) (3/4-yp.^2./ dely^2).* ...
                                exp(-1j*kk*sqrt( (xp-delx*delmn).^2 + (yp-dely*delmn).^2))...
                                ./ (4*pi*sqrt( (xp-delx*delmn).^2 + (yp-dely*delmn).^2));
                            
            ggpp3yy = @(xp,yp) (9/8-3*yp./(2*dely) +yp.^2/(2*dely.^2)).* ...
                                exp(-1j*kk*sqrt( (xp-delx*delmn).^2 + (yp-dely*delmn).^2))...
                                ./ (4*pi*sqrt( (xp-delx*delmn).^2 + (yp-dely*delmn).^2));
            
            %do integrals
            %first convolution
            TxmBxn(mm,nn) = integral2(ggfun,-delx*3/2,-delx/2,-dely/2,dely/2)...
                             +2*integral2(ggfun,-delx/2,delx/2,-dely/2,dely/2)...
                             +integral2(ggfun,delx/2,delx*3/2,-dely/2,dely/2);
            
            TxmBxnsp(mm,nn) = integral2(ggpp1xx,-delx*3/2,-delx/2,-dely/2,dely/2)...
                             +integral2(ggpp2xx,-delx/2,delx/2,-dely/2,dely/2)...
                             +integral2(ggpp3xx,delx/2,delx*3/2,-dely/2,dely/2);
            
            %second and third convolution
            %TxmByn== TynBxn expect for the 1/a vs 1/b... so only get the
            %matrix once... 
            TxmByn(mm,nn) = integral2(ggfun,-delx,0,-dely,0)-integral2(ggfun,-delx,0,0,dely)...
                           -integral2(ggfun,0,delx,-dely,0)+integral2(ggfun,0,delx,0,dely);
            %fourth convolution...if you're going to do a square you can
            %throw this out and just use the first convolution
            TymByn(mm,nn) = integral2(ggfun,-delx/2,delx/2,-dely*3/2,-dely/2)...
                         +2*integral2(ggfun,-delx/2,delx/2,-dely/2,dely/2)...
                         +  integral2(ggfun,-delx/2,delx/2,delx/2,delx*3/2);
                     
            TymBynsp(mm,nn) = integral2(ggpp1yy,-delx/2,delx/2,-dely*3/2,dely/2)...
                            +integral2(ggpp2yy,-delx/2,delx/2,-dely/2,dely/2)...
                            +integral2(ggpp2yy,-delx/2,delx/2,dely/2,dely*3/2);
                            
            

                                   
    end
end

AA = TxmBxn/delx+kk^2*delx*TxmBxnsp;
BB = TxmByn/delx;
CC = TxmByn/dely;
DD = TymByn/dely+kk^2*dely*TymBynsp;

EE = [ex;ey]
ZZ = (-netta/(1j*kk))*[AA BB; CC DD]

% %lets get quadratic spline
% qqreg1 = linspace(-3*delx/2,-delx/2,100);
% qqreg2 = linspace(-delx/2,delx/2,100);
% qqreg3 = linspace(delx/2,3*delx/2,100);
% 
% qq = [qqreg1, qqreg2, qqreg3];
% 
% f1 = 9/8 + 3*qqreg1/(2*delx) +qqreg1.^2/(2*delx.^2);
% f2 = 3/4-qqreg2.^2/delx.^2;
% f3 = 9/8 - 3*qqreg3/(2*delx) +qqreg3.^2/(2*delx.^2);
% 
% ff = [f1,f2,f3];

% figure;plot(qqreg1,f1)
% figure;plot(qqreg2,f2)
% figure;plot(qqreg3,f3)