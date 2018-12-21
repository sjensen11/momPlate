%behold the unorganized mind of sarah!
%universal constants
kk=2*pi;
netta = 377;
gamma = 1.781072418; %cme-petterson-pg39, (2.13)

NumCells = 5;
%plate info
lenx = 2.23;%lambda
leny= 2.23;%lambda
delx = lenx/NumCells;
dely = leny/NumCells;

%need to first generate x, y points
%going to have two sets depending on Ex or Ey
xxEx = zeros(1,NumCells);
yyEx = zeros(1,NumCells);
% 
% for ii=1:NumCells
%     for jj=1:NumCells
%         xxEx(1,(ii-1)*NumCells+jj) = delx*jj;
%         yyEx(1,(ii-1)*NumCells+jj) = dely*ii-dely/2;
%     end
% end
TxmByn1 = zeros(NumCells);
TxmByn2 = zeros(NumCells);
TxmByn3 = zeros(NumCells);
TxmByn4 = zeros(NumCells);

BB = zeros(NumCells);
for mm = 1:NumCells %marches with x
    for nn = 1:NumCells %marches with y
        mm
        delmn = mm-nn;


%             ggfun = @(xp,yp) exp(-1j*kk*sqrt( (xp).^2 + (yp).^2)) ./ (4*pi*sqrt( (xp).^2 + (yp).^2));
            ggfun = @(xp,yp) exp(-1j*kk*sqrt( (xp-delmn).^2 + (yp-delmn).^2)) ./ (4*pi*sqrt( (xp-delmn).^2 + (yp-delmn).^2));



%             TxmByn1(mm,nn) = integral2(ggfun,delx*(mm-nn-3/2),delx*(mm-nn-1/2),...
%                                              dely*(mm-nn-3/2),dely*(mm-nn-1/2));
% 
%             TxmByn2(mm,nn) = -integral2(ggfun,delx*(mm-nn-3/2),delx*(mm-nn-1/2),...
%                                              dely*(mm-nn-1/2),dely*(mm-nn+1/2));
% 
%             TxmByn3(mm,nn) = -integral2(ggfun,delx*(mm-nn-1/2),delx*(mm-nn+1/2),...
%                                              dely*(mm-nn-3/2),dely*(mm-nn-1/2));
% 
%             TxmByn4(mm,nn) = integral2(ggfun,delx*(mm-nn-1/2),delx*(mm-nn+1/2),...
%                                             dely*(mm-nn-1/2),dely*(mm-nn+1/2));
%             TxmByn1(mm,nn) = integral2(ggfun,delx*(delmn-3/2),delx*(delmn-1/2),...
%                                              dely*(delmn-3/2),dely*(delmn-1/2));
%                                          
%             TxmByn2(mm,nn) = -integral2(ggfun,delx*(delmn-3/2),delx*(delmn-1/2),...
%                                              dely*(delmn-1/2),dely*(delmn+1/2));
%                                          
%             TxmByn3(mm,nn) = -integral2(ggfun,delx*(delmn-1/2),delx*(delmn+1/2),...
%                                              dely*(delmn-3/2),dely*(delmn-1/2));
%                                          
%             TxmByn4(mm,nn) = integral2(ggfun,delx*(delmn-1/2),delx*(delmn+1/2),...
%                                         dely*(delmn-3/2),dely*(delmn-1/2));
            TxmByn1(mm,nn) = integral2(ggfun,-delx,0,-dely,0);
                                         
            TxmByn2(mm,nn) = -integral2(ggfun,-delx,0,0,dely);
                                         
            TxmByn3(mm,nn) = -integral2(ggfun,0,delx,-dely,0);
                                         
            TxmByn4(mm,nn) = integral2(ggfun,0,delx,0,dely);

                                   
    end
end

BB = TxmByn1 + TxmByn2+TxmByn3+TxmByn4

TxmByn2

%lets get quadratic spline
qqreg1 = linspace(-3*delx/2,-delx/2,100);
qqreg2 = linspace(-delx/2,delx/2,100);
qqreg3 = linspace(delx/2,3*delx/2,100);

qq = [qqreg1, qqreg2, qqreg3];

f1 = 9/8 + 3*qqreg1/(2*delx) +qqreg1.^2/(2*delx.^2);
f2 = 3/4-qqreg2.^2/delx.^2;
f3 = 9/8 - 3*qqreg3/(2*delx) +qqreg3.^2/(2*delx.^2);

ff = [f1,f2,f3];

% figure;plot(qqreg1,f1)
% figure;plot(qqreg2,f2)
% figure;plot(qqreg3,f3)