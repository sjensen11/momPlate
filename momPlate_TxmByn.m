%behold the unorganized mind of sarah!
%made separate script to debug TxmByn part
%universal constants
kk=2*pi;
netta = 377;
gamma = 1.781072418; %cme-petterson-pg39, (2.13)

phiInc = pi/4;
thetaInc = pi/2;
%plate info
NumCells = 4; %in one direction
NumEdges = NumCells-1;
lenx = 1%lambda
leny= 1;%lambda

delx = lenx/NumCells; %a in the book
dely = leny/NumCells; %b in the book
%-----------------------Generate points---------------------------------
Bxn_xx =  zeros(1,NumEdges^2);
Bxn_yy = zeros(1,NumEdges^2);

Byn_xx =  zeros(1,NumEdges^2);
Byn_yy = zeros(1,NumEdges^2);

for nn = 0:NumEdges

    Bxn_xx(1+nn*NumEdges:NumEdges+nn*NumEdges) = delx:delx:lenx-delx;
    Bxn_yy(1+nn*NumEdges:NumEdges+nn*NumEdges) =  dely*(nn+1/2)*ones(1,NumCells-1);
    
%     Byn_xx(1+nn*NumEdges:NumEdges+nn*NumEdges) = delx*(nn+1/2)*ones(1,NumCells-1);
%     Byn_yy(1+nn*NumEdges:NumEdges+nn*NumEdges) =  dely:dely:leny-dely;
    
end

for nn = 0:NumEdges-1
    Byn_xx(1+nn*NumCells:NumEdges+nn*NumCells+1) =  delx/2:delx:lenx-delx/2;
    Byn_yy(1+nn*NumCells:NumEdges+nn*NumCells+1) = delx*(nn+1)*ones(1,NumCells);
end

figure;plot(Bxn_xx,Bxn_yy,'-o');title('Bxn points')
axis([0 lenx 0 leny])

figure;plot(Byn_xx,Byn_yy,'-o'); title('Byn points')
axis([0 lenx 0 leny])

%--------------------------generate matrix equation-----------------------
%need to first generate x, y points
%going to have two sets depending on Ex or Ey
numElt = length(Bxn_xx);

TxmByn = zeros(NumEdges*NumCells);

TxmByn1 = zeros(NumEdges*NumCells);
TxmByn2 = zeros(NumEdges*NumCells);
TxmByn3 = zeros(NumEdges*NumCells);
TxmByn4 = zeros(NumEdges*NumCells);


for mm = 1:numElt %marches with x
    

    for nn = 1:numElt %marches with y

            %define function

            ggTxmByn = @(xp,yp) exp(-1j*kk*sqrt( (xp+Bxn_xx(mm)-Byn_xx(nn)).^2 + (yp+Bxn_yy(mm)-Byn_yy(nn)).^2))...
                                  ./ (4*pi*sqrt( (xp+Bxn_xx(mm)-Byn_xx(nn)).^2 + (yp+Bxn_yy(mm)-Byn_yy(nn)).^2));
            
            %do integrals
%             distx = Bxn_xx(mm)-Byn_xx(nn)
%             disty = Bxn_yy(mm)-Byn_yy(nn)
%             figure; 
%             subplot(1,2,1); fplot(@(xp) real(ggTxmByn(xp,0)),[-5,5]); title('y=0')
%             subplot(1,2,2); fplot(@(yp) real(ggTxmByn(0,yp)),[-5,5]); title('x=0')
            %second and third convolution
            %TxmByn== TynBxn expect for the 1/a vs 1/b... so only get the
            %matrix once... 
            TxmByn1(mm,nn) = integral2(ggTxmByn,-delx,0,-dely,0);
            TxmByn2(mm,nn) = -integral2(ggTxmByn,-delx,0,0,dely);
            TxmByn3(mm,nn) = -integral2(ggTxmByn,0,delx,-dely,0);
            TxmByn4(mm,nn) = integral2(ggTxmByn,0,delx,0,dely);
            TxmByn = TxmByn1+TxmByn2+TxmByn3+TxmByn4
%             TxmByn(mm,nn) = integral2(ggTxmByn,-delx,0,-dely,0)-integral2(ggTxmByn,-delx,0,0,dely)...
%                            -integral2(ggTxmByn,0,delx,-dely,0) +integral2(ggTxmByn,0,delx,0,dely);
%                        

                            
            

                                   
    end
end
TxmByn = TxmByn1+TxmByn2+TxmByn3+TxmByn4
% AA = TxmBxn/delx+kk^2*delx*TxmBxnsp;
% BB = TxmByn/delx;
% CC = TxmByn/dely;
% DD = TymByn/dely+kk^2*dely*TymBynsp;
% 
% EE = [ex;ey]
% ZZ = (-netta/(1j*kk))*[AA BB; CC DD];

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