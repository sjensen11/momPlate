%behold the unorganized mind of sarah!
%universal constants
kk=2*pi;
netta = 377;
gamma = 1.781072418; %cme-petterson-pg39, (2.13)

phiInc = pi/4;

%plate info
NumCells = 10; %in one direction
NumEdges = NumCells-1;
lenx = 2.23;%lambda
leny= 2.23;%lambda

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
    Byn_xx(1+nn*NumEdges:NumEdges+nn*NumEdges) = delx*(nn+1/2)*ones(1,NumCells-1);
    Byn_yy(1+nn*NumEdges:NumEdges+nn*NumEdges) =  dely:dely:leny-dely;
    
end

% for nn = 0:NumEdges-1
%     Byn_xx(1+nn*NumCells:NumEdges+nn*NumCells+1) =  delx/2:delx:lenx-delx/2;
%     Byn_yy(1+nn*NumCells:NumEdges+nn*NumCells+1) = delx*(nn+1)*ones(1,NumCells);
% end

figure;plot(Bxn_xx,Bxn_yy,'-o');title('Bxn points')
axis([0 lenx 0 leny])

figure;plot(Byn_xx,Byn_yy,'-o'); title('Byn points')
axis([0 lenx 0 leny])

%--------------------------generate matrix equation-----------------------
%need to first generate x, y points
%going to have two sets depending on Ex or Ey
numElt = length(Bxn_xx);
ex = zeros(numElt,1);
ey = zeros(numElt,1);

TxmBxn = zeros(NumEdges*NumCells);
TxmBxnsp = zeros(NumEdges*NumCells);

TxmByn = zeros(NumEdges*NumCells);

TxmByn1 = zeros(NumEdges*NumCells);
TxmByn2 = zeros(NumEdges*NumCells);
TxmByn3 = zeros(NumEdges*NumCells);
TxmByn4 = zeros(NumEdges*NumCells);

TymBxn = zeros(NumEdges*NumCells);
TymByn = zeros(NumEdges*NumCells);
TymBynsp = zeros(NumEdges*NumCells);

for mm = 1:numElt %marches with x
    
    ex(mm) = delx*sin(phiInc) * exp(-1j*kk*(Bxn_xx(mm)*cos(phiInc)+Bxn_yy(mm)*sin(phiInc)))...
                *sin(kk* cos(phiInc)*delx/2);
    ey(mm) = dely*cos(phiInc) * exp(-1j*kk*(Byn_xx(mm)*cos(phiInc)+Byn_yy(mm)*sin(phiInc)))...
                *sin(kk*sin(phiInc)*delx/2);
    for nn = 1:numElt %marches with y

            %define function
            %green's function shifted to xm-xn and ym-yn
            %Bxn points then Byn points
            ggTxnBxn = @(xp,yp) exp(-1j*kk*sqrt( (xp+Bxn_xx(mm)-Bxn_xx(nn)).^2 + (yp+Bxn_yy(mm)-Bxn_yy(nn)).^2))...
                                  ./ (4*pi*sqrt( (xp+Bxn_xx(mm)-Bxn_xx(nn)).^2 + (yp+Bxn_yy(mm)-Bxn_yy(nn)).^2));
                            
            ggTynByn = @(xp,yp) exp(-1j*kk*sqrt( (xp+Byn_xx(mm)-Byn_xx(nn)).^2 + (yp+Byn_yy(mm)-Byn_yy(nn)).^2))...
                                  ./ (4*pi*sqrt( (xp+Byn_xx(mm)-Byn_xx(nn)).^2 + (yp+Byn_yy(mm)-Byn_yy(nn)).^2));
            
            ggTxmByn = @(xp,yp) exp(-1j*kk*sqrt( (xp+Bxn_xx(mm)-Byn_xx(nn)).^2 + (yp+Bxn_yy(mm)-Byn_yy(nn)).^2))...
                                  ./ (4*pi*sqrt( (xp+Bxn_xx(mm)-Byn_xx(nn)).^2 + (yp+Bxn_yy(mm)-Byn_yy(nn)).^2));
                                  
            ggTymBxn = @(xp,yp) exp(-1j*kk*sqrt( (xp+Byn_xx(mm)-Bxn_xx(nn)).^2 + (yp+Byn_yy(mm)-Bxn_yy(nn)).^2))...
                                  ./ (4*pi*sqrt( (xp+Byn_xx(mm)-Bxn_xx(nn)).^2 + (yp+Byn_yy(mm)-Bxn_yy(nn)).^2));
            %spleen
            %Bxn points first
            ggpp1xxBxn = @(xp,yp) (9/8+xp.*3/(2*delx) + (xp.^2)/(2*delx^2)) .* ggTxnBxn(xp,yp);      
            ggpp2xxBxn = @(xp,yp) (3/4-xp.^2./ delx^2).* ggTxnBxn(xp,yp);
            ggpp3xxBxn = @(xp,yp) (9/8-xp.*3/(2*delx) + (xp.^2)/(2*delx^2)) .* ggTxnBxn(xp,yp);  
                            
            
            ggpp1yyBxn = @(xp,yp) (9/8+yp.*3/(2*dely) +yp.^2./(2*dely^2)) .*  ggTxnBxn(xp,yp);    
            ggpp2yyBxn = @(xp,yp) (3/4-yp.^2./ dely^2).* ggTxnBxn(xp,yp);
            ggpp3yyBxn = @(xp,yp) (9/8-3*yp./(2*dely) +yp.^2/(2*dely.^2)).* ggTxnBxn(xp,yp);
            
            %now Byn points
            ggpp1xxByn = @(xp,yp) (9/8+xp.*3/(2*delx) + (xp.^2)/(2*delx^2)) .* ggTynByn(xp,yp);      
            ggpp2xxByn = @(xp,yp) (3/4-xp.^2./ delx^2).* ggTynByn(xp,yp);
            ggpp3xxByn = @(xp,yp) (9/8-xp.*3/(2*delx) + (xp.^2)/(2*delx^2)) .* ggTynByn(xp,yp);  
                            
            
            ggpp1yyByn = @(xp,yp) (9/8+yp.*3/(2*dely) +yp.^2./(2*dely^2)) .*  ggTynByn(xp,yp);    
            ggpp2yyByn = @(xp,yp) (3/4-yp.^2./ dely^2).* ggTynByn(xp,yp);
            ggpp3yyByn = @(xp,yp) (9/8-3*yp./(2*dely) +yp.^2/(2*dely.^2)).* ggTynByn(xp,yp);
            
            %do integrals
            %first convolution
            TxmBxn(mm,nn) = integral2(ggTxnBxn,-delx*3/2,-delx/2,-dely/2,dely/2)...
                             +2*integral2(ggTxnBxn,-delx/2,delx/2,-dely/2,dely/2)...
                             +integral2(ggTxnBxn,delx/2,delx*3/2,-dely/2,dely/2);
            
            TxmBxnsp(mm,nn) = integral2(ggpp1xxBxn,-delx*3/2,-delx/2,-dely/2,dely/2)...
                             +integral2(ggpp2xxBxn,-delx/2,delx/2,-dely/2,dely/2)...
                             +integral2(ggpp3xxBxn,delx/2,delx*3/2,-dely/2,dely/2);
            
            %second and third convolution
            %TxmByn== TynBxn expect for the 1/a vs 1/b... so only get the
            %matrix once... 
            TxmByn1(mm,nn) = integral2(ggTxmByn,-delx,0,-dely,0);
            TxmByn2(mm,nn) = -integral2(ggTxmByn,-delx,0,0,dely);
            TxmByn3(mm,nn) = -integral2(ggTxmByn,0,delx,-dely,0);
            TxmByn4(mm,nn) = integral2(ggTxmByn,0,delx,0,dely);
%             TxmByn(mm,nn) = integral2(ggTxmByn,-delx,0,-dely,0)-integral2(ggTxmByn,-delx,0,0,dely)...
%                            -integral2(ggTxmByn,0,delx,-dely,0) +integral2(ggTxmByn,0,delx,0,dely);
%                        
%             %third convolution to test
            TymBxn(mm,nn) = integral2(ggTymBxn,-delx,0,-dely,0)-integral2(ggTymBxn,-delx,0,0,dely)...
                           -integral2(ggTymBxn,0,delx,-dely,0) +integral2(ggTymBxn,0,delx,0,dely);
%             %fourth convolution...if you're going to do a square you can
%             %throw this out and just use the first convolution
            TymByn(mm,nn) = integral2(ggTynByn,-delx/2,delx/2,-dely*3/2,-dely/2)...
                         +2*integral2(ggTynByn,-delx/2,delx/2,-dely/2,dely/2)...
                         +  integral2(ggTynByn,-delx/2,delx/2,delx/2,delx*3/2);
                     
            TymBynsp(mm,nn) = integral2(ggpp1yyByn,-delx/2,delx/2,-dely*3/2,dely/2)...
                             +integral2(ggpp2yyByn,-delx/2,delx/2,-dely/2,dely/2)...
                             +integral2(ggpp3yyByn,-delx/2,delx/2,dely/2,dely*3/2);
                            
            

                                   
    end
end
TxmByn = TxmByn1+TxmByn2+TxmByn3+TxmByn4
AA = TxmBxn/delx+kk^2*delx*TxmBxnsp;
BB = TxmByn/delx;
CC = TxmByn/dely;
DD = TymByn/dely+kk^2*dely*TymBynsp;
% 
EE = [ex;ey]
ZZ = (-netta/(1j*kk))*[AA BB; CC DD];

JJ = ZZ\EE


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