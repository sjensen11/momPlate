%generate points for mom plate
%sepearated so I could debug, but I cut and pasted this code into the main
%script already
%plate info
NumCells = 4;
NumEdges = NumCells-1;
lenx = 4;%lambda
leny= 4;%lambda

delx = lenx/NumCells; %a in the book
dely = leny/NumCells; %b in the book

%get points
Bxn_xx =  zeros(1,NumEdges^2);
Bxn_yy = zeros(1,NumEdges^2);

Byn_xx =  zeros(1,NumEdges^2);
Byn_yy = zeros(1,NumEdges^2);

for nn = 0:NumEdges

    Bxn_xx(1+nn*NumEdges:NumEdges+nn*NumEdges) = delx:delx:lenx-delx;
    Bxn_yy(1+nn*NumEdges:NumEdges+nn*NumEdges) =  dely*(nn+1/2)*ones(1,NumCells-1);
    
end

for nn = 0:NumEdges-1
    Byn_xx(1+nn*NumCells:NumEdges+nn*NumCells+1) =  delx/2:delx:lenx-delx/2;
    Byn_yy(1+nn*NumCells:NumEdges+nn*NumCells+1) = delx*(nn+1)*ones(1,NumCells);
end
% 

figure;plot(Bxn_xx,Bxn_yy,'-o')
axis([0 lenx 0 leny])
Bxn_xx
Bxn_yy

figure;plot(Byn_xx,Byn_yy,'-o')
axis([0 lenx 0 leny])
Byn_xx
Byn_yy