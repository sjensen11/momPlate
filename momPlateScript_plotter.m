%plotting values from momPlateScript, so I have x,y,J in an array 
%numElt,NumEdge,NumCell should already be defined
%Bxn so row = NumCell, col = NumEdge
XX = zeros(NumCells,NumEdges);
YY = zeros(NumCells,NumEdges);
JJMat = zeros(NumCells,NumEdges);
for row = 1:NumCells
    for col = 1:NumEdges
        XX(row,col) = Bxn_xx(col+NumEdges*(row-1));
        YY(row,col) = Bxn_yy(col+NumEdges*(row-1));
        JJMat(row,col) = JJ( col+NumEdges*(row-1));
    end
end

figure;surf(XX,YY,abs(JJMat))