%debugging momPlatescript

%this is wrong older code
% delmn = -1
% low1 =delx*delmn-3*delx/2
% high1 =  delx*delmn-delx/2
% 
% low2 = dely*delmn-dely/2
% high2 = dely*delmn+dely/2
% 
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