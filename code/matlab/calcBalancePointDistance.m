function [errorVec, errorName] = ...
  calcBalancePointDistance( balancePoint, uDirection, nDirection, ...
                            comPosition, copPosition,...
                            flag_ModeBalancePointsVsCom0VsCop1,...
                            flag_ModeAnalyzeBalanceAlong0Across1,...
                            balancePointName)

  errorVec = [];

  switch flag_ModeBalancePointsVsCom0VsCop1
    case 0
      errorVec = balancePoint - comPosition;
      errorName = ['(',balancePointName,'-Com)'];
    case 1
      errorVec = balancePoint - copPosition;
      errorName = ['(',balancePointName,'-Cop)'];
      
    otherwise assert(0)
  end              
  switch(flag_ModeAnalyzeBalanceAlong0Across1)
    case 0
      errorVec = sum(errorVec.*uDirection,2);
      errorName = [errorName,'.u'];
    case 1          
      errorVec = sum(errorVec.*nDirection,2);                  
      errorName = [errorName,'.n'];
      
    otherwise assert(0);
  end