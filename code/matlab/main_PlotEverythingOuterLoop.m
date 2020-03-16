clc; 
close all;
clear all;

flag_outerLoopMode = 1;

%%
% Non-balance point related plots
%%

%Default values
flag_ModeBalancePortrait = 0;

flag_balancePortraitSubject             = 0;    
flag_balancePortraitEnsemble            = 0;
flag_balancePortraitSeatOffEnsemble     = 0;

flag_ModeBalancePointsVsCom0VsCop1      = 0;
flag_ModeAnalyzeBalanceAlong0Across1    = 0;

  flag_fpeTimeSeriesPlotsSubject        = 0;
  flag_fpeTimeSeriesPlotsEnsemble       = 0;
  flag_fpeSeatOffEnsemble               = 0;

  flag_capTimeSeriesPlotsSubject        = 0;
  flag_capTimeSeriesPlotsEnsemble       = 0;
  flag_capSeatOffEnsemble               = 0;

flag_fpeCapErrorTimeSeriesPlotsSubject  = 0;
flag_fpeCapErrorTimeSeriesPlotsEnsemble = 0;
%flag_fpeCapErrorSeatOffEnsemble      

flag_ModeComVel0ComGpVsCop1             = 0;

  flag_comKinematicsTimeSeriesSubject   = 0;
  flag_comKinematicsTimeSeriesEnsemble  = 0;
  flag_comKinematicsSeatOffEnsemble     = 0;

flag_GrfzTimeSeriesPlotsSubject         = 0;
flag_GrfzTimeSeriesPlotsEnsemble        = 0;
flag_GrfzSeatOffEnsemble                = 0;

%%
% Non-balance point related plots
%%
flag_ModeComVel0ComGpVsCop1             = 0;

  flag_comKinematicsTimeSeriesSubject   = 1;
  flag_comKinematicsTimeSeriesEnsemble  = 1;
  flag_comKinematicsSeatOffEnsemble     = 1;

flag_GrfzTimeSeriesPlotsSubject         = 1;
flag_GrfzTimeSeriesPlotsEnsemble        = 1;
flag_GrfzSeatOffEnsemble                = 1;

main_PlotEverything;
close all;

flag_ModeComVel0ComGpVsCop1             = 1;

  flag_comKinematicsTimeSeriesSubject   = 1;
  flag_comKinematicsTimeSeriesEnsemble  = 1;
  flag_comKinematicsSeatOffEnsemble     = 1;

flag_GrfzTimeSeriesPlotsSubject         = 0;
flag_GrfzTimeSeriesPlotsEnsemble        = 0;
flag_GrfzSeatOffEnsemble                = 0;

main_PlotEverything;
close all;

flag_ModeComVel0ComGpVsCop1             = 0;

  flag_comKinematicsTimeSeriesSubject   = 0;
  flag_comKinematicsTimeSeriesEnsemble  = 0;
  flag_comKinematicsSeatOffEnsemble     = 0;


flag_ModeBalancePortrait = 0;

flag_balancePortraitSubject             = 1;    
flag_balancePortraitEnsemble            = 1;
flag_balancePortraitSeatOffEnsemble     = 0;

flag_ModeBalancePointsVsCom0VsCop1      = 0;
flag_ModeAnalyzeBalanceAlong0Across1    = 0;

  flag_fpeTimeSeriesPlotsSubject        = 1;
  flag_fpeTimeSeriesPlotsEnsemble       = 1;
  flag_fpeSeatOffEnsemble               = 1;

  flag_capTimeSeriesPlotsSubject        = 1;
  flag_capTimeSeriesPlotsEnsemble       = 1;
  flag_capSeatOffEnsemble               = 1;

flag_fpeCapErrorTimeSeriesPlotsSubject  = 1;
flag_fpeCapErrorTimeSeriesPlotsEnsemble = 1;


main_PlotEverything;
close all;

flag_fpeCapErrorTimeSeriesPlotsSubject  = 0;
flag_fpeCapErrorTimeSeriesPlotsEnsemble = 0;


flag_ModeBalancePortrait = 1;

flag_balancePortraitSubject             = 1;    
flag_balancePortraitEnsemble            = 1;
flag_balancePortraitSeatOffEnsemble     = 0;

flag_ModeBalancePointsVsCom0VsCop1      = 1;
flag_ModeAnalyzeBalanceAlong0Across1    = 0;

  flag_fpeTimeSeriesPlotsSubject        = 1;
  flag_fpeTimeSeriesPlotsEnsemble       = 1;
  flag_fpeSeatOffEnsemble               = 1;

  flag_capTimeSeriesPlotsSubject        = 1;
  flag_capTimeSeriesPlotsEnsemble       = 1;
  flag_capSeatOffEnsemble               = 1;

main_PlotEverything;
close all;
  
flag_ModeBalancePortrait                = 0;
flag_balancePortraitSubject             = 0;    
flag_balancePortraitEnsemble            = 0;
  
flag_ModeBalancePointsVsCom0VsCop1      = 0;
flag_ModeAnalyzeBalanceAlong0Across1    = 1;

  flag_fpeTimeSeriesPlotsSubject        = 1;
  flag_fpeTimeSeriesPlotsEnsemble       = 1;
  flag_fpeSeatOffEnsemble               = 1;

  flag_capTimeSeriesPlotsSubject        = 1;
  flag_capTimeSeriesPlotsEnsemble       = 1;
  flag_capSeatOffEnsemble               = 1;

main_PlotEverything;
close all;


flag_ModeBalancePointsVsCom0VsCop1      = 1;
flag_ModeAnalyzeBalanceAlong0Across1    = 1;

  flag_fpeTimeSeriesPlotsSubject        = 1;
  flag_fpeTimeSeriesPlotsEnsemble       = 1;
  flag_fpeSeatOffEnsemble               = 1;

  flag_capTimeSeriesPlotsSubject        = 1;
  flag_capTimeSeriesPlotsEnsemble       = 1;
  flag_capSeatOffEnsemble               = 1;

main_PlotEverything;
close all;
