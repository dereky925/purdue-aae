clear;clc;close all

% Max possible grade as of Saturday March 8 2025 for this class

maxgradepossible = 100 - (100-37.5)/100*10 - (100-62)/100*10 - (100-86)/100*10 - (100-99)/100*10

FourAssigmentPts = (37.5)/100*10 + (62)/100*10 + (86)/100*10 + (99)/100*10

ptsNeededForB = 83 - FourAssigmentPts

AssignmentsLeft = 6;

avgGradeNeededOnRemainingAssignments = ptsNeededForB / AssignmentsLeft * 10

% To reach at least an 80% final grade, you should aim for roughly 
% an 86% average on each of the remaining assignments and exams