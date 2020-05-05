close all
clear
clc

cd output

s=what;
MatFiles = s.mat;
OutMatFiles = MatFiles(contains(MatFiles,'physicell'));
OutMatFiles(1) = [];
OutMatFiles(1) = [];

PhysiCell_Energy = zeros(length(OutMatFiles),1);
PhysiCell_Oxygen = zeros(length(OutMatFiles),1);
PhysiCell_Glucose = zeros(length(OutMatFiles),1);

for i = 1:length(OutMatFiles)
    load(OutMatFiles{i})
    PhysiCell_Energy(i) = cells(28);
    PhysiCell_Oxygen(i) = cells(29);
    PhysiCell_Glucose(i) = cells(30);
end
% plot(PhysiCell_Energy)
% figure(2)
% plot(PhysiCell_Oxygen)
cd ..