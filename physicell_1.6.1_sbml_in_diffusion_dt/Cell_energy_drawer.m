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

for i = 1:length(OutMatFiles)
    load(OutMatFiles{i})
    PhysiCell_Energy(i) = cells(end);
end
plot(PhysiCell_Energy)
cd ..