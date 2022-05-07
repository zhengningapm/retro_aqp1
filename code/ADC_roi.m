clear all;clc;
Dir='I:\data';
subject = importdata([Dir,'\subject.txt']);
roi = importdata([Dir,'\roi.txt']);
ADC_subroi = zeros(length(subject),length(roi));
for sub =1:length(subject)
    for i =1:length(roi)
    Subroi = load_nifti([Dir,'\',subject{sub},'\',roi{i}]);
    Vol_Subroi = Subroi.vol;
    Mean_ADC = sum(sum(sum(Vol_Subroi)))/length(Vol_Subroi(Vol_Subroi>0));
    ADC_subroi(sub,i) = Mean_ADC;
    clear Mean_ADC Subroi Vol_Subroi;
    end
end

Xlsx_data=cell(length(subject)+1,length(roi)+1);

for j=1:length(roi)
Xlsx_data(1,j+1) = cellstr(roi{j});
end

for k=1:length(subject)
Xlsx_data(k+1,1) = cellstr(subject{k});
end

Xlsx_data(2:end,2:end)= num2cell(ADC_subroi);

xlswrite([Dir,'\','ADC_roi.xlsx',],Xlsx_data);