%% Direct Risk Assessment by Collagen - ECBO Data
%
% Data Repository = DataECBO_AllType.xls
% used for generating Proceedings Data
%
%% LOAD DATA
% load file data and apply these inclusion criteria:
% * S02<95
% * Density > 0
%
close all;
clear all;
SAVE_FIG=0;

TRESHOLD=15;
MAL=1;
NOR=2;
T=readtable('Data.txt','Delimiter','\t');
NumAll=height(T);
iVal=T.SO2<97;
iDens=T.Density>0;
iNotes=strcmp(T.Notes,'');
%iBen=strcmp(T.Type,'Benigna');
iNor=strcmp(T.Type,'Normal');
%iNor=iNor+iBen;
iMal=strcmp(T.Type,'Malignant');
iNor=(iNor+iDens+iVal+iNotes)==4;
iMal=(iMal+iDens+iVal+iNotes)==4;
T.Lesion=zeros(NumAll,1);
T.Lesion(iMal)=MAL;
T.Lesion(iNor)=NOR;
NumMal=sum(iMal);
NumNor=sum(iNor);
NumLes=NumMal+NumNor;
iAll=1:NumAll;

figure,
plot(T.Density(iMal),T.Collagen(iMal),'dr'), hold on;
plot(T.Density(iNor),T.Collagen(iNor),'db');
xlabel('Density (%)'), ylabel('Collagen mg/mm^{3}'), title('Not-Normalised Data');


%% NORMALISE DATA
% show linear interpolation to derive age-normalised data using the
% formula:
% * Norm = quantity/linearinterp
%

[FitColl,stats]=robustfit(T.Age,T.Collagen);
[FitDens,stats]=robustfit(T.Age,T.Density);
T.CollAge=FitColl(1)+FitColl(2)*T.Age;
T.DensAge=FitDens(1)+FitDens(2)*T.Age;
T.CollNorm=T.Collagen./T.CollAge;
T.DensNorm=T.Density./T.DensAge;
% T.CollNorm=(T.Collagen-T.CollAge)./T.CollAge;
% T.DensNorm=(T.Density-T.DensAge)./T.DensAge;
% T.CollNorm=T.Collagen;
% T.DensNorm=T.Density;



figure,
set(gca,'FontSize',12);
plot(T.Age,T.Collagen,'d'); hold on;
plot(T.Age,T.CollAge,'-r','LineWidth',2);
xlabel('Age'), ylabel('Collagen mg/mm^{3}');
if SAVE_FIG, save_figure('AgeCorrectionCollagen');
else title('Age Dependence on Collagen (AllData)'); end

figure,
set(gca,'FontSize',12);
plot(T.Age,T.Density,'d'); hold on;
plot(T.Age,T.DensAge,'-r','LineWidth',2);
xlabel('Age'), ylabel('Density (%)');
if SAVE_FIG, save_figure('AgeCorrectionDensity');
else title('Age Dependence on Density (AllData)'); end

%% SCATTER PLOTS
% show individual subjects distribution in Age-matched and not-matched
% conditions
% in the first figure, the black lines corresponds to a cluster of 15% of
% subjects

figure,
set(gca,'FontSize',12);
plot(T.DensNorm(iMal),T.CollNorm(iMal),'dr'), hold on;
plot(T.DensNorm(iNor),T.CollNorm(iNor),'db');
plot([0 7],[1.61 1.61],'k-','LineWidth',2);
plot([1.6 1.6],[0 3],'k-','LineWidth',2);
grid on;
xlabel('Normalised Density'), ylabel('Normalised Collagen');
legend('Cancer','Normal');
annotation('textbox',[0.317 0.768 0.115 0.076],'String',{'high density'},'FontWeight','bold','FontSize',14,'FitBoxToText','off','EdgeColor','none');
annotation('textbox',[0.525 0.560 0.311 0.076],'String',{'high collagen'},'FontWeight','bold','FontSize',14,'FitBoxToText','off','EdgeColor','none');
if SAVE_FIG, save_figure('ScatterPlot');
else title('Age-Normalised Data'); end

figure,
plot(T.Age(iMal),T.Collagen(iMal),'dr'), hold on;
plot(T.Age(iNor),T.Collagen(iNor),'db');
xlabel('Age'), ylabel('Collagen mg/mm^{3}'), title('Not-Normalised Data');

figure,
plot(T.Age(iMal),T.CollNorm(iMal),'dr'), hold on;
plot(T.Age(iNor),T.CollNorm(iNor),'db');
xlabel('Age'), ylabel('Collagen mg/mm^{3}'), title('Age-Normalised Data');

%% CLASSIFICATION
% identify classification scores by sorting the subjects according to the
% normalised quantity and calculating the ratio of increase in risk as a
% function of the percentage of subjects classified at high risk, according
% to the formula:
% 
% * RiskFactor = (NumMalignantClass/NumLesionClass)/(NumMalignantTor/NumLesionTot)
%
% the same is calculated starting from the lowest values
%

% OLD CLASSIFICATION
for im=1:2,

    if im==1, [CollNormSort,iSort]=sort(T.CollNorm,'descend'); end
    if im==2, [CollNormSort,iSort]=sort(T.DensNorm,'descend'); end
    
for ip=1:NumAll,
    CountMal=sum(T.Lesion(iSort(1:ip))==MAL);
    CountNor=sum(T.Lesion(iSort(1:ip))==NOR);
    CountLes=CountMal+CountNor;
    ACountMal(ip,im)=CountMal;
    ACountNor(ip,im)=CountNor;
    ACountLes(ip,im)=CountLes;
    ATreshold(ip,im)=100*CountLes/NumLes;
    if CountLes>0,
        RateMalHigh(CountLes)=100*(CountMal/CountLes)/(NumMal/NumLes);
        Treshold(CountLes)=100*CountNor/NumNor;
        %Treshold(CountLes)=100*CountLes/NumLes;
        FractNor(CountLes)=100*CountNor/CountLes;
        FractMal(CountLes)=100*CountMal/CountLes;
    end
end

% figure,
% plot(Treshold,[RateMalHigh;RateMalLow]), grid on, hold on;
% legend('Higher Value','Lower Value'), grid on;
% xlabel('Treshold (%)'), ylabel('Risk Factor');
% ylim([0 200]);
% if im==1, title('Increased Risk Factor for COLLAGEN (100% =same risk)'); end
% if im==2, title('Increased Risk Factor for DENSITY (100% =same risk)'); end
% 
% figure,
% plot(Treshold,FractMal,'-r',Treshold,FractNor,':b','LineWidth',2), grid on, hold on;
% set(gca,'FontSize',12);
% legend('Cancer','Normal'), grid on;
% xlabel('Threshold (%)'), ylabel('Incidence of Lesion Type (%)');
% ylim([0 100]);
% %if im==1, title('Increased malignant lesion occurence for COLLAGEN'); end
% %if im==2, title('Increased malignant lesion occurence for DENSITY'); end
% 
% if SAVE_FIG,
%     save_figure(['FigFract_' num2str(im)]);
% else
%     if im==1, title('Increased malignant lesion occurence for COLLAGEN'); end
%     if im==2, title('Increased malignant lesion occurence for DENSITY'); end
% end

% FIGURE RISK
if im==1, figure, plot(Treshold,FractMal,'-r',Treshold,FractNor,'-b','LineWidth',2), grid on, hold on; 
else plot(Treshold,FractMal,':r',Treshold,FractNor,':b','LineWidth',2), grid on, hold on; end
set(gca,'FontSize',12);
if im==2, legend('Cancer (Collagen)','Normal (Collagen)','Cancer (Density)','Normal (Density)'), grid on; end
xlabel('Threshold (%)'), ylabel('Incidence of Lesion Type (%)');
ylim([0 100]);
%if im==1, title('Increased malignant lesion occurence for COLLAGEN'); end
%if im==2, title('Increased malignant lesion occurence for DENSITY'); end


if SAVE_FIG,
    if im==2, save_figure('FigComparison'); end
else
    if im==1, title('Increased malignant lesion occurence for COLLAGEN'); end
    if im==2, title('Increased malignant lesion occurence for DENSITY'); end
end

end

for im=1:2,
    iTreshold(im)=max(find(ATreshold(:,im)<TRESHOLD))
    ValueMal(im)=ACountMal(iTreshold(im),im)
    ValueNor(im)=ACountNor(iTreshold(im),im)
end
    

% % OLD CLASSIFICATION
% for im=1:2,
% 
%     if im==1, [CollNormSort,iSort]=sort(T.CollNorm,'descend'); end
%     if im==2, [CollNormSort,iSort]=sort(T.DensNorm,'descend'); end
%     
% for ip=1:NumAll,
%     CountMal=sum(T.Lesion(iSort(1:ip))==MAL);
%     CountNor=sum(T.Lesion(iSort(1:ip))==NOR);
%     CountLes=CountMal+CountNor;
%     ACountMal(ip,im)=CountMal;
%     ACountNor(ip,im)=CountNor;
%     ACountLes(ip,im)=CountLes;
%     ATreshold(ip,im)=100*CountLes/NumLes;
%     if CountLes>0,
%         RateMalHigh(CountLes)=100*(CountMal/CountLes)/(NumMal/NumLes);
%         Treshold(CountLes)=100*CountNor/NumNor;
%         %Treshold(CountLes)=100*CountLes/NumLes;
%         FractNor(CountLes)=100*CountNor/CountLes;
%         FractMal(CountLes)=100*CountMal/CountLes;
%     end
% end
% 
% % FIGURE RISK
% if im==1, figure, plot(Treshold,FractMal,'-r',Treshold,FractNor,'-b','LineWidth',2), grid on, hold on; 
% else plot(Treshold,FractMal,':r',Treshold,FractNor,':b','LineWidth',2), grid on, hold on; end
% set(gca,'FontSize',12);
% if im==2, legend('Cancer (Collagen)','Normal (Collagen)','Cancer (Density)','Normal (Density)'), grid on; end
% xlabel('Threshold (%)'), ylabel('Incidence of Lesion Type (%)');
% ylim([0 100]);
% 
% 
% if SAVE_FIG,
%     if im==2, save_figure('FigComparison'); end
% else
%     if im==1, title('Increased malignant lesion occurence for COLLAGEN'); end
%     if im==2, title('Increased malignant lesion occurence for DENSITY'); end
% end
% 
% end
% 
% for im=1:2,
%     iTreshold(im)=max(find(ATreshold(:,im)<TRESHOLD))
%     ValueMal(im)=ACountMal(iTreshold(im),im)
%     ValueNor(im)=ACountNor(iTreshold(im),im)
% end
    
%% NEW CLASSIFICATION

% get sorted index (Array are not touched)
for im=1:2,

    if im==1, [CollNormSort,iSort]=sort(T.CollNorm,'descend'); end
    if im==2, [CollNormSort,iSort]=sort(T.DensNorm,'descend'); end
    
for ip=1:NumAll,
    ACountMal(ip,im)=sum(T.Lesion(iSort(1:ip))==MAL);
    ACountNor(ip,im)=sum(T.Lesion(iSort(1:ip))==NOR);
end
end
ACountLes=ACountMal+ACountNor;
ATreshold=100*ACountNor/NumLes;
ARateMalHigh=100*(ACountMal./ACountNor)./(NumMal./NumNor);
AFractNor=100*ACountNor./ACountLes;
AFractMal=100*ACountMal./ACountLes;

% FIGURE RISK
figure,
h=plot(ATreshold,AFractMal,'r','LineWidth',2); hold on;
set(h,{'LineStyle'},{'-';':'})
h=plot(ATreshold,AFractNor,'b','LineWidth',2);
set(h,{'LineStyle'},{'-';':'})
set(gca,'FontSize',12);
legend('Cancer (Collagen)','Normal (Collagen)','Cancer (Density)','Normal (Density)'), grid on;
xlabel('Threshold (%)'), ylabel('Incidence of Lesion Type (%)');
ylim([0 100]);
title('Risk Factor divided for Total Lesion');

if SAVE_FIG,
    save_figure('FigComparison');
end

% FIGURE RISK RATIO
figure,
h=plot(ATreshold,ARateMalHigh,'r','LineWidth',2); hold on;
set(h,{'LineStyle'},{'-';':'})
set(gca,'FontSize',12);
legend('Collagen','Density'), grid on;
xlabel('Threshold (%)'), ylabel('Risk Ratio (high/all) (%)');
ylim([0 600]);
title('Increased of Risk Factor (high/all)');

if SAVE_FIG,
    save_figure('FigComparison');
end



for im=1:2,
    iTreshold(im)=max(find(ATreshold(:,im)<TRESHOLD))
    ValueMal(im)=ACountMal(iTreshold(im),im)
    ValueNor(im)=ACountNor(iTreshold(im),im)
end

