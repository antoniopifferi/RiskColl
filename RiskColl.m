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

KAPPA = 8; %NORMAL/TUMOR in whole population 
TRESHOLD=15;
MAL=1;
NOR=2;
NumType=2; % MAL + NOR
COLL=1;
DENS=2;
NumMeas=2; % COLL + DENS
PERCSCALE=1:1:100; % scale of the percentile
NumPerc=length(PERCSCALE);
BIRADS=[15 50 85 100]; % percentiles for def of Birads class
NumBir=length(BIRADS);
T=readtable('Data.txt','Delimiter','\t');
NumAll=height(T);
iVal=T.SO2<97;
iDens=T.Density>0;
iNotes=strcmp(T.Notes,'');
iBen=strcmp(T.Type,'Benigna');
iNor=strcmp(T.Type,'Normal');
iMal=strcmp(T.Type,'Malignant');
%iNor=iNor+iBen;
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

% calculate percentile
Percentile(:,COLL)=prctile(T.CollNorm(T.Lesion==NOR),PERCSCALE);
Percentile(:,DENS)=prctile(T.DensNorm(T.Lesion==NOR),PERCSCALE);
for ip=1:NumPerc-1
    PCountMal(ip,COLL)=sum(T.CollNorm(T.Lesion==MAL)<=Percentile(ip,COLL));
    PCountMal(ip,DENS)=sum(T.DensNorm(T.Lesion==MAL)<=Percentile(ip,DENS));
    PCountNor(ip,COLL)=sum(T.CollNorm(T.Lesion==NOR)<=Percentile(ip,COLL));
    PCountNor(ip,DENS)=sum(T.DensNorm(T.Lesion==NOR)<=Percentile(ip,DENS));
end
PCountMal(NumPerc,COLL)=length(T.CollNorm(T.Lesion==MAL));
PCountMal(NumPerc,DENS)=length(T.DensNorm(T.Lesion==MAL));
PCountNor(NumPerc,COLL)=length(T.CollNorm(T.Lesion==NOR));
PCountNor(NumPerc,DENS)=length(T.DensNorm(T.Lesion==NOR));


% calc global terms
PCountLes=PCountMal+PCountNor;
HT=(NumMal-PCountMal);
HN=(NumNor-PCountNor);
AT=NumMal;
AN=NumNor;
%PRateMalHigh=100*((NumMal-PCountMal)./(NumNor-PCountNor))./(NumMal./NumNor); % use Num-Count because you want to invert the percentile range
PRateMalHigh=100*(HT./(HT+KAPPA*HN))./(AT./(AT+KAPPA*AN)); % use Num-Count because you want to invert the percentile range
PFractNor=100*PCountNor./PCountLes;
PFractMal=100*PCountMal./PCountLes;
PRateMalHighLow=100*((NumMal-PCountMal)./(NumNor-PCountNor))./(PCountMal./PCountNor); % use Num-Count because you want to invert the percentile range

% calc Birads
BThreshold(:,COLL)=prctile(T.CollNorm(T.Lesion==NOR),BIRADS);
BThreshold(:,DENS)=prctile(T.DensNorm(T.Lesion==NOR),BIRADS);
for ib=1:NumBir-1
    BCountMal(ib,COLL)=sum(T.CollNorm(T.Lesion==MAL)<=BThreshold(ib,COLL));
    BCountMal(ib,DENS)=sum(T.DensNorm(T.Lesion==MAL)<=BThreshold(ib,DENS));
    BCountNor(ib,COLL)=sum(T.CollNorm(T.Lesion==NOR)<=BThreshold(ib,COLL));
    BCountNor(ib,DENS)=sum(T.DensNorm(T.Lesion==NOR)<=BThreshold(ib,DENS));
end
BCountMal(NumBir,COLL)=length(T.CollNorm(T.Lesion==MAL));
BCountMal(NumBir,DENS)=length(T.DensNorm(T.Lesion==MAL));
BCountNor(NumBir,COLL)=length(T.CollNorm(T.Lesion==NOR));
BCountNor(NumBir,DENS)=length(T.DensNorm(T.Lesion==NOR));

for ib=2:NumBir
    for ibb=1:ib-1
        BCountMal(ib,:)=BCountMal(ib,:)-BCountMal(ibb,:);
        BCountNor(ib,:)=BCountNor(ib,:)-BCountNor(ibb,:);
    end
end

%BRiskRatio(1,:)=(BCountMal(4,:)./BCountNor(4,:))./(NumMal./NumNor); % Here we store Class 4 / All
for ib=2:NumBir
    Tb=BCountMal(ib,:);
    Nb=BCountNor(ib,:);
    T1=BCountMal(1,:);
    N1=BCountNor(1,:);
    BRiskRatio(ib,:)=(Tb./(Tb+KAPPA*Nb))./(T1./(T1+KAPPA*N1));
    %BRiskRatio(ib,:)=(BCountMal(ib,:)./(BCountNor(ib,:)))./(BCountMal(1,:)./BCountNor(1,:)); % Here we store Class ib / Class 1
end
BRiskRatio(1,:)=(Tb./(Tb+KAPPA*Nb))./(AT./(AT+KAPPA*AN));

display('Collagen  Density');
display(BRiskRatio);

% FIGURE RISK
figure,
XPerc=[100-PERCSCALE;100-PERCSCALE]';
h=plot(XPerc,PFractMal,'r','LineWidth',2); hold on;
set(h,{'LineStyle'},{'-';':'})
h=plot(XPerc,PFractNor,'b','LineWidth',2);
set(h,{'LineStyle'},{'-';':'})
set(gca,'FontSize',12);
legend('Cancer (Collagen)','Cancer (Density)','Normal (Collagen)','Normal (Density)'), grid on;
xlabel('Threshold (%)'), ylabel('Incidence of Lesion Type (%)');
ylim([0 100]);
title('Risk Factor divided for Total Lesion');

if SAVE_FIG,
    save_figure('FigComparison');
end

% FIGURE RISK RATIO
figure,
h=plot(XPerc,PRateMalHigh,'r','LineWidth',2); hold on;
set(h,{'LineStyle'},{'-';':'})
set(gca,'FontSize',12);
legend('Collagen','Density'), grid on;
xlabel('Percentile (%)'), ylabel('Risk Ratio (high/all) (%)');
ylim([0 600]);
title('Increased of Risk Factor (high/all)');

% FIGURE RISK RATIO HIGH/LOW NEW
figure,
h=plot(XPerc,PRateMalHighLow,'r','LineWidth',2); hold on;
set(h,{'LineStyle'},{'-';':'})
set(gca,'FontSize',12);
legend('Collagen','Density'), grid on;
xlabel('Percentile (%)'), ylabel('Risk Ratio (high/low) (%)');
ylim([0 600]);
title('Increased of Risk Factor (high/low)');

if SAVE_FIG,
    save_figure('FigComparison');
end

%% SCATTER PLOTS
% show individual subjects distribution in Age-matched and not-matched
% conditions
% in the first figure, the black lines corresponds to a cluster of 15% of
% subjects

figure,
set(gca,'FontSize',12);
plot(T.DensNorm(iMal),T.CollNorm(iMal),'dr'), hold on;
plot(T.DensNorm(iNor),T.CollNorm(iNor),'db');
bDens=BThreshold(1:NumBir-1,DENS);
bColl=BThreshold(1:NumBir-1,COLL);
maxDens=max(max(T.DensNorm(iMal)),max(T.DensNorm(iNor)));
maxColl=max(max(T.CollNorm(iMal)),max(T.CollNorm(iNor)));
X=[bDens,bDens];
Y=repmat([0 maxColl],NumBir-1,1);
line(X',Y','Color','k','LineWidth',2);
X=repmat([0 maxDens],NumBir-1,1);
Y=[bColl,bColl];
line(X',Y','Color','k','LineWidth',2);
%grid on;
xlabel('Normalised Density'), ylabel('Normalised Collagen');
legend('Cancer','Normal');
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
