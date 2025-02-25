%This scipt reads out the saved mat data structure with all animals and
%injection inlcuded and plots figures/ performs analysis 
%MATLAB 2024b needed 
%% Dependendencies (install / add before starting)
%1) uipickfiles: https://uk.mathworks.com/matlabcentral/fileexchange/10867-uipickfiles-uigetfile-on-steroids
%2) plotSpread: https://uk.mathworks.com/matlabcentral/fileexchange/37105-plot-spread-points-beeswarm-plot
%3) ploterr: https://uk.mathworks.com/matlabcentral/fileexchange/22216-ploterr
%4) 200 colormap: https://uk.mathworks.com/matlabcentral/fileexchange/120088-200-colormap 
%5) offsetAxes: https://uk.mathworks.com/matlabcentral/fileexchange/52351-offsetaxes-ax
%6) linspecer: https://uk.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-colormap
%7) cmap: https://uk.mathworks.com/matlabcentral/fileexchange/42450-custom-colormap 
%Matlab structure: Data_SW_20-02-Apr-2024.mat


%% Folder where mat structure, where save folders, csvs for boli quantificaion and GAD overalap are (CHANGE ACCORDINGLY)
%CHANGE ACCORDINGLY
clc;clear all; close all;
str   = 'C:\Users\simonw\S1-V1 interaction Dropbox\Simon_Manuel\Manuscripts\L6_dominance\Version_of_record\Data_structures';
% load structure (mat file) using uipickfiles.m 1)
folder_list = uipickfiles('FilterSpec',str);load(char(folder_list));
%folder location where to save figure panels, CHANGE ACCORDINGLY
save_folder='C:\Users\simonw\S1-V1 interaction Dropbox\Simon_Manuel\Manuscripts\L6_dominance\Version_of_record\Vectorgraphics';

%% Injection BOLI & GAD overlap 

%Injection BOLI
%CHANGE ACCORDINGLY
% Injection quantification: spillover and bolus
injection_info = readtable([str '\injection_quant.csv']);
%read out nr, abbreviation and long cortex name of ABI
inj_perc=table2cell(injection_info);
injection_bolus = readtable([str '\injection_bolus.csv']);
inj_bol=table2cell(injection_bolus);

%GAD
%CHANGE ACCORDINGLY
GAD_info = readtable([str '\GAD_overlap.csv']);
%read out GAD overlap
gad_perc=table2cell(GAD_info);

%% Perform initial calculations and extraction for all mice together and/or specific for injection 
%arial as default
set(0, 'DefaultAxesFontName', 'Arial'); 
% Color scheme for the plots
v1_color=[173 7 227]./256;%magenta V1
s1_color=[5 190 120]./256;%green S1
m1_color=[87 87 244]./256;%blue M1
% Load csv with cortex names and assign 6 larger group areas: Frontal, Lateral, Somatomotor, Visual, Medial, Auditory
cortex_info = readtable([str '\names_abbreviations.csv']);
%read out nr, abbreviation and long cortex name of ABI
cortex_names=table2cell(cortex_info);
%change here index of categories
frontal_idx=[1:8];lateral_idx=[9:16 44:45];somamo_idx=[17:26];visual_idx=[27:36];medial_idx=[37:39];aud_idx=[40:43];
%injected brain areas VISp, SSp, MOp, (AUDp)
visp_idx=[31];ssp_idx=[18];mop_idx=[25];audp_idx=[41];
%L4 missing areas
l4miss_idx=[1:11 15:16 25:26 37 39 44 45];
% Find the cells in a given area in their given layer RETRO
% select all retro injection regardeless of their injection type using
% function cell_selecter.m)
temp1=[];temp1=cell_selecter(data,'virustype',1);
retro=[];retro=find(temp1==1);
%call function anatomy_cellnr.m)
[i_animal c_animal i_areas_animal c_areas_animal i_areas_animalf c_areas_animalf  i_areas_animalcd c_areas_animalcd] = anatomy_cellnr(data,retro,cortex_names);
% Extract fraction per injeciton area M1 / S1 / V1
temp2=[];temp2=cell_selecter(data,'virustype',1,'area',1);
temp3=[];temp3=cell_selecter(data,'virustype',1,'area',2);
temp4=[];temp4=cell_selecter(data,'virustype',1,'area',3);
%MOs M2 area
temp5=[];temp5=cell_selecter(data,'virustype',10,'area',10);
v1r=[];v1r=find(temp2==1);s1r=[];s1r=find(temp3==1);m1r=[];m1r=find(temp4==1);m2r=[];m2r=find(temp5==1);
%call function anatomy_cellnr.m
%V1, S1 and M1
[i_v1a c_v1a i_v1aa c_v1aa i_v1aaf c_v1aaf i_v1aacd c_v1aacd] = anatomy_cellnr(data,v1r,cortex_names);
idx_nanv1=find(sum(isnan(i_v1aa(3,:,:)),3)>0 & sum(isnan(i_v1aa(3,:,:)),3)<size(i_v1aa,3));
i_v1aaf(3,idx_nanv1,:)=NaN;i_v1aa(3,idx_nanv1,:)=NaN;i_v1aacd(3,idx_nanv1,:)=NaN;
c_v1aaf(3,idx_nanv1,:)=NaN;c_v1aa(3,idx_nanv1,:)=NaN;c_v1aacd(3,idx_nanv1,:)=NaN;
[i_s1a c_s1a i_s1aa c_s1aa i_s1aaf c_s1aaf i_s1aacd c_s1aacd] = anatomy_cellnr(data,s1r,cortex_names);
idx_nans1=find(sum(isnan(i_s1aa(3,:,:)),3)>0 & sum(isnan(i_s1aa(3,:,:)),3)<size(i_s1aa,3));
i_s1aaf(3,idx_nans1,:)=NaN;i_s1aa(3,idx_nans1,:)=NaN;i_s1aacd(3,idx_nans1,:)=NaN;
c_s1aaf(3,idx_nans1,:)=NaN;c_s1aa(3,idx_nans1,:)=NaN;c_v1aacd(3,idx_nans1,:)=NaN;
[i_m1a c_m1a i_m1aa c_m1aa i_m1aaf c_m1aaf i_m1aacd c_m1aacd] = anatomy_cellnr(data,m1r,cortex_names); 
idx_nanm1=find(sum(isnan(i_m1aa(3,:,:)),3)>0 & sum(isnan(i_m1aa(3,:,:)),3)<size(i_m1aa,3));
i_m1aaf(3,idx_nanm1,:)=NaN;i_m1aa(3,idx_nanm1,:)=NaN;i_m1aacd(3,idx_nanm1,:)=NaN;
c_m1aaf(3,idx_nanm1,:)=NaN;c_m1aa(3,idx_nanm1,:)=NaN;c_m1aacd(3,idx_nanm1,:)=NaN;

[i_m2a c_m2a i_m2aa c_m2aa i_m2aaf c_m2aaf i_m2aacd c_m2aacd] = anatomy_cellnr(data,m2r,cortex_names); 
% calculate indexes using anatomy_indexcalc.m : ILN, L6ab, L6a, hindex ALL INJECTIONS TOGETHER
    %1= ILN
    %2= L6ab dominance index
    %3= L6a dominacne index
    %4= h-index
    %5= L6 pure index 
%overal indexes
[i_index] = anatomy_indexcalc(i_animal);
[c_index] = anatomy_indexcalc(c_animal);
%indexes for specific injection area: V1, S1, M1
[iv1_index] = anatomy_indexcalc(i_v1a);
[cv1_index] = anatomy_indexcalc(c_v1a);
[is1_index] = anatomy_indexcalc(i_s1a);
[cs1_index] = anatomy_indexcalc(c_s1a);
[im1_index] = anatomy_indexcalc(i_m1a);
[cm1_index] = anatomy_indexcalc(c_m1a);

%% criteria when selecting cellcounts per 45 brain areas, and mice 
%Decided with TWM on the 11.3.2024: 
%when dividing into 45 areas we exclude an area with 10 or less cells and
%this need to be present at least in 3 animals
%cell threshold 
thr_count=10;
%Animal threshold 
an_count=3;
%first output: absolute numbers per layer per area per animal
%second output: relative numbers per layer per area per animal
[i_v1aam i_v1aafm]=zero_area_layers(i_v1aa,i_v1aaf,thr_count,an_count);
[c_v1aam c_v1aafm]=zero_area_layers(c_v1aa,c_v1aaf,thr_count,an_count);
[i_s1aam i_s1aafm]=zero_area_layers(i_s1aa,i_s1aaf,thr_count,an_count);
[c_s1aam c_s1aafm]=zero_area_layers(c_s1aa,c_s1aaf,thr_count,an_count);
[i_m1aam i_m1aafm]=zero_area_layers(i_m1aa,i_m1aaf,thr_count,an_count);
[c_m1aam c_m1aafm]=zero_area_layers(c_m1aa,c_m1aaf,thr_count,an_count);

%absolute across animals per injection after criteria 
v1a_tm=[nanmean(nansum(i_v1aam(:,:,:)),3)' nanmean(nansum(c_v1aam(:,:,:)),3)'];
s1a_tm=[nanmean(nansum(i_s1aam(:,:,:)),3)' nanmean(nansum(c_s1aam(:,:,:)),3)'];
m1a_tm=[nanmean(nansum(i_m1aam(:,:,:)),3)' nanmean(nansum(c_m1aam(:,:,:)),3)'];

%absolute numbers per animal for ipsi and contra 
v1a_tcam=squeeze(nansum(c_v1aam(:,:,:)));
v1a_tiam=squeeze(nansum(i_v1aam(:,:,:)));
s1a_tcam=squeeze(nansum(c_s1aam(:,:,:)));
s1a_tiam=squeeze(nansum(i_s1aam(:,:,:)));
m1a_tcam=squeeze(nansum(c_m1aam(:,:,:)));
m1a_tiam=squeeze(nansum(i_m1aam(:,:,:)));

%in case ine usese the homotopi contralteral set to zero 
%set homotopic to NaN
v1a_tcam2=v1a_tcam;
s1a_tcam2=s1a_tcam;
m1a_tcam2=m1a_tcam;
%absolute with NaN
v1a_tcam2(visp_idx,:)=ones(1,size(v1a_tcam2,2))*NaN;
s1a_tcam2(ssp_idx,:)=ones(1,size(s1a_tcam2,2))*NaN;
m1a_tcam2(mop_idx,:)=ones(1,size(m1a_tcam2,2))*NaN;
%relatives setting contra homotopic to NaN
v1a_rcam2=v1a_tcam2./nansum(v1a_tcam2);
s1a_rcam2=s1a_tcam2./nansum(s1a_tcam2);
m1a_rcam2=m1a_tcam2./nansum(m1a_tcam2);

% relatives numbers per animal for ipsi and contra NOT setting homotopic
% contalteral to NAN
v1a_rcam=squeeze(nansum(c_v1aam(:,:,:))./nansum(nansum(c_v1aam(:,:,:))));
v1a_riam=squeeze(nansum(i_v1aam(:,:,:))./nansum(nansum(i_v1aam(:,:,:))));
s1a_rcam=squeeze(nansum(c_s1aam(:,:,:))./nansum(nansum(c_s1aam(:,:,:))));
s1a_riam=squeeze(nansum(i_s1aam(:,:,:))./nansum(nansum(i_s1aam(:,:,:))));
m1a_rcam=squeeze(nansum(c_m1aam(:,:,:))./nansum(nansum(c_m1aam(:,:,:))));
m1a_riam=squeeze(nansum(i_m1aam(:,:,:))./nansum(nansum(i_m1aam(:,:,:))));

%IPSI: 
%replace the zeros with NaNs for the injection areas IPSI
%ipsi absolute
v1a_tiam(visp_idx,:)=ones(1,size(v1a_tiam,2))*NaN;
s1a_tiam(ssp_idx,:)=ones(1,size(s1a_tiam,2))*NaN;
m1a_tiam(mop_idx,:)=ones(1,size(m1a_tiam,2))*NaN;
%ispi relative
v1a_riam(visp_idx,:)=ones(1,size(v1a_riam,2))*NaN;
s1a_riam(ssp_idx,:)=ones(1,size(s1a_riam,2))*NaN;
m1a_riam(mop_idx,:)=ones(1,size(m1a_riam,2))*NaN;
%layers abs
for i=1:size(i_v1aam,3)
    i_v1aam(:,visp_idx,i)=ones(1,size(i_v1aam,1))*NaN;
end
for i=1:size(i_s1aam,3)
    i_s1aam(:,ssp_idx,i)=ones(1,size(i_s1aam,1))*NaN;
end
for i=1:size(i_m1aam,3)
    i_m1aam(:,mop_idx,i)=ones(1,size(i_m1aam,1))*NaN;
end
%layers relative
for i=1:size(i_v1aafm,3)
    i_v1aafm(:,visp_idx,i)=ones(1,size(i_v1aafm,1))*NaN;
end
for i=1:size(i_s1aam,3)
    i_s1aafm(:,ssp_idx,i)=ones(1,size(i_s1aafm,1))*NaN;
end
for i=1:size(i_m1aafm,3)
    i_m1aafm(:,mop_idx,i)=ones(1,size(i_m1aafm,1))*NaN;
end

%% add total numbers in a strcuture/csv -> for heat maps in flat map
cd(str)
ipsi_contra_cortex_homo_included =table(cortex_names(:,2),v1a_tm(:,1),v1a_tm(:,2),s1a_tm(:,1),s1a_tm(:,2),m1a_tm(:,1),m1a_tm(:,2),...
    'variablenames',{'Area','V1ipsi','V1contra','S1ipsi','S1contra','M1ipsi','M1contra'});
writetable(ipsi_contra_cortex_homo_included,'ipsi_contra_cortex_homo_included.csv')
%% add total numbers in a strcuture/csv -> for heat maps in flat map
%THIS IS USED FOR Figure 1E heatmaps 
v1a_tm2=[];s1a_tm2=[];m1a_tm2=[];
v1a_tm2=v1a_tm;s1a_tm2=s1a_tm;m1a_tm2=m1a_tm;
v1a_tm2(visp_idx,:)=0;s1a_tm2(ssp_idx,:)=0;m1a_tm2(mop_idx,:)=0;
ipsi_contra_cortex_homo_excluded =table(cortex_names(:,2),v1a_tm2(:,1),v1a_tm2(:,2),s1a_tm2(:,1),s1a_tm2(:,2),m1a_tm2(:,1),m1a_tm2(:,2),...
    'variablenames',{'Area','V1ipsi','V1contra','S1ipsi','S1contra','M1ipsi','M1contra'});
cd(str)
writetable(ipsi_contra_cortex_homo_excluded,'ipsi_contra_cortex_homo_excluded.csv')

% %% sort contra 
% contra_normalized_visp=v1a_tm(:,2)/sum(v1a_tm(:,2));
% ipsi_normalized_visp=v1a_tm(:,1)/sum(v1a_tm(:,1));
% [kk ll]=sort(contra_normalized_visp,'descend')
% cortex_names(ll,2);
% fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 200, 800]);set(gcf,'color','w');
% imagesc(kk);
% hh=cmap(v1_color,100,10,5);colorbar;
% colormap(flipud(hh));box off;yticks([1:1:45]);yticklabels(cortex_names(ll,2));
% set(gca,'TickDir','out');set(gca,'FontSize',11);xticklabels({[]});h = gca;h.XAxis.Visible = 'off';title('contra')


%% %%%%% FIGURES START HERE%%%%

%% Figure 1: Widespread and bilaterally symmetrical cortical projections to VISp, SSp-bfd and MOp
%PANEL A top: illustrator
%PANEL A bottom: code from Manueel Teichert 
%PANEL B: cellfinder software and illustrator  
%PANEL E: Use above created csv ipsi_contra_cortex_homo_excluded.csv. The
%final version was created manually using GraphPad Prism 10, Please contact
%Manuel Teichert directly
%% PANEL C: Plot all cell numbers ipsi contra colour-coded per injection type
%all super imposed
clc;
temp_color=[v1_color ;s1_color; m1_color];module_names={'VISp','SSp-bfd','MOp'};
all_dati={};all_dati={squeeze(nansum(nansum(i_v1aam)))./(squeeze(nansum(nansum(i_v1aam)))+squeeze(nansum(nansum(c_v1aam))))...
     squeeze(nansum(nansum(i_s1aam)))./(squeeze(nansum(nansum(i_s1aam)))+squeeze(nansum(nansum(c_s1aam))))...
     squeeze(nansum(nansum(i_m1aam)))./(squeeze(nansum(nansum(i_m1aam)))+squeeze(nansum(nansum(c_m1aam))))};
all_datc={};all_datc={squeeze(nansum(nansum(c_v1aam)))./(squeeze(nansum(nansum(i_v1aam)))+squeeze(nansum(nansum(c_v1aam))))...
     squeeze(nansum(nansum(c_s1aam)))./(squeeze(nansum(nansum(i_s1aam)))+squeeze(nansum(nansum(c_s1aam))))...
     squeeze(nansum(nansum(c_m1aam)))./(squeeze(nansum(nansum(i_m1aam)))+squeeze(nansum(nansum(c_m1aam))))};
%plot
fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 200, 250]);set(gcf,'color','w');
for j=1:3
dat=[];dat=[all_dati{:,j} all_datc{:,j}];
%lines between paired data points
for i=1:length(data)
     pl=plot([1,2],[dat(:,1),dat(:,2)],'color',[0.5 0.5 0.5]);    
end
%single data points
hold on;pS=plotSpread([dat(:,1),dat(:,2)],'categoryIdx',[ones(1,length(dat(:,1)))' ones(1,length(dat(:,2)))'*2],...
    'categoryMarkers',{'o','o'},'categoryColors',{temp_color(j,:), temp_color(j,:)});hold on;
%mean +- SEM
hold on;er1=errorbar([0.2+j*0.14],nanmean(dat(:,1)),nanstd(dat(:,1),[],1)/sqrt(length(dat(:,1))),'ok','MarkerFaceColor',temp_color(j,:),'Markersize',8);
hold on;er2=errorbar([2.2+j*0.14],nanmean(dat(:,2)),nanstd(dat(:,2),[],1)/sqrt(length(dat(:,1))),'ok','MarkerFaceColor',temp_color(j,:),'Markersize',8);
%legend colour
text(2.2,1.1-0.07*j,module_names{j},'Color',temp_color(j,:));
%stats
%test for normality
[hnorm_i(:,j) pnorm_i(:,j)]=adtest(dat(:,1));
[hnorm_c(:,j) pnorm_c(:,j)]=adtest(dat(:,2));
%[u p1]=ttest(dat(:,1),dat(:,2));
[p1 u]=signrank(dat(:,1),dat(:,2))
end
xticklabels({'ipsi','contra'});ylabel('Fraction');hold on;title([]);set(gca,'FontSize',11);set(gca,'TickDir','out');box off;
hold on;text(1.25,1,['*'],'FontSize',18);offsetAxes;h = gca;h.XAxis.Visible = 'off';
t1=text(0.5,-0.1,'ipsi','FontSize',11,'Color',[0.7 0.7 0.7]);set(t1,'Rotation',45);t1=text(1.5,-0.1,'contra','FontSize',11,'Color',[0.3 0.3 0.3]);set(t1,'Rotation',45); ax=gca;ax.LineWidth=1;
% save
cd(save_folder);saveas(gcf, 'ipsi_contra_fraction.pdf');
%% PANEL D: Bar plots displaying number of projection areas onto VISp, SSp-bfd and MOp
%here contra homotopic needs to be removed so use v1a_tcam2 etc
clc;
temp_color=[v1_color ;s1_color; m1_color];
%ipsi_total numbers area
p1=[];p1=sum(v1a_tiam>0);p2=[];p2=sum(s1a_tiam>0);p3=[];p3=sum(m1a_tiam>0);
%contra_total numbers area
p4=[];p4=sum(v1a_tcam2>0);p5=[];p5=sum(s1a_tcam2>0);p6=[];p6=sum(m1a_tcam2>0);
temp_p= {p1' p2' p3'};temp_p2= {p4' p5' p6'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOT
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 300, 270]);set(gcf,'color','w');
for i=1:3
    b2=bar([i],[nanmean(temp_p{:,i}) nanmean(temp_p2{:,i})],0.7);
     hold on;b2(1).FaceColor=[0.8 0.8 0.8];b2(1).EdgeColor=temp_color(i,:);b2(1).LineWidth=1.4;set(b2,'ShowBaseLine','off')
     hold on;b2(2).FaceColor=[0.3 0.3 0.3];b2(2).EdgeColor=temp_color(i,:);b2(2).LineWidth=1.4;set(b2,'ShowBaseLine','off')
     hold on;errorbar([i-0.15],[nanmean(temp_p{:,i})],[nanstd(temp_p{:,i})/sqrt(length(temp_p{:,i}))]...
      , 'LineStyle', 'none', ... 
          'Color', 'k', 'LineWidth', 1.2,'CapSize',0);hold on;
      hold on;errorbar([i+0.15],[nanmean(temp_p2{:,i})],[nanstd(temp_p2{:,i})/sqrt(length(temp_p2{:,i}))]...
      , 'LineStyle', 'none', ... 
          'Color', 'k', 'LineWidth', 1.2,'CapSize',0);hold on;
      disp(['the range ipsi is ' num2str(min(temp_p{:,i})) ' - ' num2str(max(temp_p{:,i}))])
      disp(['the range contra is ' num2str(min(temp_p2{:,i})) ' - ' num2str(max(temp_p2{:,i}))])
end
t1=text(0.85,18,[num2str(min(p1)) ' - ' num2str(max(p1))],'FontSize',10,'Color','k');set(t1,'Rotation',90);
t1=text(1.15,18,[num2str(min(p4)) ' - ' num2str(max(p4))],'FontSize',10,'Color','w');set(t1,'Rotation',90);t1=text(0.82,45,'ipsi','FontSize',11);set(t1,'Rotation',90);
t1=text(1.15,45,'contra','FontSize',11);set(t1,'Rotation',90);t1=text(1.85,18,[num2str(min(p2)) ' - ' num2str(max(p2))],'FontSize',10,'Color','k');set(t1,'Rotation',90);
t1=text(2.15,18,[num2str(min(p5)) ' - ' num2str(max(p5))],'FontSize',10,'Color','w');set(t1,'Rotation',90);t1=text(1.82,45,'ipsi','FontSize',11);set(t1,'Rotation',90);
t1=text(2.15,45,'contra','FontSize',11);set(t1,'Rotation',90);t1=text(2.85,9,[num2str(min(p3)) ' - ' num2str(max(p3))],'FontSize',10,'Color','k');set(t1,'Rotation',90);
t1=text(3.15,9,[num2str(min(p6)) ' - ' num2str(max(p6))],'FontSize',10,'Color','w');set(t1,'Rotation',90);
t1=text(2.82,45,'ipsi','FontSize',11);set(t1,'Rotation',90);t1=text(3.15,45,'contra','FontSize',11);set(t1,'Rotation',90);set(gca,'FontSize',11);set(gca,'TickDir','out');box off;xtickangle(45);ylim([-1 52]);axis off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TEST
%Normality ipsi
hnorm_i=[];pnorm_i=[];
x1=[];x2=[];x3=[];x1 = temp_p{:,1};x2 = temp_p{:,2};x3 = temp_p{:,3};
[hnorm_i(:,1) pnorm_i(:,1)]=adtest(x1);
[hnorm_i(:,2) pnorm_i(:,2)]=adtest(x2);
[hnorm_i(:,3) pnorm_i(:,3)]=adtest(x3);
%Normality contra
hnorm_c=[];pnorm_c=[];
x1=[];x2=[];x3=[];x1 = temp_p2{:,1};x2 = temp_p2{:,2};x3 = temp_p2{:,3};
[hnorm_c(:,1) pnorm_c(:,1)]=adtest(x1);
[hnorm_c(:,2) pnorm_c(:,2)]=adtest(x2);
[hnorm_c(:,3) pnorm_c(:,3)]=adtest(x3);
% 
% dat_test = [x1' x2' x3']; %// Create row vector with your data
% group = {'G1','G1','G1','G1','G1','G1','G2','G2','G2','G2','G2','G2','G3','G3','G3','G3','G3','G3'}; %// set the groups according to the data above
% %[p,tbl,stats]  = anova1(dat_test, group) %// Use the 'off' option to prevent the table/box plot from showing up.
% [p,tbl,stats] = kruskalwallis(dat_test, group)
% presults = multcompare(stats)
%dat_test = [x1' x2' x3']; %// Create row vector with your data
% %[p,tbl,stats]  = anova1(dat_test, group)
% [p,tbl,stats] = kruskalwallis(dat_test, group)
% presults = multcompare(stats)

%Sensory vs motor (n=12 vs n=6)
%ipsi only
disp('ipsi sensory vs motor');
%[u p11]=ttest2([temp_p{:,1} ;temp_p{:,2}],temp_p{:,3})
[p11 u]=ranksum([temp_p{:,1} ;temp_p{:,2}],temp_p{:,3})
nanmean([temp_p{:,1} ;temp_p{:,2}])
nanstd([temp_p{:,1} ;temp_p{:,2}]/sqrt(length([temp_p{:,1} ;temp_p{:,2}])))
nanmean(temp_p{:,3})
nanstd([temp_p{:,3}]/sqrt(length([temp_p{:,3}])))
% %contra
% disp('contra sensory vs motor');
% %[u p11]=ttest2([temp_p2{:,1} ;temp_p2{:,2}],temp_p2{:,3})
% [p11 u]=ranksum([temp_p2{:,1} ;temp_p2{:,2}],temp_p2{:,3})
% nanmean([temp_p2{:,1} ;temp_p2{:,2}])
% nanstd([temp_p2{:,1} ;temp_p2{:,2}]/sqrt(length([temp_p2{:,1} ;temp_p2{:,2}])))
% nanmean(temp_p2{:,3})
% nanstd([temp_p2{:,3}]/sqrt(length([temp_p2{:,3}])))

%IPSI vs CONTRA for resepctive target area
disp('VISp test ipsi vs contra');
%[u p11]=ttest(temp_p{:,1},temp_p2{:,1})
[p11 u]=signrank(temp_p{:,1},temp_p2{:,1})
nanmean(temp_p{:,1})
nanstd([temp_p{:,1}]/sqrt(length([temp_p{:,1}])))
nanmean(temp_p2{:,1})
nanstd([temp_p2{:,1}]/sqrt(length([temp_p2{:,1}])))
disp('SSp test ipsi vs contra');
%[u p11]=ttest(temp_p{:,2},temp_p2{:,2})
[p11 u]=signrank(temp_p{:,2},temp_p2{:,2})
nanmean(temp_p{:,2})
nanstd([temp_p{:,2}]/sqrt(length([temp_p{:,2}])))
nanmean(temp_p2{:,2})
nanstd([temp_p2{:,2}]/sqrt(length([temp_p2{:,2}])))
disp('MOp test ipsi vs contra');
%[u p11]=ttest(temp_p{:,3},temp_p2{:,3})
[p11 u]=signrank(temp_p{:,3},temp_p2{:,3})
nanmean(temp_p{:,3})
nanstd([temp_p{:,3}]/sqrt(length([temp_p{:,3}])))
nanmean(temp_p2{:,3})
nanstd([temp_p2{:,3}]/sqrt(length([temp_p2{:,3}])))
% save
cd(save_folder);saveas(gcf, 'ipsi_contra_are_numbercounts.pdf');

%% PANEL E bottom right: Contralateral fraction cells in homotopic vs fraction in heterotpic areas for VISp, SSp-bfd and MOp
clc;
temp_color=[v1_color ;s1_color; m1_color];
%contra fraction homo
p1=[];p1=v1a_rcam(visp_idx,:);p2=[];p2=s1a_rcam(ssp_idx,:);p3=[];p3=m1a_rcam(mop_idx,:);
%contra fraction hetero
p4=[];p4=sum(v1a_rcam([1:visp_idx-1 visp_idx+1:45] ,:));p5=[];p5=sum(s1a_rcam([1:ssp_idx-1 ssp_idx+1:45] ,:));p6=[];p6=sum(m1a_rcam([1:mop_idx-1 mop_idx+1:45] ,:));
temp_p= [p1' p2' p3'];temp_p2= [p4' p5' p6'];
%plot
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 300, 270]);set(gcf,'color','w');
for i=1:3
    b1=bar([i*2],[nanmean(temp_p2(:,i)) nanmean(temp_p(:,i))],0.5,'stacked');hold on;b1(2).FaceColor=[1 1 1];b1(1).FaceColor=temp_color(i,:);
    b1(1).FaceAlpha=[0.6];b1(2).EdgeColor=temp_color(i,:);b1(1).LineWidth=1;b1(2).LineWidth=1;b1(1).EdgeColor=temp_color(i,:);
     hold on;errorbar([i*2],[nanmean(temp_p2(:,i)) nanmean(temp_p(:,i))+nanmean(temp_p2(:,i))],[nanstd(temp_p2(:,i))/sqrt(length(temp_p2(:,i))) nanstd(temp_p(:,i))/sqrt(length(temp_p(:,i)))]...
     , 'LineStyle', 'none', ... 
         'Color', 'k', 'LineWidth', 1.2,'CapSize',0);hold on;
     set(b1,'ShowBaseLine','off')
end
xticks([2:2:6]);hold on;box off;xticklabels({[]});%ylabel('Nr of areas');
set(gca,'FontSize',11);set(gca,'TickDir','out');box off;xtickangle(45);axis off
text(1.8,0.9,[num2str(round(nanmean(p1*100)),2) '%'],'FontSize',8,'Color','k');text(1.8,0.45,[num2str(round(nanmean(p4*100)),2) '%'],'FontSize',8,'Color','k');
text(3.8,0.9,[num2str(round(nanmean(p2*100)),2) '%'],'FontSize',8,'Color','k');text(3.8,0.45,[num2str(round(nanmean(p5*100)),2) '%'],'FontSize',8,'Color','k');
text(5.8,0.75,[num2str(round(nanmean(p3*100)),2) '%'],'FontSize',8,'Color','k');text(5.8,0.27,[num2str(round(nanmean(p6*100)),2) '%'],'FontSize',8,'Color','k');
 t1=text(1.3,0.23,'Other areas','FontSize',11);set(t1,'Rotation',90);t1=text(1.3,0.85,'#','FontSize',11);set(t1,'Rotation',90);
 t1=text(3.3,0.23,'Other areas','FontSize',11);set(t1,'Rotation',90);t1=text(3.3,0.85,'#','FontSize',11);set(t1,'Rotation',90);
 t1=text(5.3,0.07,'Other areas','FontSize',11);set(t1,'Rotation',90);t1=text(5.3,0.65,'#','FontSize',11);set(t1,'Rotation',90);
% save
cd(save_folder);saveas(gcf, 'homo_hetero_fraction.pdf');
%% PANEL F: Correlation between hemisphere areas, R2 square of fraction per area avergae across injection areas
clc;
%Calculate R2 and slope of correlation 
%V1
for i=1:size(v1a_riam,2)
idx = find(~isnan(v1a_riam(:,i)));
mdl = fitlm(v1a_riam(idx,i),v1a_rcam2(idx,i));
coefs = polyfit(v1a_riam(idx,i)', v1a_rcam2(idx,i)', 1);
slope_v1(i)=coefs(1);
r2_sqv1(i)=mdl.Rsquared.Ordinary; 
end
%S1
for i=1:size(s1a_riam,2)
idx = find(~isnan(s1a_riam(:,i)));
mdl = fitlm(s1a_riam(idx,i),s1a_rcam2(idx,i));
coefs = polyfit(s1a_riam(idx,i)', s1a_rcam2(idx,i)', 1);
slope_s1(i)=coefs(1);
r2_sqs1(i)=mdl.Rsquared.Ordinary;
end
%M1
for i=1:size(m1a_riam,2)
idx = find(~isnan(m1a_riam(:,i)));
mdl = fitlm(m1a_riam(idx,i),m1a_rcam2(idx,i));
coefs = polyfit(m1a_riam(idx,i)', m1a_rcam2(idx,i)', 1);
slope_m1(i)=coefs(1);
r2_sqm1(i)=mdl.Rsquared.Ordinary;
end

%R2
tm1=[];tm2=[];tm3=[];tm1=r2_sqv1;tm2=r2_sqs1;tm3=r2_sqm1;
%%%%%%%%%%%%%%%
%plot
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 170, 250]);set(gcf,'color','w');tiledlayout("horizontal")
b1=bar(1,nanmean(tm1),0.7);hold on;b1.FaceColor=v1_color;b2=bar(2,nanmean(tm2));hold on;b2.FaceColor=s1_color;
b3=bar(3,nanmean(tm3),0.7);b3.FaceColor=m1_color;
set(b1,'ShowBaseLine','off');set(b2,'ShowBaseLine','off');set(b3,'ShowBaseLine','off');
r=1;rng(1);r1 = r-0.01 + (0.2)*rand(length(tm1),1);
sc1=scatter(r1,tm1,15,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
 r=2;rng(1);r1 = r-0.01 + (0.2)*rand(length(tm2),1);
 sc1=scatter(r1,tm2,15,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
 r=3;rng(1);r1 = r-0.01 + (0.2)*rand(length(tm3),1);
 sc1=scatter(r1,tm3,15,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];xticks([1:3]);hold on;
hold on;errorbar([1 2 3],[nanmean(tm1) nanmean(tm2) nanmean(tm3)],[nanstd(tm1)/sqrt(size(tm1,2)) nanstd(tm2)/sqrt(size(tm2,2)) nanstd(tm3)/sqrt(size(tm3,2))]...
    , 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1.2,'CapSize',0);hold on;xticklabels({'VISp','SSpbf','MOp'});
ylabel({'Fractional count ipsi vs contra' ,'Correlation R2'});set(gca,'FontSize',11);set(gca,'TickDir','out');box off;xtickangle(45);
 h = gca;h.XAxis.Visible = 'off';offsetAxes;
t1=text(0.4,-0.2,'VISp','FontSize',11);set(t1,'Rotation',45);t1=text(0.8,-0.3,'SSp-bfd','FontSize',11);set(t1,'Rotation',45);t1=text(2.4,-0.2,'MOp','FontSize',11);set(t1,'Rotation',45);
disp('VISp R2')
nanmean(tm1)
nanstd(tm1)/sqrt(size(tm1,2))
 disp('SSp R2')
nanmean(tm2)
nanstd(tm2)/sqrt(size(tm2,2))
 disp('MOp R2')
nanmean(tm3)
nanstd(tm3)/sqrt(size(tm3,2))
% save
cd(save_folder);saveas(gcf, 'r2_correlation_hemi.pdf');
%% PANEL G: Bilateral symmetry test and logarthmic plot
clc;
%creat HEMI Index that will be tested against 0
% v1_hemi=(v1a_rcam2-v1a_riam)./(v1a_rcam2+v1a_riam);
% s1_hemi=(s1a_rcam2-s1a_riam)./(s1a_rcam2+s1a_riam);
% m1_hemi=(m1a_rcam2-m1a_riam)./(m1a_rcam2+m1a_riam);
%use different 
v1_hemi=(v1a_rcam2-v1a_riam);
s1_hemi=(s1a_rcam2-s1a_riam);
m1_hemi=(m1a_rcam2-m1a_riam);
% Plot 
temp_color=[v1_color ;s1_color; m1_color];
hemi_all={v1_hemi s1_hemi m1_hemi};

temp_data_i={v1a_riam s1a_riam m1a_riam};
temp_data_c={v1a_rcam2 s1a_rcam2 m1a_rcam2};
panel_title={'VISp','SSp-bfd','MOp'};crt_all1={};
%%%%%%%%%%%%%%%%%%%
%Plot
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 1000, 500]);set(gcf,'color','w');tiledlayout("horizontal")
for i=1:3
nexttile
%CALL FUNCTION anatomy_testarea
[p_bin p_bin_c] = anatomy_testarea(hemi_all{i},temp_data_i{i},temp_data_c{i});

rm1=[];rm2=[];
rm1=temp_data_i{i};rm2=temp_data_c{i};
dir_diff=[];dir_diff=(nanmean(rm2,2)-nanmean(rm1,2))';
disp(['bilateral' num2str(i)])
sum(p_bin>0 | nanmean(rm2,2)'==0)
%1)
errx_l=[];errx_l=[nanstd(rm1(find(p_bin==0 & p_bin_c==0),:),[],2)]/sqrt(6);
erry_l=[];erry_l=[nanstd(rm2(find(p_bin==0 & p_bin_c==0),:),[],2)]/sqrt(6);
temp_datax=[];temp_datay=[];diferrx=[];diferry=[];
temp_datax=nanmean(rm1(find(p_bin==0 & p_bin_c==0),:),2);
temp_datay=nanmean(rm2(find(p_bin==0 & p_bin_c==0),:),2);
diferrx=temp_datax-errx_l;
diferry=temp_datay-erry_l;
%call util function ploterr
h=ploterr(temp_datax,temp_datay,errx_l,erry_l,'logx','logy','hhxy',0);
h(1).LineStyle='none';h(2).Color=[0.9 0.9 0.9];h(3).Color=[0.9 0.9 0.9];hold on;
crti=find(diferrx<0 | diferry<0);
crt_all1{i}=crti;
if isempty(crti)==0
for k=1:length(crti)
l1=line([temp_datax(crti(k)) temp_datax(crti(k))+errx_l(crti(k))],[temp_datay(crti(k)) temp_datay(crti(k))],'Color',[0.9 0.9 0.9]);
l1=line([temp_datax(crti(k)) temp_datax(crti(k))],[temp_datay(crti(k)) temp_datay(crti(k))+erry_l(crti(k))],'Color',[0.9 0.9 0.9]);
end
end
% if i==1
%     l1=line([0.0001 temp_datax(crti(1))],[temp_datay(crti(1)) temp_datay(crti(1))],'Color',[0.9 0.9 0.9]);
% l1=line([temp_datax(crti(1)) temp_datax(crti(1))],[0.0001 temp_datay(crti(1))],'Color',[0.9 0.9 0.9]);
% end
%2)
errx_l=[];errx_l=[nanstd(rm1(find(p_bin>0 & dir_diff<0),:),[],2)]/sqrt(6);
erry_l=[];erry_l=[nanstd(rm2(find(p_bin>0 & dir_diff<0),:),[],2)]/sqrt(6);
temp_datax=[];temp_datay=[];diferrx=[];diferry=[];
temp_datax=nanmean(rm1(find(p_bin>0 & dir_diff<0),:),2);
temp_datay=nanmean(rm2(find(p_bin>0 & dir_diff<0),:),2);
diferrx=temp_datax-errx_l;
diferry=temp_datay-erry_l;
h=ploterr(temp_datax,temp_datay,errx_l,erry_l,'logx','logy','hhxy',0);
h(1).LineStyle='none';h(2).Color=[0.9 0.9 0.9];h(3).Color=[0.9 0.9 0.9];hold on;
crti=find(diferrx<0 | diferry<0);
if isempty(crti)==0
for k=1:length(crti)
l1=line([temp_datax(crti(k)) temp_datax(crti(k))+errx_l(crti(k))],[temp_datay(crti(k)) temp_datay(crti(k))],'Color',[0.9 0.9 0.9]);
l1=line([temp_datax(crti(k)) temp_datax(crti(k))],[temp_datay(crti(k)) temp_datay(crti(k))+erry_l(crti(k))],'Color',[0.9 0.9 0.9]);
end
end



%3) 
errx_l=[];errx_l=[nanstd(rm1(find(p_bin>0 & dir_diff>0),:),[],2)]/sqrt(6);
erry_l=[];erry_l=[nanstd(rm2(find(p_bin>0 & dir_diff>0),:),[],2)]/sqrt(6);
temp_datax=[];temp_datay=[];diferrx=[];diferry=[];
temp_datax=nanmean(rm1(find(p_bin>0 & dir_diff>0),:),2);
temp_datay=nanmean(rm2(find(p_bin>0 & dir_diff>0),:),2);
diferrx=temp_datax-errx_l;
diferry=temp_datay-erry_l;
h=ploterr(temp_datax,temp_datay,errx_l,erry_l,'logx','logy','hhxy',0);
h(1).LineStyle='none';h(2).Color=[0.9 0.9 0.9];h(3).Color=[0.9 0.9 0.9];hold on;
crti=find(diferrx<0 | diferry<0);
if isempty(crti)==0
for k=1:length(crti)
l1=line([temp_datax(crti(k)) temp_datax(crti(k))+errx_l(crti(k))],[temp_datay(crti(k)) temp_datay(crti(k))],'Color',[0.9 0.9 0.9]);
l1=line([temp_datax(crti(k)) temp_datax(crti(k))],[temp_datay(crti(k)) temp_datay(crti(k))+erry_l(crti(k))],'Color',[0.9 0.9 0.9]);
end
end

loglog(nanmean(rm1(find(p_bin>0 & dir_diff<0),:),2),nanmean(rm2(find(p_bin>0 & dir_diff<0),:),2),'o','MarkerFaceColor',[0.8 0.8 0.8],'MarkerEdgeColor','k');hold on;
loglog(nanmean(rm1(find(p_bin>0 & dir_diff>0),:),2),nanmean(rm2(find(p_bin>0 & dir_diff>0),:),2),'o','MarkerFaceColor',[0.3 0.3 0.3],'MarkerEdgeColor','k');hold on;
loglog(nanmean(rm1(find(p_bin==0 & p_bin_c==0),:),2),nanmean(rm2(find(p_bin==0 & p_bin_c==0),:),2),'o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor','k');hold on;axis square



xlabel('Ipsi fraction (log)');ylabel('Contra fraction (log)');set(gca,'FontSize',11);set(gca,'TickDir','out');box off;title(panel_title{i},'Color',temp_color(i,:),'FontWeight','normal');
if i==1
xlim([0.0001 0.2]);
ylim([0.0001 0.2]);
set(gca,'XLim',[0.0001 0.2],'XTick',10.^(-4:-1), ...
        'YLim',[0.0001 0.2],'YTick',10.^(-4:-1))
hold on;rf=refline(1,0);rf.Color='k';hold on;rf.LineStyle=':';
elseif i==2
 xlim([0.001 0.25]);
 ylim([0.0001 0.25]);
 set(gca,'XLim',[0.0001 0.25],'XTick',10.^(-4:-1), ...
        'YLim',[0.0001 0.25],'YTick',10.^(-4:-1))
 hold on;rf=refline(1,0);rf.Color='k';hold on;rf.LineStyle=':';
else
xlim([0.0001 0.6]);
ylim([0.0001 0.6]);
 set(gca,'XLim',[0.0001 0.6],'XTick',10.^(-4:-1), ...
        'YLim',[0.0001 0.6],'YTick',10.^(-4:-1))
 hold on;rf=refline(1,0);rf.Color='k';hold on;rf.LineStyle=':';
end

set(gca, 'XScale','log', 'YScale','log');
axis square;
end

cd(save_folder);saveas(gcf, 'correlation_log_fraction.pdf');


%% Figure 2: Modular organization of the intra and interhemispheric cortical projection onto VISp, SSp-bfd and MOp. 
%Fig 2A: Illustrator
%% PANEL B: ipsi vs contra shown in modules, sorted by largest module
clc;
%dat_all={};dat_all={v1a_riam v1a_rcam; s1a_riam s1a_rcam; m1a_riam m1a_rcam};
dat_all={};dat_all={v1a_riam v1a_rcam2; s1a_riam s1a_rcam2; m1a_riam m1a_rcam2};
frontal_idx=[1:8];lateral_idx=[9:16 44:45];somamo_idx=[17:26];visual_idx=[27:36];medial_idx=[37:39];aud_idx=[40:43];
all_ipsim=[];all_contram=[];temp_color=[v1_color ;s1_color; m1_color];panel_title={'VISp','SSp-bfd','MOp'};name_sub={'modules_VISp.pdf', 'modules_SSp.pdf' ,'modules_MOp.pdf'};
order_modules={};order_modules_names={};
order_modules={visual_idx lateral_idx medial_idx aud_idx frontal_idx somamo_idx;...
    somamo_idx visual_idx lateral_idx frontal_idx aud_idx medial_idx;...
    somamo_idx lateral_idx frontal_idx medial_idx visual_idx aud_idx};
order_modules_names={'Vis','Lat','Med','Aud','Prefron','SoMo';...
    'SoMo','Vis','Lat','Prefron','Aud','Med';...
    'SoMo','Lat','Prefron','Med','Vis','Aud'};
%%%%%%%%%%%Plot
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [200, 500, 800, 240]);set(gcf,'color','w');t=tiledlayout("horizontal",'TileSpacing','Compact');
m=[];
for m=1:3
    figure(fig7);nexttile
rm1=[];rm2=[];rm1=dat_all{m,1};rm2=dat_all{m,2};
p1=[];p1=nansum(rm1(order_modules{m,1},:));p2=[];p2=nansum(rm1(order_modules{m,2},:));p3=[];p3=nansum(rm1(order_modules{m,3},:));
p4=[];p4=nansum(rm1(order_modules{m,4},:));p5=[];p5=nansum(rm1(order_modules{m,5},:));p6=[];p6=nansum(rm1(order_modules{m,6},:));
p1_1=[];p1_1=nansum(rm2(order_modules{m,1},:));p2_1=[];p2_1=nansum(rm2(order_modules{m,2},:));p3_1=[];p3_1=nansum(rm2(order_modules{m,3},:));
p4_1=[];p4_1=nansum(rm2(order_modules{m,4},:));p5_1=[];p5_1=nansum(rm2(order_modules{m,5},:));p6_1=[];p6_1=nansum(rm2(order_modules{m,6},:));
temp_p=[];temp_p= [p1' p2' p3' p4' p5' p6'];
temp_p2=[];temp_p2= [p1_1' p2_1' p3_1' p4_1' p5_1' p6_1'];
all_ipsim(:,:,m)=temp_p;
all_contram(:,:,m)=temp_p2;
% [uuj uui]=sort(temp_p)
% [uuj uuc]=sort(temp_p2)
% [pw ,ee, stats]=signrank(uui(:),uuc(:))
title(panel_title{m},'FontWeight','normal','Color',temp_color(m,:));

    for i=1:6
        b1=bar(i,[nanmean(temp_p(:,i));nanmean(temp_p2(:,i))],1);hold on;b1(1).FaceColor=[0.8 0.8 0.8];b1(2).FaceColor=[0.3 0.3 0.3];set(b1,'ShowBaseLine','off');
        b1(1).EdgeColor=[1 1 1];b1(2).EdgeColor=[1 1 1];
        % disp(['stats ipsi vs contra' panel_title{m}])
        % [ufig2(i) p11]=ttest(temp_p(:,i),temp_p2(:,i))
         r=i;rng(i);r1 = r-0.01 + (0.2)*rand(length(temp_p(:,i)),1);
         r2=i;rng(i);r2 = r2-0.01 + (0.2)*rand(length(temp_p2(:,i)),1);
         sc1=scatter(r1-0.2,temp_p(:,i),6,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
         sc1=scatter(r2+0.15,temp_p2(:,i),6,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];

         hold on;h=errorbar([i-0.15],[nanmean(temp_p(:,i))],[nanstd(temp_p(:,i))/sqrt(length(temp_p(:,i)))], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
         hold on;h=errorbar([i+0.15],[nanmean(temp_p2(:,i))],[nanstd(temp_p2(:,i))/sqrt(length(temp_p2(:,i)))], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
         %disp(['module 1 vs others stats' panel_title{m}]);
         %[u p12]=ttest2([temp_p(:,1) ;temp_p2(:,1)],temp_p2(:,i))
       
        % hold on;errorbar([i],[nanmean(temp_p(:,i))],[nanstd(temp_p(:,i))/sqrt(length(temp_p(:,i)))]...
        % , 'LineStyle', 'none', ... 
        % 'Color', 'k', 'LineWidth', 1.2);hold on;   
    end

   
 if m==1
yticks([0:0.25:1]);
 xticks([1:6]);hold on;box off;xticklabels(order_modules_names(m,:));ylim([0 1]);set(gca,'FontSize',11);set(gca,'TickDir','out');box off;xtickangle(45);
 ylabel('Fraction of cells');
 text(3.5,0.9,'ipsi (normalized)','Color',[0.7 0.7 0.7]);
 text(3.5,0.8,'contra (normalized)','Color',[0.3 0.3 0.3]);
 else
 xticks([1:6]);xticklabels(order_modules_names(m,:));ylim([0 1]);set(gca,'FontSize',11);set(gca,'TickDir','out');box off;xtickangle(45);
 h = gca;h.YAxis.Visible = 'off';  
 end
   xlim([0.5 6.5]);
offsetAxes;

%save
%cd(save_folder);saveas(gcf, name_sub{m});
%figure;
 %disp(['modules ipsi+contra against each other ' panel_title{m}])
      % [p,tbl,stats] = anova1([[temp_p(:,1) ;temp_p2(:,1)] [temp_p(:,2) ;temp_p2(:,2)]...
      %     [temp_p(:,3) ;temp_p2(:,3)] [temp_p(:,4) ;temp_p2(:,4)] [temp_p(:,5) ;temp_p2(:,5)] [temp_p(:,6) ;temp_p2(:,6)]])
      % presults = multcompare(stats)
  disp(['modules ipsi against next higher ' panel_title{m}])  
  nanmean(temp_p(:,1))
  [nanstd(temp_p(:,1))/sqrt(length(temp_p(:,1)))]
nanmean(temp_p(:,2))
 [nanstd(temp_p(:,2))/sqrt(length(temp_p(:,2)))]
       [p u]=ttest2(temp_p(:,1),temp_p(:,2))
 disp(['modules contra against next higher ' panel_title{m}])  

 if m==1
       nanmean(temp_p2(:,1))
  [nanstd(temp_p2(:,1))/sqrt(length(temp_p2(:,1)))]
nanmean(temp_p2(:,2))
 [nanstd(temp_p2(:,2))/sqrt(length(temp_p2(:,2)))]
       [p u]=ttest2(temp_p2(:,1),temp_p2(:,2))

 elseif m==2
        nanmean(temp_p2(:,1))
  [nanstd(temp_p2(:,1))/sqrt(length(temp_p2(:,1)))]
nanmean(temp_p2(:,3))
 [nanstd(temp_p2(:,3))/sqrt(length(temp_p2(:,3)))]
       [p u]=ttest2(temp_p2(:,1),temp_p2(:,3))
 else m==3
       nanmean(temp_p2(:,1))
  [nanstd(temp_p2(:,1))/sqrt(length(temp_p2(:,1)))]
nanmean(temp_p2(:,2))
 [nanstd(temp_p2(:,2))/sqrt(length(temp_p2(:,2)))]
       [p u]=ttest2(temp_p2(:,1),temp_p2(:,2))

 end
disp(['lateral module ipsi vs contra ' panel_title{m}])  
 if m==1
     [p u]=ttest(temp_p(:,3),temp_p2(:,3),'Tail','left')
     nanmean(temp_p(:,3))
     [nanstd(temp_p(:,3))/sqrt(length(temp_p(:,3)))]
      nanmean(temp_p2(:,3))
     [nanstd(temp_p2(:,3))/sqrt(length(temp_p2(:,3)))]
 elseif m==2
     [p u]=ttest(temp_p(:,4),temp_p2(:,4),'Tail','left')
      nanmean(temp_p(:,4))
     [nanstd(temp_p(:,4))/sqrt(length(temp_p(:,4)))]
      nanmean(temp_p2(:,4))
     [nanstd(temp_p2(:,4))/sqrt(length(temp_p2(:,4)))]
 else m==3
     [p u]=ttest(temp_p(:,3),temp_p2(:,3),'Tail','left')
      nanmean(temp_p(:,3))
     [nanstd(temp_p(:,3))/sqrt(length(temp_p(:,3)))]
      nanmean(temp_p2(:,3))
     [nanstd(temp_p2(:,3))/sqrt(length(temp_p2(:,3)))]
 end


end
%% PANEL B: Rank order sorted and connect ipsi and contra areas for VISp, SSp-bfd and MOp
clc;
temp_color=linspecer(6);idx_modules={frontal_idx lateral_idx somamo_idx visual_idx medial_idx aud_idx};
module_names={'Pfron','Lat','SoMo','Vis','Med','Aud'};
frontal_idx=[1:8];lateral_idx=[9:16];somamo_idx=[17:26];visual_idx=[27:36];medial_idx=[37:39];aud_idx=[40:43];lateral_idx2=[44:45]
areas_idx=ones(45,1);
areas_modul=[];areas_modul=[areas_idx(frontal_idx)*1; areas_idx(lateral_idx)*2 ; areas_idx(somamo_idx)*3; ...
    areas_idx(visual_idx)*4; areas_idx(medial_idx)*5; areas_idx(aud_idx)*6 ; areas_idx(lateral_idx2)*2];
panel_tit={'VISp','SSp-bfd','MOp'};
dat_all={};dat_all={v1a_riam s1a_riam m1a_riam};
%excluding homotopic
dat_all2={};dat_all2={v1a_rcam2 s1a_rcam2 m1a_rcam2};

fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [100, 300, 1400, 700]);set(gcf,'color','w');tiledlayout("horizontal")
for j=1:3    
    temp_cortex1=[];temp_cortex1=nanmean(dat_all{j},2);
    temp_cortex2=[];temp_cortex2=nanmean(dat_all2{j},2);
 so_v1a_riam=[];idxi=[];cs1=[];so_v1a_rcam=[];idxc=[];tcon=[];
 %tcon=find(temp_cortex2>0 & temp_cortex1>0);
 tcon=[];
 tcon=find(~isnan(temp_cortex1));
[so_v1a_riam idxi]=sort(temp_cortex1(tcon),'descend');
[so_v1a_rcam idxc]=sort(temp_cortex2(tcon),'descend');
temp1=[];temp2=[];areas_coo1=[];areas_coo2=[];
temp1=areas_modul(tcon);temp2=areas_modul(tcon);temp3=cortex_names(tcon,2);temp4=cortex_names(tcon,2);
areas_coo1=temp1(idxi);areas_coo2=temp2(idxc);names_sel=temp3(idxi);names_sel2=temp4(idxc);
mom1=[length(tcon):-1:1];color_matrix=[];rank_data=[];
nexttile
for k=1:length(idxi)
    i2=[];c1=[];
    i2=find(idxi(k)==idxc);
    p1=plot([2 3],[mom1(k) mom1(i2)]);hold on;
    p3=plot([2],[mom1(k)],'o');hold on;
    p4=plot([3],[mom1(i2)],'o');
    if so_v1a_riam(k)==0 
        p1.LineStyle='none';
        p3.MarkerSize=0.01;
        values_dot(k)=0;
    else
    p3.MarkerSize=10+log(so_v1a_riam(k));
    values_dot(k)=(so_v1a_riam(k));
    end
      if so_v1a_rcam(i2)==0 
          p4.MarkerSize=0.01;
          values_dot(k)=0;
      else
           p4.MarkerSize=10+log(so_v1a_rcam(i2));
           values_dot(k)=(so_v1a_riam(k));
      end
   hold on;t1=text(1.5,mom1(k),names_sel{k});
      hold on;
     t2=text(3.2,mom1(i2),names_sel2{i2});
    rank_data(k,:)=[mom1(k) mom1(i2)];
    if areas_coo1(k)==1
    col1=temp_color(1,:)
    p1.MarkerFaceColor=col1;p1.MarkerEdgeColor=col1;p3.MarkerFaceColor=col1;p4.MarkerFaceColor=col1;p3.MarkerEdgeColor=col1;p4.MarkerEdgeColor=col1;
    p1.Color=col1;t1.Color=col1;t2.Color=col1;color_matrix(k,:)=col1;
   elseif areas_coo1(k)==2
    col1=temp_color(2,:);
    p1.MarkerFaceColor=col1;p1.MarkerEdgeColor=col1;p3.MarkerFaceColor=col1;p4.MarkerFaceColor=col1;p3.MarkerEdgeColor=col1;p4.MarkerEdgeColor=col1;
    p1.Color=col1;t1.Color=col1;t2.Color=col1;
    color_matrix(k,:)=col1;
        
    elseif areas_coo1(k)==3
          col1=temp_color(3,:);
    p1.MarkerFaceColor=col1;
    p1.MarkerEdgeColor=col1;
      p3.MarkerFaceColor=col1;
    p4.MarkerFaceColor=col1;
       p3.MarkerEdgeColor=col1;
    p4.MarkerEdgeColor=col1;
    p1.Color=col1;
    t1.Color=col1;
    t2.Color=col1;
    color_matrix(k,:)=col1;
        
    elseif areas_coo1(k)==4
         col1=temp_color(4,:);
    p1.MarkerFaceColor=col1;
    p1.MarkerEdgeColor=col1;
      p3.MarkerFaceColor=col1;
    p4.MarkerFaceColor=col1;
    p3.MarkerEdgeColor=col1;
    p4.MarkerEdgeColor=col1;
    p1.Color=col1;
    t1.Color=col1;
    t2.Color=col1;
    color_matrix(k,:)=col1;
    elseif areas_coo1(k)==5
          col1=temp_color(5,:);
    p1.MarkerFaceColor=col1;
    p1.MarkerEdgeColor=col1;
      p3.MarkerFaceColor=col1;
    p4.MarkerFaceColor=col1;
    p3.MarkerEdgeColor=col1;
    p4.MarkerEdgeColor=col1;
    p1.Color=col1;
    t1.Color=col1;
    t2.Color=col1;
    color_matrix(k,:)=col1;
          
    else areas_coo1(k)==6
          col1=temp_color(6,:);
    p1.MarkerFaceColor=col1;
    p1.MarkerEdgeColor=col1;
      p3.MarkerFaceColor=col1;
    p4.MarkerFaceColor=col1;
    p3.MarkerEdgeColor=col1;
    p4.MarkerEdgeColor=col1;
    p1.Color=col1;
    t1.Color=col1;
    t2.Color=col1;
    color_matrix(k,:)=col1;
       
    end
hold on;
box off;

xlim([1.5 3.5])
%offsetAxes;
 h = gca;h.XAxis.Visible = 'off';
 axis off
end
title(panel_tit{j})
diff_rank{j}=rank_data(:,1)-rank_data(:,2);
% figure(1)
% histogram(rank_data(:,1)-rank_data(:,2));hold on
text(1.8,47,'ipsi','Color',[0.7 0.7 0.7],'FontSize',11);
text(2.9,47,'contra','Color',[0.3 0.3 0.3],'FontSize',11);
set(gca,'FontSize',11)
if j==1
  
end
end
temp_scale=[32:2:40];
temp_thick=[2:2:10];
temp_thicktext={'0.0004','','','','0.4'};
nexttile

    for oo=1:length(temp_scale)
   pt= plot(1,temp_scale(oo),'o','Color','k');
   pt.MarkerSize=[temp_thick(oo)];hold on;
pt.MarkerFaceColor='k';
text(1.1,temp_scale(oo),temp_thicktext{oo},'FontSize',11);
    end
       text(0.9,43,'Fractional count','FontSize',11)
       text(0.9,27,'Modules','FontSize',11,'Color','k');
       text(0.9,25,'Prefrontal','FontSize',11,'Color',temp_color(1,:));
       text(0.9,23.5,'Lateral','FontSize',11,'Color',temp_color(2,:));
       text(0.9,22,'Somatomotor','FontSize',11,'Color',temp_color(3,:));
       text(0.9,20.5,'Visual','FontSize',11,'Color',temp_color(4,:));
       text(0.9,19,'Medial','FontSize',11,'Color',temp_color(5,:));
       text(0.9,17.5,'Auditory','FontSize',11,'Color',temp_color(6,:));  
    ylim([1 44])
    box off;axis off;
    set(gca,'FontSize',11)

cd(save_folder);saveas(gcf, 'rank_order_connectome.pdf');

%% FIGURE 3:  Cortical layer 6 is a major source of input to VISp, SSp-bfd and MOp
%Figure 3A: illustrator
%% Panel B: Plot all fraction laminar same plot differently sorted for layers (L2/3, L4, L5, L6a, L6b)
clc;
frontal_idx=[1:8];lateral_idx=[9:16 44:45];somamo_idx=[17:26];visual_idx=[27:36];medial_idx=[37:39];aud_idx=[40:43];
anim_incl={1:6 7:12 13:18};
labels_b={'L2/3','L4','L5','L6a' 'L6b'};temp_color=[v1_color ;s1_color; m1_color];panel_tit={'VISp','SSp-bfd','MOp'};

%%%Plot 
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 1200, 250]);set(gcf,'color','w');tiledlayout("horizontal");
for j=1:3
incl1=[];incl1=anim_incl{j};
%means
means_layers=[];err_layers=[];
means_layers=[nanmean(i_animal(2,incl1)) nanmean(i_animal(3,incl1)) nanmean(i_animal(4,incl1)) nanmean(i_animal(5,incl1)) nanmean(i_animal(6,incl1));...
    nanmean(c_animal(2,incl1)) nanmean(c_animal(3,incl1)) nanmean(c_animal(4,incl1)) nanmean(c_animal(5,incl1)) nanmean(c_animal(6,incl1)) ];
semd=sqrt(length(i_animal(:,incl1)));
%error
err_layers=[nanstd(i_animal(2,incl1))/semd nanstd(i_animal(3,incl1))/semd nanstd(i_animal(4,incl1))/semd nanstd(i_animal(5,incl1))/semd nanstd(i_animal(6,incl1))/semd; ...
    nanstd(c_animal(2,incl1))/semd nanstd(c_animal(3,incl1))/semd nanstd(c_animal(4,incl1))/semd nanstd(c_animal(5,incl1))/semd nanstd(c_animal(6,incl1))/semd];
nexttile
%bars
b1=bar(1,means_layers(1,:),0.7);
b1(1).FaceColor=[0.8 0.8 0.8];b1(2).FaceColor=[0.8 0.8 0.8];b1(3).FaceColor=[0.8 0.8 0.8];b1(4).FaceColor=[0.8 0.8 0.8];b1(5).FaceColor=[0.8 0.8 0.8];
b1(1).EdgeColor=[1 1 1];b1(2).EdgeColor=[1 1 1];b1(3).EdgeColor=[1 1 1];b1(4).EdgeColor=[1 1 1];b1(5).EdgeColor=[1 1 1];hold on;set(b1,'ShowBaseLine','off');
b2=bar(2,means_layers(2,:),0.7);
b2(1).FaceColor=[0.3 0.3 0.3];b2(2).FaceColor=[0.3 0.3 0.3];b2(3).FaceColor=[0.3 0.3 0.3];b2(4).FaceColor=[0.3 0.3 0.3];b2(5).FaceColor=[0.3 0.3 0.3];
b2(1).EdgeColor=[1 1 1];b2(2).EdgeColor=[1 1 1];b2(3).EdgeColor=[1 1 1];b2(4).EdgeColor=[1 1 1];b2(5).EdgeColor=[1 1 1];hold on;set(b2,'ShowBaseLine','off');
%indivdual animals
r=1;
rng(1);r1 = r-0.35 + (0.1)*rand(length(incl1),1);
rng(1);r2 = r -0.2+ (0.1)*rand(length(incl1),1);
rng(1);r3 = r-0.05 + (0.1)* rand(length(incl1),1);
rng(1);r4 = r+0.12 + (0.1)*rand(length(incl1),1);
rng(1);r5 = r+0.27 + (0.1)*rand(length(incl1),1);
sc1=scatter([r1 r2 r3 r4 r5],i_animal([2 3 4 5 6],incl1)',5,'k','filled');hold on;
for oo=1:5
sc1(oo).MarkerEdgeColor=[0.5 0.5 0.5];sc1(oo).MarkerFaceColor=[0.2 0.2 0.2];sc1(oo).MarkerFaceAlpha=[0.5];sc1(oo).MarkerEdgeAlpha=[0.5];
end
hold on;
r=2;
rng(1);r1 = r-0.35 + (0.1)*rand(length(incl1),1);
rng(1);r2 = r -0.2+ (0.1)*rand(length(incl1),1);
rng(1);r3 = r-0.05 + (0.1)* rand(length(incl1),1);
rng(1);r4 = r+0.12 + (0.1)*rand(length(incl1),1);
rng(1);r5 = r+0.27 + (0.1)*rand(length(incl1),1);
sc1=scatter([r1 r2 r3 r4 r5],c_animal([2 3 4 5 6],incl1)',5,'k','filled');hold on;
for oo=1:5
sc1(oo).MarkerEdgeColor=[0.5 0.5 0.5];sc1(oo).MarkerFaceColor=[0.2 0.2 0.2];sc1(oo).MarkerFaceAlpha=[0.5];sc1(oo).MarkerEdgeAlpha=[0.5];
end
  hold on;h=errorbar([1-0.3],[means_layers(1,1)],[err_layers(1,1)], 'LineStyle', 'none','Color', 'k', 'LineWidth', 2,'CapSize',0);
   hold on;h=errorbar([1-0.15],[means_layers(1,2)],[err_layers(1,2)], 'LineStyle', 'none','Color', 'k', 'LineWidth', 2,'CapSize',0);
 hold on;h=errorbar([1],[means_layers(1,3)],[err_layers(1,3)], 'LineStyle', 'none','Color', 'k', 'LineWidth', 2,'CapSize',0);
 hold on;h=errorbar([1+0.15],[means_layers(1,4)],[err_layers(1,4)], 'LineStyle', 'none','Color', 'k', 'LineWidth', 2,'CapSize',0);
 hold on;h=errorbar([1+0.3],[means_layers(1,5)],[err_layers(1,5)], 'LineStyle', 'none','Color', 'k', 'LineWidth', 2,'CapSize',0);
hold on;h=errorbar([2-0.3],[means_layers(2,1)],[err_layers(2,1)], 'LineStyle', 'none','Color', 'k', 'LineWidth', 2,'CapSize',0);
hold on;h=errorbar([2-0.15],[means_layers(2,2)],[err_layers(2,2)], 'LineStyle', 'none','Color', 'k', 'LineWidth', 2,'CapSize',0);
 hold on;h=errorbar([2],[means_layers(2,3)],[err_layers(2,3)], 'LineStyle', 'none','Color', 'k', 'LineWidth', 2,'CapSize',0);
 hold on;h=errorbar([2+0.15],[means_layers(2,4)],[err_layers(2,4)], 'LineStyle', 'none','Color', 'k', 'LineWidth', 2,'CapSize',0);
 hold on;h=errorbar([2+0.3],[means_layers(2,5)],[err_layers(2,5)], 'LineStyle', 'none','Color', 'k', 'LineWidth', 2,'CapSize',0);
 box off;xlim([0.5 2.5]);
ylabel('Fraction of neurons');set(gca,'FontSize',11);set(gca,'TickDir','out');box off;
%[p1]=signrank(i_animal(6,:),c_animal(6,:))
%[u p1]=ttest(i_animal(1,:),i_animal(5,:))

%text(0.75,0.6,['***'],'FontSize',18);hold on;text(3.75,0.6,['***'],'FontSize',18);
%legend({'ipsi','contra'},"Location","northeast");legend boxoff;
%offsetAxes;
h = gca;h.XAxis.Visible = 'off';  
text(0.68,0.62,'ipsi','FontSize',11,'Color',[0.7 0.7 0.7]);text(0.68,0.58,'contra','FontSize',11,'Color',[0.3 0.3 0.3]);
t1=text(0.68,-0.09,'L2/3','FontSize',11);set(t1,'Rotation',90);t1=text(0.84,-0.075,'L4','FontSize',11);set(t1,'Rotation',90);t1=text(1,-0.075,'L5','FontSize',11);set(t1,'Rotation',90);
t1=text(1.14,-0.09,'L6a','FontSize',11);set(t1,'Rotation',90);t1=text(1.30,-0.09,'L6b','FontSize',11);set(t1,'Rotation',90);
t1=text(1.68,-0.09,'L2/3','FontSize',11);set(t1,'Rotation',90);t1=text(1.84,-0.075,'L4','FontSize',11);set(t1,'Rotation',90);
t1=text(2,-0.075,'L5','FontSize',11);set(t1,'Rotation',90);t1=text(2.14,-0.09,'L6a','FontSize',11);set(t1,'Rotation',90);
t1=text(2.30,-0.09,'L6b','FontSize',11);set(t1,'Rotation',90);title(panel_tit{j},'FontWeight','normal','Color',temp_color(j,:));
end
cd(save_folder);saveas(gcf, 'layerFraction_injection_areas.pdf');

%% Panel C1 and C2 : Heatmaps, layer dominance L2/3, L5 and L6 sorted for ipsi and contra
name_sub={'VISpi_heatmaps.pdf', 'SSpi_heatmaps.pdf' ,'MOpi_heatmaps.pdf'};
%IPSI
%dat_all={};dat_all={c_v1aam c_s1aam c_m1aam};
%dat_all2={};dat_all2={v1a_rcam s1a_rcam m1a_rcam};
dat_all={};dat_all={i_v1aafm i_s1aafm i_m1aafm};
dat_all2={};dat_all2={v1a_riam s1a_riam m1a_riam};
temp_color=[v1_color ;s1_color; m1_color];
%panel_tit={'contra VISp','contra SSp-bfd','contra MOp'};
panel_tit={'VISp','SSp-bfd','MOp'};
cortex_abb=[];cortex_abb=cortex_names(:,2);

for j=1:3
p1=[];p1=nanmean(dat_all{j},3);
%p2=[];p2=nanmean(dat_all3{j},3)
idx_ke=[];idx_ke=find(sum(dat_all2{j}>0,2)>1);
% idx_ke2=find(sum(dat_all4{j}>0,2)>1);
p1e=[];p1e=p1(:,idx_ke);
% p2e=[];p2e=p2(:,idx_ke);
cortex_abb_e=[];
cortex_abb_e=cortex_abb(idx_ke);
ind=[];
[X, ind]=nanmax(p1e,[],1);
l6_an=[];l5_an=[];l23_an=[];
l6_an=find(ind==5);
l5_an=find(ind==4);
l23_an=find(ind==2);
kk1=[];kk2=[];kk3=[];
[l23_p1 kk1 ]=sort(p1e(2,l23_an));
[l5_p1 kk2]=sort(p1e(4,l5_an));
[l6_p1 kk3]=sort(p1e(5,l6_an));

tem_cm=[];tem_cm=[p1e(2:end,l23_an(kk1)) ones(5,2)*1 p1e(2:end,l5_an(kk2)) ones(5,2)*1 p1e(2:end,l6_an(kk3))];
%tem_cm2=[];tem_cm2=[p2e(2:end,l23_an(kk1)) zeros(5,2)*NaN p2e(2:end,l5_an(kk2)) zeros(5,2)*NaN p2e(2:end,l6_an(kk3))];
p3=[];p3=[cortex_abb_e(l23_an(kk1)) ;'.';'.';cortex_abb_e(l5_an(kk2));'.';'.' ;cortex_abb_e(l6_an(kk3))];

fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 1000-250*j, 750, 160]);set(gcf,'color','w');tiledlayout("horizontal")
nexttile
h=imagesc(tem_cm);
set(h, 'AlphaData', ~isnan(tem_cm))
cmap(temp_color(j,:),100,0.5,0.5);colorbar;
[hh oo]=find(isnan(tem_cm));
coord_b=[oo hh];
for u=1:length(coord_b)
    hold on
    text(coord_b(u,1)-0.3,coord_b(u,2)-0.15,'x','FontSize',14,'Color','k')
end
offsetAxes;

box off;xticks([1:length(tem_cm)]);
%ylim([1 5.5])
hold on;xticklabels(p3);yticklabels({'L2/3','L4','L5','L6a','L6b'})
title(panel_tit{j})
%text(-1,-0.3,'L2/3 dom');
text(4,-0.3,'L2/3 dom');

if j==1 
%text(10,-0.3,'L5 dom');
text(15,-0.3,'L5 dom');
%text(31,-0.3,'L6 dom'); 
text(35,-0.3,'L6 dom'); 
elseif j==2
  %text(14,-0.3,'L5 dom');  
  %text(36,-0.3,'L6 dom'); 
  text(17,-0.3,'L5 dom');
  text(40,-0.3,'L6 dom');
else
  %text(6,-0.3,'L5 dom');  
  text(12,-0.3,'L5 dom');  
  text(28,-0.3,'L6 dom');  
end
set(gca,'TickDir','out');box off;set(gca,'FontSize',10)
% fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 1400, 160]);set(gcf,'color','w');tiledlayout("horizontal")
% nexttile
% h=imagesc(tem_cm2);
% set(h, 'AlphaData', ~isnan(tem_cm2))
% cmap(temp_color(j,:),100,5,5);colorbar;box off;xticks([1:length(tem_cm)])
% hold on;xticklabels(p3);
% title(panel_tit{j})
cd(save_folder);saveas(gcf, name_sub{j});
end

%Contra
name_sub={'VISpc_heatmaps.pdf', 'SSpc_heatmaps.pdf' ,'MOpc_heatmaps.pdf'};
dat_all={};dat_all={c_v1aafm c_s1aafm c_m1aafm};
dat_all2={};dat_all2={v1a_rcam s1a_rcam m1a_rcam};
temp_color=[v1_color ;s1_color; m1_color];
panel_tit={'VISp','SSp-bfd','MOp'};
cortex_abb=[];cortex_abb=cortex_names(:,2);

for j=1:3
p1=[];p1=nanmean(dat_all{j},3);
%p2=[];p2=nanmean(dat_all3{j},3)
idx_ke=find(sum(dat_all2{j}>0,2)>1);
% idx_ke2=find(sum(dat_all4{j}>0,2)>1);
p1e=[];p1e=p1(:,idx_ke);
% p2e=[];p2e=p2(:,idx_ke);
cortex_abb_e=[];
cortex_abb_e=cortex_abb(idx_ke);
ind=[];
[X, ind]=nanmax(p1e,[],1);
l6_an=[];l5_an=[];l23_an=[];
l6_an=find(ind==5);
l5_an=find(ind==4);
l23_an=find(ind==2);
kk1=[];kk2=[];kk3=[];
[l23_p1 kk1 ]=sort(p1e(2,l23_an));
[l5_p1 kk2]=sort(p1e(4,l5_an));
[l6_p1 kk3]=sort(p1e(5,l6_an));

tem_cm=[];tem_cm=[p1e(2:end,l23_an(kk1)) ones(5,2)*1 p1e(2:end,l5_an(kk2)) ones(5,2)*1 p1e(2:end,l6_an(kk3))];
%tem_cm2=[];tem_cm2=[p2e(2:end,l23_an(kk1)) zeros(5,2)*NaN p2e(2:end,l5_an(kk2)) zeros(5,2)*NaN p2e(2:end,l6_an(kk3))];
p3=[];p3=[cortex_abb_e(l23_an(kk1)) ;'.';'.';cortex_abb_e(l5_an(kk2));'.';'.' ;cortex_abb_e(l6_an(kk3))];

fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [1300, 1000-250*j, 750, 160]);set(gcf,'color','w');tiledlayout("horizontal")
nexttile
h=imagesc(tem_cm);
set(h, 'AlphaData', ~isnan(tem_cm))
cmap(temp_color(j,:),100,0.5,0.5);colorbar;
[hh oo]=find(isnan(tem_cm));
coord_b=[oo hh];
for u=1:length(coord_b)
    hold on
    text(coord_b(u,1)-0.3,coord_b(u,2)-0.15,'x','FontSize',14,'Color','k')
end
offsetAxes;

box off;;xticks([1:length(tem_cm)]);

hold on;xticklabels(p3);xtickangle(90);yticklabels({'L2/3','L4','L5','L6a','L6b'})
title(panel_tit{j})
text(-1,-0.3,'L2/3 dom');

if j==1 
text(10,-0.3,'L5 dom');
text(31,-0.3,'L6 dom'); 

elseif j==2
  text(16,-0.3,'L5 dom');  
  text(36,-0.3,'L6 dom'); 
 
else
  text(5,-0.3,'L5 dom');  
  text(30,-0.3,'L6 dom');  
   
end
set(gca,'TickDir','out');box off;set(gca,'FontSize',10)
cd(save_folder);saveas(gcf, name_sub{j});
end
%% Panel D2 and D1 : Based on area VISp, SSpbf, MOp dominance of L2/3 L5 and L6a 
%CONTRA 
clc;
dat_all={};dat_all={c_v1aafm c_s1aafm c_m1aafm};
temp_color=[v1_color ;s1_color; m1_color];
panel_tit={'','contralateral',''};
all_ind={};
%%%%plot
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 350, 350, 250]);set(gcf,'color','w');t=tiledlayout("horizontal");
t.TileSpacing = 'compact';t.Padding = 'compact';
for j=1:3
    figure(fig7)
ind=[];temp_c=temp_color(j,:);temp_area=[];temp_area=dat_all{j};
    for i=1:size(temp_area,3)
    [X, ind(:,i)]=nanmax(temp_area(:,:,i),[],1);
    end
    all_ind{j}=ind;
   
    p1=[];p1=(sum(ind==2)./sum(ind>1))*100;
    p2=[];p2=(sum(ind==4)./sum(ind>1))*100;
    p3=[];p3=(sum(ind==5)./sum(ind>1))*100;
    temp_p= [p1' p2' p3'];
nexttile
title(panel_tit{j},'FontWeight','normal','Color',[0.3 0.3 0.3]);
    for i=1:3
    hold on;
    b1=bar(i,nanmean(temp_p(:,i)),0.7);hold on;b1.FaceColor=temp_c;b1.EdgeColor='none';
      r=i;rng(i);r1 = r-0.01 + (0.2)*rand(length(temp_p(:,i)),1);
    sc1=scatter(r1,temp_p(:,i),12,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
    hold on;errorbar([i],[nanmean(temp_p(:,i))],[nanstd(temp_p(:,i))/sqrt(length(temp_p(:,i)))]...
    , 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1,'CapSize',0);hold on;set(b1,'ShowBaseLine','off');
    end
    if j==1
xticks([1:3]);hold on;xticklabels({'L2/3','L5','L6a'});ylabel('Layer dominance per area (%)');
set(gca,'FontSize',11);set(gca,'TickDir','out');box off;
    else
        h = gca;h.YAxis.Visible = 'off';xticks([1:3]);hold on;xticklabels({'L2/3','L5','L6a'});
        set(gca,'FontSize',11);set(gca,'TickDir','out');box off;xtickangle(45);
    end

xlim([0.5 3.5]);
offsetAxes;h = gca;h.XAxis.Visible = 'off'
t1=text(0.9,-17,'L2/3','FontSize',11);set(t1,'Rotation',90);
t1=text(1.9,-17,'L5','FontSize',11);set(t1,'Rotation',90);
t1=text(2.9,-17,'L6a','FontSize',11);set(t1,'Rotation',90);
if j==2
line([2 3],[60 60],'Color','k');
text(2.05,65,'n.s.','Color','k','FontSize',11);
end
ylim([0 80]);
cd(save_folder);saveas(gcf, 'dominance_per_injection.pdf');
%stats
disp(['contra test all ' num2str(j)])
 [p,tbl,stats] = anova1([temp_p])
 presults=[];
  presults = multcompare(stats)
end

%% IPSI
clc;
dat_all={};dat_all={i_v1aafm i_s1aafm i_m1aafm};
temp_color=[v1_color ;s1_color; m1_color];
panel_tit={'','ipsilateral',''};
all_ind={};
%plot
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 350, 350, 250]);set(gcf,'color','w');t=tiledlayout("horizontal");
t.TileSpacing = 'compact';t.Padding = 'compact';
for j=1:3
    figure(fig7)
ind=[];temp_c=temp_color(j,:);temp_area=[];temp_area=dat_all{j};
    for i=1:size(temp_area,3)
    [X, ind(:,i)]=nanmax(temp_area(:,:,i),[],1);
    end
    all_ind{j}=ind;
   
    p1=[];p1=(sum(ind==2)./sum(ind>1))*100;
    p2=[];p2=(sum(ind==4)./sum(ind>1))*100;
    p3=[];p3=(sum(ind==5)./sum(ind>1))*100;
    temp_p= [p1' p2' p3'];
nexttile
title(panel_tit{j},'FontWeight','normal','Color',[0.7 0.7 0.7]);
    for i=1:3
    hold on;
    b1=bar(i,nanmean(temp_p(:,i)),0.7);hold on;b1.FaceColor=temp_c;b1.EdgeColor='none';
      r=i;rng(i);r1 = r-0.01 + (0.2)*rand(length(temp_p(:,i)),1);
    sc1=scatter(r1,temp_p(:,i),12,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
    hold on;errorbar([i],[nanmean(temp_p(:,i))],[nanstd(temp_p(:,i))/sqrt(length(temp_p(:,i)))]...
    , 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1,'CapSize',0);hold on;set(b1,'ShowBaseLine','off');
    end
    if j==1
xticks([1:3]);hold on;xticklabels({'L2/3','L5','L6a'});ylabel('Layer dominance per area (%)');
set(gca,'FontSize',11);set(gca,'TickDir','out');box off;
line([1 2],[40 40],'Color','k');
text(1.05,45,'n.s.','Color','k','FontSize',11);
    else
        h = gca;h.YAxis.Visible = 'off';xticks([1:3]);hold on;xticklabels({'L2/3','L5','L6a'});
        set(gca,'FontSize',11);set(gca,'TickDir','out');box off;xtickangle(45);
    end
if j==2
  line([1 3],[62 62],'Color','k');
text(1.6,67,'n.s.','Color','k','FontSize',11);  
end
if j==3
  line([1 2],[30 30],'Color','k');
text(1.05,35,'n.s.','Color','k','FontSize',11);  
end
xlim([0.5 3.5]);
offsetAxes;h = gca;h.XAxis.Visible = 'off'
t1=text(0.9,-17,'L2/3','FontSize',11);set(t1,'Rotation',90);
t1=text(1.9,-17,'L5','FontSize',11);set(t1,'Rotation',90);
t1=text(2.9,-17,'L6a','FontSize',11);set(t1,'Rotation',90);
ylim([0 80]);

disp(['ipsi test all ' num2str(j)])
 [p,tbl,stats] = anova1([temp_p])
 presults=[];
  presults = multcompare(stats)

end

%cd(save_folder);saveas(gcf, 'ipsi dominance_per_injection.pdf');
%% Panel E Plot the numbers as fraction per area where L2/3, L5 or L6a is dominant ACROSS ALL
clc
%all areas stacked
%IPSI
cat_temp=[];cat_temp=cat(3,i_v1aafm,i_s1aafm,i_m1aafm);
ind=[];
for i=1:size(cat_temp,3)
[X, ind(:,i)]=nanmax(cat_temp(:,:,i),[],1);
end
p1=[];p1=sum(ind==2)./sum(ind>1);
p2=[];p2=sum(ind==4)./sum(ind>1);
p3=[];p3=sum(ind==5)./sum(ind>1);
temp_p=[];temp_p= [p1' p2' p3'];


fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 200, 270]);set(gcf,'color','w');
b1=bar([1],[nanmean(temp_p(:,3)) nanmean(temp_p(:,2)) nanmean(temp_p(:,1))],'stacked');hold on;
b1(3).FaceColor=[1 1 1];b1(2).FaceColor=[0.7 0.7 0.7];b1(1).FaceColor=[0.55 0.55 0.55];
b1(3).LineWidth=1.2;b1(2).LineWidth=1.2;b1(1).LineWidth=1.2;

hold on;errorbar([1],[nanmean(temp_p(:,3)) nanmean(temp_p(:,2))+nanmean(temp_p(:,3)) nanmean(temp_p(:,1))+nanmean(temp_p(:,2))+nanmean(temp_p(:,3))],[nanstd(temp_p(:,3))/sqrt(length(temp_p(:,3))) nanstd(temp_p(:,2))/sqrt(length(temp_p(:,2))) nanstd(temp_p(:,2))/sqrt(length(temp_p(:,1)))]...
      , 'LineStyle', 'none', ... 
         'Color', 'k', 'LineWidth', 1.2,'CapSize',0);hold on;

text(0.665,0.85,num2str(round(nanmean(temp_p(:,1)),3)*100),'FontSize',8,'Color','k');
text(0.665,0.6,num2str(round(nanmean(temp_p(:,2)),3)*100),'FontSize',8,'Color','w');
text(0.665,0.25,num2str(round(nanmean(temp_p(:,3)),3)*100),'FontSize',8,'Color','w');

disp(['ipsi test all across '])
 [p,tbl,stats] = anova1([temp_p],[],'off')
 presults=[];
  presults = multcompare(stats,"Display","off")

%CONTRA
cat_temp=[];cat_temp=cat(3,c_v1aafm,c_s1aafm,c_m1aafm);
ind=[];
for i=1:size(cat_temp,3)
[X, ind(:,i)]=nanmax(cat_temp(:,:,i),[],1);
end
p1=[];p1=sum(ind==2)./sum(ind>1);
p2=[];p2=sum(ind==4)./sum(ind>1);
p3=[];p3=sum(ind==5)./sum(ind>1);
temp_p=[];temp_p= [p1' p2' p3'];

hold on;
b1=bar([3],[nanmean(temp_p(:,3)) nanmean(temp_p(:,2)) nanmean(temp_p(:,1))],'stacked');hold on;
b1(3).FaceColor=[1 1 1];b1(2).FaceColor=[0.7 0.7 0.7];b1(1).FaceColor=[0.55 0.55 0.55];
b1(3).LineWidth=1.2;b1(2).LineWidth=1.2;b1(1).LineWidth=1.2;

hold on;errorbar([3],[nanmean(temp_p(:,3)) nanmean(temp_p(:,2))+nanmean(temp_p(:,3)) nanmean(temp_p(:,1))+nanmean(temp_p(:,2))+nanmean(temp_p(:,3))],[nanstd(temp_p(:,3))/sqrt(length(temp_p(:,3))) nanstd(temp_p(:,2))/sqrt(length(temp_p(:,2))) nanstd(temp_p(:,1))/sqrt(length(temp_p(:,1)))]...
      , 'LineStyle', 'none', ... 
         'Color', 'k', 'LineWidth', 1.2,'CapSize',0);hold on;

hold on;box off;xticks([1:2:3]);xticklabels({'ipsi','contra'});
t1=text(0,0.2,{'Layer dominance ','per cortical area (%)'},'FontSize',11);set(t1,'Rotation',90);
set(gca,'FontSize',11);set(gca,'TickDir','out');
h = gca;h.YAxis.Visible = 'off';  

text(2.665,0.95,num2str(round(nanmean(temp_p(:,1)),3)*100),'FontSize',8,'Color','k');
text(2.665,0.71,num2str(round(nanmean(temp_p(:,2)),3)*100),'FontSize',8,'Color','w');
text(2.665,0.25,num2str(round(nanmean(temp_p(:,3)),3)*100),'FontSize',8,'Color','w');


t1=text(3.7,0.98,'L2/3','FontSize',11);set(t1,'Rotation',270);
t1=text(3.7,0.75,'L5','FontSize',11);set(t1,'Rotation',270);
t1=text(3.7,0.35,'L6a','FontSize',11);set(t1,'Rotation',270);

%Overall
% hold on;text(1.8,1.15,'All','FontSize',12);
% ax = gca; ax.XColor = 'w'; % Red
% hold on;line([0.55 1.45],[0 0],'Color','k','LineWidth', 1);
% hold on;line([2.55 3.45],[0 0],'Color','k','LineWidth', 1);
h = gca;h.XAxis.Visible = 'off';set(b1,'ShowBaseLine','off');
hold on;text(0.65,-0.065,'ipsi','FontSize',11);
hold on;text(2.4,-0.065,'contra','FontSize',11);
hold on;text(1.88,0.68,'*','FontSize',18);
line([2.3 2.6],[0.95 0.95],'Color','k');
line([2.3 2.6],[0.7 0.7],'Color','k');
line([2.3 2.6],[0.25 0.25],'Color','k');
line([2.3 2.3],[0.25 0.95],'Color','k');

line([1.7 1.4],[0.85 0.85],'Color','k');
line([1.7 1.4],[0.6 0.6],'Color','k');
line([1.7 1.4],[0.2 0.2],'Color','k');
line([1.7 1.7],[0.2 0.85],'Color','k');

cd(save_folder);saveas(gcf,'overall_dominance_2.pdf');

disp(['contra test all across '])
 [p,tbl,stats] = anova1([temp_p],[],'off')
 presults=[];
  presults = multcompare(stats,"Display","off")





%% Figure 4: Contralateral input onto VISp, SSp-bfd and MOp displays high cortical hierarchy
%Figure 4A: illustrator 
%Figure 4B: illustrator 
%% Panel B Next to each other hemispheric fILN ipsi vs contra 
clc;
idx=1;%fILN
temp_color=[v1_color ;s1_color; m1_color];
module_names={'VISp','SSp-bfd','MOp'};
all_dati={};all_dati={iv1_index(:,idx) is1_index(:,idx) im1_index(:,idx)};
all_datc={};all_datc={cv1_index(:,idx) cs1_index(:,idx) cm1_index(:,idx)};
%plot
fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 400, 250]);set(gcf,'color','w');t=tiledlayout("horizontal");
t.TileSpacing = 'compact';t.Padding = 'compact';
for j=1:3
    nexttile
 
dat=[];dat=[all_dati{:,j} all_datc{:,j}];
for i=1:length(data)
     pl=plot([1,2],[dat(:,1),dat(:,2)],'color',[0.5 0.5 0.5]);    
end

hold on;pS=plotSpread([dat(:,1),dat(:,2)],'categoryIdx',[ones(1,length(dat(:,1)))' ones(1,length(dat(:,2)))'*2],...
    'categoryMarkers',{'o','o'},'categoryColors',{temp_color(j,:), temp_color(j,:)});hold on;
hold on;er1=errorbar([0.75],nanmean(dat(:,1)),nanstd(dat(:,1),[],1)/sqrt(length(dat(:,1))),'ok','MarkerFaceColor',temp_color(j,:),'Markersize',8);
hold on;er2=errorbar([2.25],nanmean(dat(:,2)),nanstd(dat(:,2),[],1)/sqrt(length(dat(:,1))),'ok','MarkerFaceColor',temp_color(j,:),'Markersize',8);
%stats
disp(['test ipsi vs contra ' module_names{j}])
[u p1]=ttest(dat(:,1),dat(:,2),'Tail','left')
nanmean(dat(:,1))
nanstd(dat(:,1))/(sqrt(length(dat(:,1))))
nanmean(dat(:,2))
nanstd(dat(:,2))/(sqrt(length(dat(:,2))))

offsetAxes
xticklabels({'ipsi','contra'});ylabel('fILN');hold on;title([]);xtickangle(45);set(gca,'FontSize',11);set(gca,'TickDir','out');box off;
hold on;text(1.5,0.9,['*'],'FontSize',18);ylim([0.4 0.9])
 h = gca;h.XAxis.Visible = 'off'
 t1=text(0.75,0.4,'ipsi','FontSize',11);set(t1,'Rotation',45);
  t1=text(1.75,0.4,'contra','FontSize',11);set(t1,'Rotation',45);
  line([-0.1 3],[0.5 0.5],'Color','k','LineStyle','--')
  if j==2 | j==3
     h = gca;h.YAxis.Visible = 'off'  
  end
end
cd(save_folder);saveas(gcf, 'ILN_global.pdf');



%% supplement 2 to Fig4: Sorted areas based on index for supplement 
temp_color2=linspecer(6);
dat_all1={};dat_all1={i_v1aafm i_s1aafm i_m1aafm};
dat_all2={};dat_all2={c_v1aafm c_s1aafm c_m1aafm};
title_pan={'VISp','SSp-bfd','MOp'};
cl_idx=[visp_idx ssp_idx mop_idx];
frontal_idx=[1:8];lateral_idx=[9:16];somamo_idx=[17:26];visual_idx=[27:36];medial_idx=[37:39];aud_idx=[40:43];lateral_idx2=[44:45];
areas_idx=ones(45,1);
areas_modul=[];
areas_modul=[areas_idx(frontal_idx)*1; areas_idx(lateral_idx)*2 ; areas_idx(somamo_idx)*3; ...
    areas_idx(visual_idx)*4; areas_idx(medial_idx)*5; areas_idx(aud_idx)*6 ; areas_idx(lateral_idx2)*2];
temp_c1=[0.8 0.8 0.8];
temp_c2=[0.3 0.3 0.3];
temp_color=[v1_color ;s1_color; m1_color];
idx=[];idx=1;
aILNi=[];aILNc=[];

for k=1:3
    rm1=[];rm2=[];rm1=dat_all1{k};rm2=dat_all2{k};
%calculate and plot 
i_metric=[];c_metric=[];
%calculate indexes across areas
for i=1:length(rm1)
temp_l=[];temp_l2=[];
temp_l=squeeze(rm1(:,i,:));
temp_l2=squeeze(rm2(:,i,:));
[i_metric(:,i,:)] = anatomy_indexcalc(temp_l);
[c_metric(:,i,:)] = anatomy_indexcalc(temp_l2);
end 
temp_metric1=i_metric(:,:,idx);temp_metric2=c_metric(:,:,idx);
rm1=[];rm2=[];
rm1=temp_metric1;rm2=temp_metric2;
aILNi=[aILNi; rm1];
aILNc=[aILNc; rm2];

temp_rm1=[];temp_rm1=nanmean(rm1);
temp_rm2=[];temp_rm2=nanmean(rm2);
temp_rm1(cl_idx(k))=[];
temp_rm2(cl_idx(k))=[];
cortex_abb=[];
cortex_abb=cortex_names(:,2);
cortex_abb(cl_idx(k))=[];
yerr=[];yerr=nanstd(rm1)./(sqrt(sum(~isnan(rm1))));
yerr2=[];yerr2=nanstd(rm2)./(sqrt(sum(~isnan(rm2))));

yerr(cl_idx(k))=[];
yerr2(cl_idx(k))=[];
label_all=[];
label_all=areas_modul;
label_all(cl_idx(k))=[];
sort_d=[];kk=[];
[sort_d kk]=sort(temp_rm1,'ascend');

data_ipsi=[];data_contra=[];data_contra2=[];yerr_ipsi=[];yerr_contra=[];yerr_ipsi2=[];yerr_contra2=[];
cortex_abb2=[];
cortex_abb2=cortex_abb(kk);

data_contra=temp_rm2(kk);
both_pre=[];
both_pre=find(~isnan(data_contra));
data_ipsi=sort_d(both_pre);
data_contra2=data_contra(both_pre);
label_a=[];label_a=label_all(kk);
label_module=[];label_module=label_a(both_pre);
yerr_ipsi=yerr(kk);
yerr_ipsi2=yerr_ipsi(both_pre);
yerr_contra=yerr2(kk);
yerr_contra2=yerr_contra(both_pre);
cortex_fnames=[];
cortex_fnames=cortex_abb2(both_pre);


alpha1=1;
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, -200+k*300, 850, 350]);set(gcf,'color','w');

title(title_pan{k},'FontWeight','normal','Color',temp_color(k,:));ylabel('fILN');
hold on;h=errorbar([1:length(both_pre)],[data_ipsi],yerr_ipsi2, 'LineStyle', 'none', ... 
        'Color', 'k' ,'LineWidth', 0.5,'CapSize',0);set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha1]);
hold on;h=errorbar([1:length(both_pre)],[data_contra2],yerr_contra2...
   , 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 0.5,'CapSize',0);
%hold on;h=errorbar([45],[nanmean(temp_metric2(:,cl_idx(k)))],[yerr3]...
  %  , 'LineStyle', 'none', ... 
 %       'Color', temp_color(k,:), 'LineWidth', 0.5,'CapSize',2);
li=line([0 length(both_pre)],[0.5 0.5],'Color',[0.5 0.5 0.5],'LineStyle',':');
 xticklabels([cortex_fnames]);
 xtickangle(90)
for m=1:length(both_pre)
    hold on;pp1=scatter([m],[data_ipsi(m)],50,'o','LineWidth',2);
    hold on;pp2=scatter([[m]],[data_contra2(m)],50,'o','LineWidth',2);%pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor=temp_c2;box off;
    if label_module(m)==1
    pp1.MarkerEdgeColor=temp_color2(1,:)
    pp1.MarkerFaceColor=temp_color2(1,:);
    pp2.MarkerEdgeColor=temp_color2(1,:);
    pp2.MarkerFaceColor='w';
    elseif label_module(m)==2
    pp1.MarkerEdgeColor=temp_color2(2,:)
    pp1.MarkerFaceColor=temp_color2(2,:);
    pp2.MarkerEdgeColor=temp_color2(2,:);
    pp2.MarkerFaceColor='w';
    elseif label_module(m)==3
     pp1.MarkerEdgeColor=temp_color2(3,:)
    pp1.MarkerFaceColor=temp_color2(3,:);
    pp2.MarkerEdgeColor=temp_color2(3,:);
    pp2.MarkerFaceColor='w';
    elseif label_module(m)==4
    pp1.MarkerEdgeColor=temp_color2(4,:)
    pp1.MarkerFaceColor=temp_color2(4,:);
    pp2.MarkerEdgeColor=temp_color2(4,:);
    pp2.MarkerFaceColor='w';
    elseif label_module(m)==5
    pp1.MarkerEdgeColor=temp_color2(5,:)
    pp1.MarkerFaceColor=temp_color2(5,:);
    pp2.MarkerEdgeColor=temp_color2(5,:);
    pp2.MarkerFaceColor='w';
    else label_module(m)==6
    pp1.MarkerEdgeColor=temp_color2(6,:)
    pp1.MarkerFaceColor=temp_color2(6,:);
    pp2.MarkerEdgeColor=temp_color2(6,:);
    pp2.MarkerFaceColor='w';
    end
end
xticks([1:length(both_pre)]);
ylim([0 1]);%xticks([0:0.2:1.2])
%xticklabels([cortex_abb(kk) ;{title_pan{k}}])
offsetAxes;

     %h = gca;h.XAxis.Visible = 'off'  
   
     text(42,0.4,'FF','FontSize',11);
     text(42,0.6,'FB','FontSize',11);
    pp3=scatter(5,0.1,50); pp3.MarkerEdgeColor='k';pp3.MarkerFaceColor='k';
    pp3=scatter(5,0.05,50); pp3.MarkerEdgeColor='k';pp3.MarkerFaceColor='w';
     text(6,0.1,'ipsi','Color',[0.7 0.7 0.7],'FontSize',11);
text(6,0.05,'contra','Color',[0.3 0.3 0.3],'FontSize',11);

%for m=1:length(cortex_abb(kk))
 %yticklabels([cortex_abb(kk)]);
%ax = gca;

% Simply color an XTickLabel



%end
set(gca,'FontSize',11);set(gca,'TickDir','out');
offsetAxes;
%view([90 -270]) %// instead of normal view, which is view([0 90])
if k==1
cd(save_folder);saveas(gcf, 'ILNipsi_sorted_areasVISP.pdf');
end
if k==2
    cd(save_folder);saveas(gcf, 'ILNipsi_sorted_areasSSP.pdf');
end
if k==3

      %text(10,-0.1,'Cortical Areas (sorted by ipsi fILN)','FontSize',11);
      hold on;
          %text(33,0.45,'Modules','FontSize',11,'Color','k');
       text(33,0.55,'Prefrontal','FontSize',11,'Color',temp_color2(1,:));
       text(33,0.45,'Lateral','FontSize',11,'Color',temp_color2(2,:));
       text(33,0.35,'Somatomotor','FontSize',11,'Color',temp_color2(3,:));
       text(33,0.25,'Visual','FontSize',11,'Color',temp_color2(4,:));
       text(33,0.15,'Medial','FontSize',11,'Color',temp_color2(5,:));
       text(33,0.05,'Auditory','FontSize',11,'Color',temp_color2(6,:));  
    cd(save_folder);saveas(gcf, 'ILNipsi_sorted_areasMOP.pdf');
end
end

%%  Panel D: ipsi sorted, 18 animals average 
clc;
temp_rm1=[];temp_rm1=nanmean(aILNi);
temp_rm2=[];temp_rm2=nanmean(aILNc);
temp_c1=[0.8 0.8 0.8];
temp_c2=[0.3 0.3 0.3];
cortex_abb=[];cortex_abb=cortex_names(:,2);

yerr=[];yerr=nanstd(aILNi)./(sqrt(sum(~isnan(aILNi))));
yerr2=[];yerr2=nanstd(aILNc)./(sqrt(sum(~isnan(aILNc))));

%yerr3=[];yerr3=nanstd(rm2(:,cl_idx(k)))/sqrt(length(rm2(:,cl_idx(k))));
sort_d=[];kk=[];
[sort_d kk]=sort(temp_rm1,'ascend');
cortex_lab=[];
for oo=1:45
cortex_lab{oo}=[cortex_abb{kk(oo)} ' ' num2str(kk(oo))];
end
alpha1=1;
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 300, 900, 300]);set(gcf,'color','w');

title('','FontWeight','normal','Color','k');
ylabel('fILN')
li=line([1 45],[0.5 0.5],'Color',[0.5 0.5 0.5],'LineStyle',':');
hold on;plot(sort_d,'Color',temp_c1);
hold on;plot(temp_rm2(kk),'Color',temp_c2);
hold on;h=errorbar([1:45],[sort_d],[yerr(kk)]...
    , 'LineStyle', 'none', ... 
        'Color', temp_c1, 'LineWidth', 0.5,'CapSize',0);set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha1]);
hold on;h=errorbar([1:45],[temp_rm2(kk)],[yerr2(kk)]...
    , 'LineStyle', 'none', ... 
        'Color', temp_c2, 'LineWidth', 0.5,'CapSize',0);

hold on;pp1=scatter(1:45,sort_d,50,'o');pp1.MarkerEdgeColor=[1 1 1];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor=temp_c1;
hold on;pp1=scatter(1:45,temp_rm2(kk),50,'o');pp1.MarkerEdgeColor=[1 1 1];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor=temp_c2;box off;

xticks([1:1:45]);ylim([0.3 1])
xticklabels(cortex_lab)
ylabel('fILN');xlabel('Cortical Areas (sorted by ipsi fILN)')
set(gca,'FontSize',11);set(gca,'TickDir','out');

offsetAxes;

     %h = gca;h.XAxis.Visible = 'off'  
     %text(20,-0.1,'Cortical Areas (sorted by ipsi fILN)','FontSize',11);
     text(42,0.4,'FF','FontSize',11);
     text(42,0.6,'FB','FontSize',11);
     text(3,0.95,'ipsi','Color',[0.7 0.7 0.7],'FontSize',11);
text(3,1,'contra','Color',[0.3 0.3 0.3],'FontSize',11);
     cd(save_folder);saveas(gcf, 'ILNipsi_sorted_areasALL.pdf');
  disp('ipsi ILN range')
  sort_d(1)
  tg=yerr(kk);tg(1)
  sort_d(end)
  tg(end)
    disp('contra ILN range')
nanmin(temp_rm2)
tg2=yerr2;tg2(find(nanmin(temp_rm2)==temp_rm2))
nanmax(temp_rm2)
tg2(find(nanmax(temp_rm2)==temp_rm2))


cd(save_folder);saveas(gcf, 'fiLNipsi_contra_sorted.pdf');
%% Panel E: correlation between indexes per area  ILN
  %all correlations superimposed on each other 
  clc;
idx=[];idx=1;
temp_all={i_v1aafm c_v1aafm; i_s1aafm c_s1aafm; i_m1aafm c_m1aafm};
temp_color=[v1_color ;s1_color; m1_color];
%plot
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 350, 350]);set(gcf,'color','w');
for i=1:3
  
    rm1=[];rm1=temp_all{i,1};rm2=[];rm2=temp_all{i,2};
[vali_all valdci_all val_i val_dci xerr yerr r_a p_a r_avg p_avg] = anatomy_correlation(rm1, rm2, idx);
 temp_c=temp_color(i,:);
 col_ra2(:,i)=r_a;
hold on;errorbar(val_i,val_dci,yerr,yerr,xerr,xerr, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 0.5,'CapSize',0);hold on;
sc1=scatter(val_i,val_dci,50,'filled');sc1.MarkerFaceColor=temp_c;sc1.MarkerEdgeColor='w';
hold on;
hold on;text(0.65,0.85-(0.08*i),...
    ['r= ' num2str(round(nanmean(r_a),2)) ' +- ' num2str(round(nanstd(r_a)/sqrt(length(r_a)),2))],'FontSize',11,'Color',temp_c);
%hold on;text(0.5,0.5-(0.05*i),['p= ' num2str(round(nanmean(p_a),2))],'FontSize',11,'Color',temp_c);
disp(['pvalues' num2str(i)])
p_a
end
%ylim([-0.3 0.8])
hold on;line([0.5 0.5],[-0.4 0.8],'Color','k','LineStyle','--')
xlim([0.15 1])
%ylabel('L6d_c - L6d_i');xlabel('L6d_i');
ylabel('fILN contra - fILN ipsi');xlabel('fILN ipsi');
axis square;set(gca,'FontSize',11);set(gca,'TickDir','out');box off;
offsetAxes
cd(save_folder);saveas(gcf, 'correlation ILNdelta.pdf');
%% Panel F: fILn contra vs IPSI, added now to Figure 4 
temp_color2=linspecer(6);
dat_all1={};dat_all1={i_v1aafm i_s1aafm i_m1aafm};
dat_all2={};dat_all2={c_v1aafm c_s1aafm c_m1aafm};
title_pan={'VISp','SSp-bfd','MOp'};

temp_c1=[0.8 0.8 0.8];
temp_c2=[0.3 0.3 0.3];
temp_color=[v1_color ;s1_color; m1_color];
idx=[];idx=1;
aILNi=[];
aILNc=[];


%Plot total number of cells vs L6 fraction per area 
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 800, 400]);set(gcf,'color','w');t=tiledlayout('horizontal');
t.TileSpacing = 'compact';t.Padding = 'compact';
for k=1:3
    nexttile
    rm1=[];rm2=[];rm1=dat_all1{k};rm2=dat_all2{k};
    temp_c=temp_color(k,:);
%calculate and plot 
i_metric=[];c_metric=[];
%calculate indexes across areas
for i=1:length(rm1)
temp_l=[];temp_l2=[];
temp_l=squeeze(rm1(:,i,:));
temp_l2=squeeze(rm2(:,i,:));
[i_metric(:,i,:)] = anatomy_indexcalc(temp_l);
[c_metric(:,i,:)] = anatomy_indexcalc(temp_l2);
end 
temp_metric1=i_metric(:,:,idx);temp_metric2=c_metric(:,:,idx);
rm1=[];rm2=[];
rm1=temp_metric1;rm2=temp_metric2;
aILNi=[aILNi; rm1];
aILNc=[aILNc; rm2];

temp_rm1=[];temp_rm1=nanmean(rm1);
temp_rm2=[];temp_rm2=nanmean(rm2);

yerr=[];yerr=nanstd(rm1)./(sqrt(sum(~isnan(rm1))));
yerr2=[];yerr2=nanstd(rm2)./(sqrt(sum(~isnan(rm2))));

hold on;errorbar(temp_rm1,temp_rm2,yerr2,yerr2,yerr,yerr, 'LineStyle', 'none', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5,'CapSize',0);hold on;
r_a=[];p_a=[];
for i=1:6
r=[];p=[];
[r p] = corr(rm1(i,:)',rm2(i,:)','Type','Spearman','Rows','complete');
r_a(i)=r(1);
p_a(i)=p(1);
end
p_a
combined_p = harmonic_mean_pvalue(p_a)
hold on;text(0.2,0.4,...
    ['r= ' num2str(round(nanmean(r_a),2)) ' +- ' num2str(round(nanstd(r_a)/sqrt(length(r_a)),2))],'FontSize',11,'Color',temp_c);
hold on;
sc1=scatter(temp_rm1,temp_rm2,50,'filled');sc1.MarkerFaceColor=temp_c;sc1.MarkerEdgeColor='w';
hold on;
axis square;
xlim([0.1 1]);ylim([0.1 1])
ylabel('fILN contra');xlabel('fILN ipsi')
 box off;set(gca,'FontSize',11);set(gca,'TickDir','out');
 offsetAxes;
 %set(gca, 'XScale', 'log');set(gca, 'YScale', 'log')
end
cd(save_folder);saveas(gcf, 'ipsiILN_vs_contraILN.pdf');
%% Panel G: Plot circles  
%temp_color2=linspecer(6);
dat_all1={};dat_all1={i_v1aafm i_s1aafm i_m1aafm};
dat_all2={};dat_all2={c_v1aafm c_s1aafm c_m1aafm};
title_pan={'VISp','SSp-bfd','MOp'};
cl_idx=[visp_idx ssp_idx mop_idx];
frontal_idx=[1:8];lateral_idx=[9:16];somamo_idx=[17:26];visual_idx=[27:36];medial_idx=[37:39];aud_idx=[40:43];lateral_idx2=[44:45];
areas_idx=ones(45,1);
areas_modul=[];
areas_modul=[areas_idx(frontal_idx)*1; areas_idx(lateral_idx)*2 ; areas_idx(somamo_idx)*3; ...
    areas_idx(visual_idx)*4; areas_idx(medial_idx)*5; areas_idx(aud_idx)*6 ; areas_idx(lateral_idx2)*2];
temp_c1=[0.8 0.8 0.8];
temp_c2=[0.3 0.3 0.3];
temp_color=[v1_color ;s1_color; m1_color];
idx=[];idx=1;
aILNi=[];
aILNc=[];

for k=1:3
    rm1=[];rm2=[];rm1=dat_all1{k};rm2=dat_all2{k};
%calculate and plot 
i_metric=[];c_metric=[];
%calculate indexes across areas
for i=1:length(rm1)
temp_l=[];temp_l2=[];
temp_l=squeeze(rm1(:,i,:));
temp_l2=squeeze(rm2(:,i,:));
[i_metric(:,i,:)] = anatomy_indexcalc(temp_l);
[c_metric(:,i,:)] = anatomy_indexcalc(temp_l2);
end 
temp_metric1=i_metric(:,:,idx);temp_metric2=c_metric(:,:,idx);
rm1=[];rm2=[];
rm1=temp_metric1;rm2=temp_metric2;
aILNi=[aILNi; rm1];
aILNc=[aILNc; rm2];

temp_rm1=[];temp_rm1=nanmean(rm1);
temp_rm2=[];temp_rm2=nanmean(rm2);
% temp_rm1(cl_idx(k))=[];
% temp_rm2(cl_idx(k))=[];
cortex_abb=[];
cortex_abb=cortex_names(:,2);
fILN_visualize(cell2mat(cortex_names(:,4)),cell2mat(cortex_names(:,5)),temp_rm1,temp_rm2,temp_color(k,:),title_pan{k})
if k==1
cd(save_folder);saveas(gcf, 'cortex_iln_VISp.pdf');
elseif k==2
    cd(save_folder);saveas(gcf, 'cortex_iln_SSp.pdf');
else k==3
    cd(save_folder);saveas(gcf, 'cortex_iln_MOp.pdf');
end
% cortex_abb(cl_idx(k))=[];
end

%% Figure 5:  Sensory and motor but not lateral, prefrontal and medial projections explain asymmetric cortical hierarchy
%% Panel A: ILN ipsi and ILN contra across all animals across 6 modules  
clc;
frontal_idx=[1:8];lateral_idx=[9:16 44:45];somamo_idx=[17:26];visual_idx=[27:36];medial_idx=[37:39];aud_idx=[40:43];
idx_modules={frontal_idx lateral_idx somamo_idx visual_idx medial_idx aud_idx};
dat_all1={};dat_all1={i_v1aafm i_s1aafm i_m1aafm};
dat_all2={};dat_all2={c_v1aafm c_s1aafm c_m1aafm};
senso_mo=[idx_modules{3} idx_modules{4} idx_modules{6}];
higher_a=[idx_modules{1} idx_modules{2} idx_modules{5}];

module_names={'Pfron','Lat','SoMo','Vis','Med','Aud'};
pan_title={'VISp','SSp-bfd','MOp'};
temp_color=[v1_color ;s1_color; m1_color];
temp_c1='r';
temp_c2=[0.3 0.3 0.3];
idx=[];idx=1;

all_metric1=[];all_metric2=[];
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 250, 450]);set(gcf,'color','w');
%plot
for j=1:3
rm1=[];rm2=[];rm1=dat_all1{j};rm2=dat_all2{j};
%calculate and plot 
i_metric=[];c_metric=[];
%calculate indexes across areas
for i=1:length(rm1)
temp_l=[];temp_l2=[];
temp_l=squeeze(rm1(:,i,:));
temp_l2=squeeze(rm2(:,i,:));
[i_metric(:,i,:)] = anatomy_indexcalc(temp_l);
[c_metric(:,i,:)] = anatomy_indexcalc(temp_l2);
end 
temp_metric1=i_metric(:,:,idx);temp_metric2=c_metric(:,:,idx);
all_metric1=[all_metric1; temp_metric1];
all_metric2=[all_metric2; temp_metric2];
end

    hold on;
for i=1:6
    ymean1=[];ymean1=nanmean(all_metric1(:,idx_modules{i}));
    ymean2=[];ymean2=nanmean(all_metric2(:,idx_modules{i}));
     pp1=scatter(1,nanmean(ymean1),80,'filled');     pp1.MarkerEdgeColor=[1 1 1]
     if i==1 | i==2 | i==5
pp1.MarkerFaceColor=[0.7 0.7 0.7];

     else
         pp1.MarkerFaceColor='r';
     end
     if i==2
         hold on;text(1.25,nanmean(ymean1)+0.015,module_names{i});
     elseif i==1
          hold on;text(1.25,nanmean(ymean1)+0.006,module_names{i});
     else
hold on;text(1.25,nanmean(ymean1),module_names{i});
     end
     
        pp1=scatter(3,nanmean(ymean2),80,'filled');pp1.MarkerEdgeColor=[1 1 1];
         if i==1 | i==2 | i==5 
        pp1.MarkerFaceColor=temp_c2;
        pp1.MarkerFaceColor=[0.7 0.7 0.7];
     else
         pp1.MarkerFaceColor='r';
         end
         if i==2
    hold on;text(3.25,nanmean(ymean2)+0.015,module_names{i});
      elseif i==1
          hold on;text(3.25,nanmean(ymean2)-0.006,module_names{i});
         else
hold on;text(3.25,nanmean(ymean2),module_names{i});
         end
    xlim([0 6]);
    ylim([0.55 0.95]);
end
    box off;
    ylabel('fILN');set(gca,'FontSize',11);set(gca,'TickDir','out');
    yticks([0.55:0.1:0.95]);xticks([1:2:3]);xticklabels({'ipsi','contra'})
    offsetAxes
h = gca;h.XAxis.Visible = 'off';
t1=text(0.7,0.5,'ipsi','FontSize',11,'Color',[0.7 0.7 0.7]);
 set(t1,'Rotation',45);
 t1=text(2.8,0.5,'contra','FontSize',11,'Color',[0.3 0.3 0.3]);
 set(t1,'Rotation',45); 
line([4.9 4.9],[0.89 0.95],'Color',[0.5 0.5 0.5]);
line([4.9 4.9],[0.75 0.82],'Color','r');
line([4.9 5.3],[0.7850 0.7850],'Color','r');
line([4.9 5.3],[0.92 0.92],'Color',[0.5 0.5 0.5]);
text(5.35,0.7850,'v-a-sm','Color','r');
text(5.35,0.92,'pf-m-l','Color',[0.5 0.5 0.5]);
cd(save_folder);saveas(gcf, 'ILN_modules_hierarchy.pdf');

%% Panel B: two groups ipsi vs contra ILN absolute
clc;
frontal_idx=[1:8];lateral_idx=[9:16 44:45];somamo_idx=[17:26];visual_idx=[27:36];medial_idx=[37:39];aud_idx=[40:43];
idx=1;
rm1=[];rm1=cat(3,i_v1aafm,i_s1aafm,i_m1aafm);
rm2=[];rm2=cat(3,c_v1aafm,c_s1aafm,c_m1aafm);
temp_color=[v1_color ;s1_color; m1_color];
senso_mo=[idx_modules{3} idx_modules{4} idx_modules{6}];
higher_a=[idx_modules{1} idx_modules{2} idx_modules{5}];
[vali_all valdci_all val_i val_dci xerr yerr r_a p_a r_avg p_avg valc_all xerr_c] = anatomy_correlation(rm1, rm2, idx);


% 
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 350, 220]);set(gcf,'color','w');
b2=bar(2,[nanmean(nanmean(vali_all(:,senso_mo),2)) ;nanmean(nanmean(valc_all(:,senso_mo),2))],0.7);
b2(1).FaceColor=[0.8 0.8 0.8];hold on;b2(1).EdgeColor='r';b2(2).FaceColor=[0.55 0.55 0.55];hold on;b2(2).EdgeColor='r';set(b2,'ShowBaseLine','off');
b1=bar(1,[nanmean(nanmean(vali_all(:,higher_a),2)) ;nanmean(nanmean(valc_all(:,higher_a),2))],0.7);set(b1,'ShowBaseLine','off');set(b1,'ShowBaseLine','off');
b1(1).FaceColor=[0.8 0.8 0.8];hold on;b1(1).EdgeColor=[0.3 0.3 0.3];b1(2).FaceColor=[0.55 0.55 0.55];hold on;b1(2).EdgeColor=[0.3 0.3 0.3];

%indivdual animals

r=1;
rng(1);r1 = r-0.2 + (0.1)*rand(length(nanmean(valc_all(:,senso_mo),2)),1);
rng(1);r2 = r + 0.05 + (0.25)*rand(length(nanmean(valc_all(:,senso_mo),2)),1);
sc1=scatter([r1],[nanmean(vali_all(:,higher_a),2)],5,'k','filled');hold on;
sc1(1).MarkerEdgeColor=[0.5 0.5 0.5];sc1(1).MarkerFaceColor=[0.2 0.2 0.2];sc1(1).MarkerFaceAlpha=[0.5];sc1(1).MarkerEdgeAlpha=[0.5];
sc1=scatter([r2],[nanmean(valc_all(:,higher_a),2)],5,'k','filled');hold on;
sc1(1).MarkerEdgeColor=[0.5 0.5 0.5];sc1(1).MarkerFaceColor=[0.2 0.2 0.2];sc1(1).MarkerFaceAlpha=[0.5];sc1(1).MarkerEdgeAlpha=[0.5];
r=2;
rng(1);r1 = r-0.2 + (0.1)*rand(length(nanmean(valc_all(:,senso_mo),2)),1);
rng(1);r2 = r+0.05 +  (0.25)*rand(length(nanmean(valc_all(:,senso_mo),2)),1);
sc1=scatter([r1],[nanmean(vali_all(:,senso_mo),2)],5,'k','filled');hold on;
sc1(1).MarkerEdgeColor=[0.5 0.5 0.5];sc1(1).MarkerFaceColor=[0.2 0.2 0.2];sc1(1).MarkerFaceAlpha=[0.5];sc1(1).MarkerEdgeAlpha=[0.5];
sc1=scatter([r2],[nanmean(valc_all(:,senso_mo),2)],5,'k','filled');hold on;
sc1(1).MarkerEdgeColor=[0.5 0.5 0.5];sc1(1).MarkerFaceColor=[0.2 0.2 0.2];sc1(1).MarkerFaceAlpha=[0.5];sc1(1).MarkerEdgeAlpha=[0.5];
hold on;   
xpos = b2(1).XData + b2(1).XOffset;
errorbar(1.85, b2(1).YData, nanstd(nanmean(vali_all(:,senso_mo),2))/sqrt((length(nanmean(vali_all(:,senso_mo),2)))), 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1,'CapSize',0);
 text(0.77,-0.04,'ipsi','FontSize',11,'Color',[0.7 0.7 0.7]);
xpos = b2(2).XData + b2(2).XOffset;
errorbar(2.15, b2(2).YData, nanstd(nanmean(valc_all(:,senso_mo),2))/sqrt((length(nanmean(valc_all(:,senso_mo),2)))), 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1,'CapSize',0);
text(1,-0.04,'contra','FontSize',11,'Color',[0.3 0.3 0.3]);

xpos = b1(1).XData + b1(1).XOffset;
errorbar(0.85, b1(1).YData, nanstd(nanmean(vali_all(:,higher_a),2))/sqrt((length(nanmean(vali_all(:,higher_a),2)))), 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1,'CapSize',0);
 text(1.77,-0.04,'ipsi','FontSize',11,'Color',[0.7 0.7 0.7])
xpos = b1(2).XData + b1(2).XOffset;
errorbar(1.15, b1(2).YData, nanstd(nanmean(valc_all(:,higher_a),2))/sqrt((length(nanmean(valc_all(:,higher_a),2)))), 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1,'CapSize',0);
    text(2,-0.04,'contra','FontSize',11,'Color',[0.3 0.3 0.3]);
    % [ha pa] = adtest(nanmean(vali_all(:,senso_mo),2))
    % [hac pac] = adtest(nanmean(valc_all(:,senso_mo),2))
    % [ha pa] = adtest(nanmean(vali_all(:,higher_a),2))
    % [hac pac] = adtest(nanmean(valc_all(:,higher_a),2))
 disp('test for senso mo ipsi vs contra')
[u p1]=ttest(nanmean(vali_all(:,senso_mo),2),nanmean(valc_all(:,senso_mo),2))
nanmean(nanmean(vali_all(:,senso_mo),2))
nanstd(nanmean(vali_all(:,senso_mo),2))/sqrt((length(nanmean(vali_all(:,senso_mo),2))))
nanmean(nanmean(valc_all(:,senso_mo),2))
nanstd(nanmean(valc_all(:,senso_mo),2))/sqrt((length(nanmean(valc_all(:,senso_mo),2))))
text(1.9,1.05,'***','FontSize',18,'Color','k')
[u p1]=ttest(nanmean(vali_all(:,higher_a),2),nanmean(valc_all(:,higher_a),2))
text(0.95,1.05,'','FontSize',11,'Color','k')

   text(1.85,1.2,'v-a-sm','FontSize',11,'Color','r')
   text(0.85,1.2,'pf-m-l','FontSize',11,'Color',[0.5 0.5 0.5])
   ylabel('fILN')
box off;set(gca,'FontSize',11);set(gca,'TickDir','out')
h = gca;h.XAxis.Visible = 'off';
ylim([0 1.2])
cd(save_folder);saveas(gcf, 'ILN_large_groups_ipsivscontra.pdf');

%%  Panel C: Average in two major areas delta
clc;
frontal_idx=[1:8];lateral_idx=[9:16 44:45];somamo_idx=[17:26];visual_idx=[27:36];medial_idx=[37:39];aud_idx=[40:43];
idx=1;
temp_all={i_v1aafm c_v1aafm; i_s1aafm c_s1aafm; i_m1aafm c_m1aafm};
temp_color=[v1_color ;s1_color; m1_color];
senso_mo=[idx_modules{3} idx_modules{4} idx_modules{6}];
higher_a=[idx_modules{1} idx_modules{2} idx_modules{5}];
pan_title={'VISp','SSp-bfd','MOp'};

fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 450, 500, 220]);set(gcf,'color','w');t=tiledlayout('horizontal');
t.TileSpacing = 'compact';t.Padding = 'compact';  
for i=1:3
    nexttile
rm1=[];rm1=temp_all{i,1};rm2=[];rm2=temp_all{i,2};

[vali_all valdci_all val_i val_dci xerr yerr r_a p_a r_avg p_avg valc_all] = anatomy_correlation(rm1, rm2, idx);
disp('ipsi contra delta against each other')
% [ha pa] = adtest(nanmean(valdci_all(:,higher_a),2))
% [ha pa] = adtest(nanmean(valdci_all(:,senso_mo),2))
[u p1]=ttest(nanmean(valdci_all(:,higher_a),2),nanmean(valdci_all(:,senso_mo),2))
[ua0 p1a0]=ttest(nanmean(valdci_all(:,higher_a),2),0)

b1=bar(1,nanmean(nanmean(valdci_all(:,senso_mo),2)),0.5);b1.FaceColor='r';hold on;b1.EdgeColor='none';
b2=bar(2,nanmean(nanmean(valdci_all(:,higher_a),2)),0.5);b2.FaceColor=[0.9 0.9 0.9];set(b1,'ShowBaseLine','off');set(b2,'ShowBaseLine','off');b2.EdgeColor='none';
  r=1;rng(1);r1 = r-0.01 + (0.2)*rand(length(nanmean(valdci_all(:,senso_mo),2)),1);
sc1=scatter(r1,nanmean(valdci_all(:,senso_mo),2),15,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];

r=2;rng(1);r1 = r-0.01 + (0.2)*rand(length(nanmean(valdci_all(:,higher_a),2)),1);
sc1=scatter(r1,nanmean(valdci_all(:,higher_a),2),15,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
hold on;errorbar([1],[nanmean(nanmean(valdci_all(:,senso_mo),2))],[nanstd(nanmean(valdci_all(:,senso_mo),2))/sqrt(length(nanmean(valdci_all(:,senso_mo),2)))]...
     , 'LineStyle', 'none', ... 
          'Color', 'k', 'LineWidth', 1,'CapSize',0);hold on;set(b1,'ShowBaseLine','off');
    hold on;errorbar([2],[nanmean(nanmean(valdci_all(:,higher_a),2))],[nanstd(nanmean(valdci_all(:,higher_a),2))/sqrt(length(nanmean(valdci_all(:,higher_a),2)))]...
     , 'LineStyle', 'none', ... 
          'Color', 'k', 'LineWidth', 1,'CapSize',0);hold on;set(b1,'ShowBaseLine','off');

hold on;box off
xticks([1:2])
xlim([0 2.5])
   
     set(gca,'FontSize',11);set(gca,'TickDir','out');box off;
  
xticklabels({'v-a-sm','pf-m-l'})
xtickangle(45)
 ylim([-0.025 0.3]);
 ax=gca
ax.XTickLabel{1} = ['\color{red}' ax.XTickLabel{1}];
ax.XTickLabel{2} = ['\color{gray}' ax.XTickLabel{2}];
 ylabel('fILN contra - fILN ipsi');
 title(pan_title{i},'Color',temp_color(i,:),'FontWeight','normal','FontSize',11)
 if i==1
%h = gca;h.XAxis.Visible = 'off';
text(1.3,0.22,'**','FontSize',18);

 elseif i==2
h = gca;h.YAxis.Visible = 'off';
text(1.3,0.22,'***','FontSize',18);

 else i==3
  h = gca;h.YAxis.Visible = 'off';   
text(1.3,0.22,'**','FontSize',18);
 end

offsetAxes
end
cd(save_folder);saveas(gcf, 'ILN_delta_modules_2largegroups.pdf');

%% FIGURE 6: L6 dominates cortical hierarchy in sensory and motor areas within and between the two hemispheres
%Figure 6A: illustrator 
%% Data for flat maps shown in panel C
clc;
frontal_idx=[1:8];lateral_idx=[9:16 44:45];somamo_idx=[17:26];visual_idx=[27:36];medial_idx=[37:39];aud_idx=[40:43];
idx_modules={frontal_idx lateral_idx somamo_idx visual_idx medial_idx aud_idx};
dat_all1={};dat_all1={i_v1aafm i_s1aafm i_m1aafm};
dat_all2={};dat_all2={c_v1aafm c_s1aafm c_m1aafm};
senso_mo=[idx_modules{3} idx_modules{4} idx_modules{6}];
higher_a=[idx_modules{1} idx_modules{2} idx_modules{5}];
ymean1=[];ymean2=[];ymean3=[];
all_l23p=[];all_l5p=[];all_l6p=[];
for j=1:3
rm1=[];rm2=[];rm1=dat_all1{j};rm2=dat_all2{j};

l23c=squeeze(rm2(2,:,:));
l5c=squeeze(rm2(4,:,:));
l6c=squeeze(rm2(5,:,:));
l23i=squeeze(rm1(2,:,:));
l5i=squeeze(rm1(4,:,:));
l6i=squeeze(rm1(5,:,:));
ymean1(:,j)=nanmean(l23c-l23i,2);
ymean2(:,j)=nanmean(l5c-l5i,2);
ymean3(:,j)=nanmean(l6c-l6i,2);
p_sl23=[];p_sl5=[];p_sl6=[];hl23=[];hl5=[];hl6=[];
    for m=1:45
    [u p_sl23(m)]=ttest(l23c(m,:),l23i(m,:));
    [u p_sl5(m)]=ttest(l5c(m,:),l5i(m,:));
    [u p_sl6(m)]=ttest(l6c(m,:),l6i(m,:));
    try
    [hl23(m)] = adtest(l23i(m,:));
    [hl5(m)] = adtest(l5i(m,:));
    [hl6(m)] = adtest(l6i(m,:));
    catch
    end
 
    
    end
 all_l23p(:,j)=p_sl23;
 all_l5p(:,j)=p_sl5;
  all_l6p(:,j)=p_sl6;
end
%% save for MT flatmaps/table
pthresh1=[];pthresh2=[];pthresh3=[];
pthresh1=all_l23p<0.05;
pthresh2=all_l5p<0.05;
pthresh3=all_l6p<0.05;
cd(str)
delta_layers_cortex =table(cortex_names(:,2),ymean1(:,1),ymean2(:,1),ymean3(:,1),ymean1(:,2),ymean2(:,2),ymean3(:,2),ymean1(:,3),ymean2(:,3),ymean3(:,3),...
    pthresh1(:,1),pthresh2(:,1),pthresh3(:,1),pthresh1(:,2),pthresh2(:,2),pthresh3(:,2),pthresh1(:,3),pthresh2(:,3),pthresh3(:,3),...
    'variablenames',{'Area','V1dL2/3','V1dL5','V1dL6','S1dL2/3','S1dL5','S1dL6','M1dL2/3','M1dL5','M1dL6','V1L2/3sig','V1L5sig','V1L6sig'...
    'S1L2/3sig','S1L5sig','S1L6sig','M1L2/3sig','M1L5sig','M1L6sig'});
writetable(delta_layers_cortex,'delta_layers_cortex.csv')
%% Panel B: show fraction of layers L2/3, L5 and L6 only for the senso- motor module
clc;
frontal_idx=[1:8];lateral_idx=[9:16 44:45];somamo_idx=[17:26];visual_idx=[27:36];medial_idx=[37:39];aud_idx=[40:43];
dat_all1={};dat_all1={i_v1aafm i_s1aafm i_m1aafm};
dat_all2={};dat_all2={c_v1aafm c_s1aafm c_m1aafm};
idx_modules={frontal_idx lateral_idx somamo_idx visual_idx medial_idx aud_idx};
temp_color=[v1_color ;s1_color; m1_color];
senso_mo=[idx_modules{3} idx_modules{4} idx_modules{6}];
higher_a=[idx_modules{1} idx_modules{2} idx_modules{5}];
pan_title={'VISp','SSp-bfd','MOp'};
temp_color=[v1_color ;s1_color; m1_color];
temp_c1=[0.8 0.8 0.8];
temp_c2=[0.3 0.3 0.3];
allsm_l23i=[];allsm_l5i=[];allsm_l6i=[];
allsm_l23c=[];allsm_l5c=[];allsm_l6c=[];
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 100, 250, 600]);set(gcf,'color','w');t=tiledlayout('vertical');
t.TileSpacing = 'compact';t.Padding = 'compact';
for j=1:3
    
rm1=[];rm2=[];rm1=dat_all1{j};rm2=dat_all2{j};

l23smi=squeeze(nanmean(rm1(2,senso_mo,:),2));
l5smi=squeeze(nanmean(rm1(4,senso_mo,:),2));
l6smi=squeeze(nanmean(rm1(5,senso_mo,:),2));
l23smc=squeeze(nanmean(rm2(2,senso_mo,:),2));
l5smc=squeeze(nanmean(rm2(4,senso_mo,:),2));
l6smc=squeeze(nanmean(rm2(5,senso_mo,:),2));
allsm_l23i=[allsm_l23i; l23smi];
allsm_l5i=[allsm_l5i ;l5smi];
allsm_l6i=[allsm_l6i ;l6smi];
allsm_l23c=[allsm_l23c ;l23smc];
allsm_l5c=[allsm_l5c; l5smc];
allsm_l6c=[allsm_l6c; l6smc];
%ymean1=[];ymean1=nanmean(squeeze(rm1(2,idx_modules{i},:)));
disp(['difference ipsi contra ' pan_title{j}])
%[hac pac] = adtest(l23smi)
[uu kk]=ttest(l23smi,l23smc,'Tail','right')
[uu5 kk5]=ttest(l5smi,l5smc)
[uu6 kk6]=ttest(l6smi,l6smc,'Tail','left')
nexttile

b1=bar(1,[mean(l23smi); mean(l23smc)]);hold on;b1(1).FaceColor=temp_c1;b1(2).FaceColor=temp_c2;
set(b1,'ShowBaseLine','off');b1(1).EdgeColor='none';b1(2).EdgeColor='none';
b1=bar(2,[mean(l5smi); mean(l5smc)]);hold on;b1(1).FaceColor=temp_c1;b1(2).FaceColor=temp_c2;
set(b1,'ShowBaseLine','off');b1(1).EdgeColor='none';b1(2).EdgeColor='none';
b1=bar(3,[mean(l6smi); mean(l6smc)]);hold on;b1(1).FaceColor=temp_c1;b1(2).FaceColor=temp_c2;
set(b1,'ShowBaseLine','off');b1(1).EdgeColor='none';b1(2).EdgeColor='none';
r=1;
rng(1);r1 = r-0.2 + (0.1)*rand(length(l23smi),1);
rng(1);r2 = r +  (0.35)*rand(length(l23smc),1);
sc1=scatter([r1],[l23smi],5,'k','filled');hold on;
sc1(1).MarkerEdgeColor=[0.5 0.5 0.5];sc1(1).MarkerFaceColor=[0.2 0.2 0.2];sc1(1).MarkerFaceAlpha=[0.5];sc1(1).MarkerEdgeAlpha=[0.5];
sc1=scatter([r2],[l23smc],5,'k','filled');hold on;
sc1(1).MarkerEdgeColor=[0.5 0.5 0.5];sc1(1).MarkerFaceColor=[0.2 0.2 0.2];sc1(1).MarkerFaceAlpha=[0.5];sc1(1).MarkerEdgeAlpha=[0.5];
r=2;
rng(1);r1 = r-0.2 + (0.1)*rand(length(l5smi),1);
rng(1);r2 = r +  (0.35)*rand(length(l5smc),1);
sc1=scatter([r1],[l5smi],5,'k','filled');hold on;
sc1(1).MarkerEdgeColor=[0.5 0.5 0.5];sc1(1).MarkerFaceColor=[0.2 0.2 0.2];sc1(1).MarkerFaceAlpha=[0.5];sc1(1).MarkerEdgeAlpha=[0.5];
sc1=scatter([r2],[l5smc],5,'k','filled');hold on;
sc1(1).MarkerEdgeColor=[0.5 0.5 0.5];sc1(1).MarkerFaceColor=[0.2 0.2 0.2];sc1(1).MarkerFaceAlpha=[0.5];sc1(1).MarkerEdgeAlpha=[0.5];
r=3;
rng(1);r1 = r-0.2 + (0.1)*rand(length(l6smi),1);
rng(1);r2 = r +  (0.35)*rand(length(l6smc),1);
sc1=scatter([r1],[l6smi],5,'k','filled');hold on;
sc1(1).MarkerEdgeColor=[0.5 0.5 0.5];sc1(1).MarkerFaceColor=[0.2 0.2 0.2];sc1(1).MarkerFaceAlpha=[0.5];sc1(1).MarkerEdgeAlpha=[0.5];
sc1=scatter([r2],[l6smc],5,'k','filled');hold on;
sc1(1).MarkerEdgeColor=[0.5 0.5 0.5];sc1(1).MarkerFaceColor=[0.2 0.2 0.2];sc1(1).MarkerFaceAlpha=[0.5];sc1(1).MarkerEdgeAlpha=[0.5];
  hold on;h=errorbar([1-0.15],[mean(l23smi)],[nanstd(l23smi)/sqrt(length(l23smi))], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
  hold on;h=errorbar([1+0.15],[mean(l23smc)],[nanstd(l23smc)/sqrt(length(l23smc))], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
    hold on;h=errorbar([2-0.15],[mean(l5smi)],[nanstd(l5smi)/sqrt(length(l5smi))], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
  hold on;h=errorbar([2+0.15],[mean(l5smc)],[nanstd(l5smc)/sqrt(length(l5smc))], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
    hold on;h=errorbar([3-0.15],[mean(l6smi)],[nanstd(l6smi)/sqrt(length(l6smi))], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
  hold on;h=errorbar([3+0.15],[mean(l6smc)],[nanstd(l6smc)/sqrt(length(l6smc))], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
  box off;

 ylim([0 0.5]);
 ylabel('Fraction of neurons');
 xlim([0.5 3.5]);
 xticks([1:3]);xticklabels({'L2/3','L5','L6a'});set(gca,'FontSize',11);set(gca,'TickDir','out')
offsetAxes
 text(1.6,0.5,pan_title{j},'Color',temp_color(j,:),'FontWeight','normal');
 
 text(0.8,0.48,'***','Color','k','FontSize',18);
 if j==1 | j==2
 text(2.8,0.48,'***','Color','k','FontSize',18);
 end

end
cd(save_folder);saveas(gcf, 'fraction_layer_change_semo.pdf');

%% Panel D: Decrease increase percentage of sensory motor areas
clc;
pan_title={'VISp','SSp-bfd','MOp'};
temp_color=[v1_color ;s1_color; m1_color];
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 100, 230, 400]);set(gcf,'color','w');t=tiledlayout('vertical');
t.TileSpacing = 'compact';t.Padding = 'compact';
for j=1:3
nexttile
temp_don1=sum(ymean1(senso_mo(pthresh1(senso_mo,j)),j)>0)/24;
temp_don2=sum(ymean2(senso_mo(pthresh2(senso_mo,j)),j)>0)/24
temp_don3=sum(ymean3(senso_mo(pthresh3(senso_mo,j)),j)>0)/24;

temp_don4=sum(ymean1(senso_mo(pthresh1(senso_mo,j)),j)<0)/24;
temp_don5=sum(ymean2(senso_mo(pthresh2(senso_mo,j)),j)<0)/24
temp_don6=sum(ymean3(senso_mo(pthresh3(senso_mo,j)),j)<0)/24;

%donut([temp_don4,temp_don5,temp_don6;temp_don1,temp_don2,temp_don3],{''},{[ .945 .345 .329],[ .376 .741 .408],[ .365 .647 .855 ]})
%xlim([-4 4]);ylim([-4 4]);

b1=barh([temp_don3,temp_don2,temp_don1],0.5);b1.FaceColor=[246 136 68]/256;b1.EdgeColor='none';
hold on
b1=barh([-temp_don6,-temp_don5,-temp_don4],0.5);b1.FaceColor=[48 65 154]/256;b1.EdgeColor='none';
box off;
xlim([-0.6 0.6])
offsetAxes
  h = gca;h.YAxis.Visible = 'off'; 
  if j==1
      text(-0.4,4,'decrease','FontSize',11);
      text(0.05,4,'increase','FontSize',11);
       h = gca;h.XAxis.Visible = 'off'; 
      % text(-0.5,6,'contra - ipsi fraction','FontSize',11)
     %ylim([0 6])
  end
  if j==2
       h = gca;h.XAxis.Visible = 'off'; 
  end
  if j==3
     xlabel('Percentage areas (%)')
  end
   text(-0.78,3,'L2/3','FontSize',11);
      text(-0.75,2,'L5','FontSize',11);
      text(-0.75,1,'L6a','FontSize',11);
  set(gca,'FontSize',11);set(gca,'TickDir','out');
  xticks([-0.6:0.2:0.6]);xticklabels({'60','40','20','0','20','40','60'});
  title(pan_title{j},'Color',temp_color(j,:),'FontWeight','normal')
end
cd(save_folder);saveas(gcf, 'percentage_layer_change_ILN.pdf');
%% Remove L6 and then L5 for ILN ipsi and contra
frontal_idx=[1:8];lateral_idx=[9:16 44:45];somamo_idx=[17:26];visual_idx=[27:36];medial_idx=[37:39];aud_idx=[40:43];
idx_modules={frontal_idx lateral_idx somamo_idx visual_idx medial_idx aud_idx};
senso_mo=[idx_modules{3} idx_modules{4} idx_modules{6}];
higher_a=[idx_modules{1} idx_modules{2} idx_modules{5}];

 rm2=[];rm2=cat(3,c_v1aafm(:,senso_mo,:), c_s1aafm(:,senso_mo,:),c_m1aafm(:,senso_mo,:));
 rm1=[];rm1=cat(3,i_v1aafm(:,senso_mo,:), i_s1aafm(:,senso_mo,:),i_m1aafm(:,senso_mo,:));
 % rm2=[];rm2=c_m1aafm(:,senso_mo,:);
 % rm1=[];rm1=i_m1aafm(:,senso_mo,:);
idx=1;
[ilni ilnc ilni_e ilnc_e order_r order_rc all_i all_c] = iln_change(rm1,rm2,idx);
%exL6
idx=7;
[ilni_l6 ilnc_l6 ilni_el6 ilnc_el6 order_l6 order_l6c all_il6 all_cl6] = iln_change(rm1,rm2,idx);
%ex L5
idx=8;
[ilni_l5 ilnc_l5 ilni_el5 ilnc_el5 order_l5 order_l5c all_il5 all_cl5] = iln_change(rm1,rm2,idx);
%% 
for i=1:24
    [pi1(i) u2]=signrank(all_i(:,i),all_il5(:,i));
    [pi2(i) u2]=signrank(all_i(:,i),all_il6(:,i));
    [pi3(i) u2]=signrank(all_il5(:,i),all_il6(:,i));
    [pc1(i) u2]=signrank(all_c(:,i),all_cl5(:,i))
    [pc2(i) u2]=signrank(all_c(:,i),all_cl6(:,i))
    [pc3(i) u2]=signrank(all_cl5(:,i),all_cl6(:,i))
end
%0.05/24

%% Panel E / F: IPSI ILN excl L5 and L6
clc;
temp_c1=[0.8 0.8 0.8];
temp_c2=[0.3 0.3 0.3];


fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 300, 450, 300]);set(gcf,'color','w');
alpha1=1;
title('','FontWeight','normal','Color','k');
ylabel('fILN');
li=line([1 size(rm1,2)],[0.5 0.5],'Color',[0.5 0.5 0.5],'LineStyle',':');
 hold on;plot(ilni(order_r),'Color',[0.5 0.5 0.5]);
 hold on;plot(ilni_l5(order_l5),'Color',[0.5 0.5 0.5]);
 hold on;plot(ilni_l6(order_l6),'Color','k');
hold on;h=errorbar([1:size(rm1,2)],[ilni(order_r)],[ilni_e(order_r)]...
    , 'LineStyle', 'none', ... 
        'Color', [0.5 0.5 0.5], 'LineWidth', 0.5,'CapSize',0);set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha1]);
hold on;h=errorbar([1:size(rm1,2)],[ilni_l5(order_l5)],[ilni_el5(order_l5)]...
    , 'LineStyle', 'none', ... 
        'Color', [0.5 0.5 0.5], 'LineWidth', 0.5,'CapSize',0);
hold on;h=errorbar([1:size(rm1,2)],[ilni_l6(order_l6)],[ilni_el6(order_l6)]...
    , 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 0.5,'CapSize',0);
hold on;pp1=scatter(1:size(rm1,2),ilni(order_r),30,'o');pp1.MarkerEdgeColor=[0.5 0.5 0.5];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor=[0.5 0.5 0.5];
hold on;pp1=scatter(1:size(rm1,2),ilni_l5(order_l5),30,'o');pp1.MarkerEdgeColor=[0.5 0.5 0.5];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor='w';box off;
hold on;pp1=scatter(1:size(rm1,2),ilni_l6(order_l6),30,'o');pp1.MarkerEdgeColor=[0.5 0.5 0.5];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor='k';
xticks([1:size(rm1,2)]);ylim([0 1])
ylabel('Ipsi fILN');xlabel('Cortical areas')
set(gca,'FontSize',11);set(gca,'TickDir','out');
text(size(rm1,2)-2,0.35,'Default','Color',[0.4 0.4 0.4]);
text(size(rm1,2)-2,0.28,'Excl L5','Color',[0.7 0.7 0.7]);
text(size(rm1,2)-2,0.21,'Excl L6','Color','k');
offsetAxes;

     h = gca;h.XAxis.Visible = 'off'  
     text(size(rm1,2)/2-7,-0.1,'Cortical Areas (resorted by ipsi fILN)','FontSize',11);
     %cd(save_folder);saveas(gcf, 'ILNipsi_sorted_areasALL.pdf');
 % [p,tbl,stats]=anova1([ilni' ilni_l5' ilni_l6'])
 % presults = multcompare(stats)
cd(save_folder);saveas(gcf, 'ipsi_ILN_removeL5L6.pdf');

%% Contra ILN excl L5 and L6
temp_c1=[0.8 0.8 0.8];
temp_c2=[0.3 0.3 0.3];
alpha1=1;
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 300, 450, 300]);set(gcf,'color','w');
title('','FontWeight','normal','Color','k');
ylabel('fILN')
li=line([1 size(rm2,2)],[0.5 0.5],'Color',[0.5 0.5 0.5],'LineStyle',':');
 hold on;plot(ilnc(order_rc),'Color',[0.5 0.5 0.5]);
 hold on;plot(ilnc_l5(order_l5c),'Color',[0.5 0.5 0.5]);
 hold on;plot(ilnc_l6(order_l6c),'Color','k');
hold on;h=errorbar([1:size(rm2,2)],[ilnc(order_rc)],[ilnc_e(order_rc)]...
    , 'LineStyle', 'none', ... 
        'Color', [0.5 0.5 0.5], 'LineWidth', 0.5,'CapSize',0);set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha1]);
hold on;h=errorbar([1:size(rm2,2)],[ilnc_l5(order_l5c)],[ilnc_el5(order_l5c)]...
    , 'LineStyle', 'none', ... 
        'Color', [0.5 0.5 0.5], 'LineWidth', 0.5,'CapSize',0);
hold on;h=errorbar([1:size(rm2,2)],[ilnc_l6(order_l6c)],[ilnc_el6(order_l6c)]...
    , 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 0.5,'CapSize',0);

hold on;pp1=scatter(1:size(rm2,2),ilnc(order_rc),30,'o');pp1.MarkerEdgeColor=[0.5 0.5 0.5];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor=[0.5 0.5 0.5];
hold on;pp1=scatter(1:size(rm2,2),ilnc_l5(order_l5c),30,'o');pp1.MarkerEdgeColor=[0.5 0.5 0.5];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor='w';box off;
hold on;pp1=scatter(1:size(rm2,2),ilnc_l6(order_l6c),30,'o');pp1.MarkerEdgeColor=[0.5 0.5 0.5];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor='k';
xticks([1:size(rm2,2)]);ylim([0 1])
ylabel('Contra fILN');xlabel('Cortical areas')
set(gca,'FontSize',11);set(gca,'TickDir','out');
text(size(rm2,2)-2,0.35,'Default','Color',[0.5 0.5 0.5]);
text(size(rm2,2)-2,0.28,'Excl L5','Color',[0.7 0.7 0.7]);
text(size(rm2,2)-2,0.21,'Excl L6','Color','k');
offsetAxes;

     h = gca;h.XAxis.Visible = 'off'  
     text(size(rm2,2)/2-7,-0.1,'Cortical Areas (resorted by contra fILN)','FontSize',11);
     %cd(save_folder);saveas(gcf, 'ILNipsi_sorted_areasALL.pdf');
 % [p,tbl,stats]=anova1([ilni' ilni_l5' ilni_l6'])
 % presults = multcompare(stats)
  % [p,tbl,stats]=kruskalwallis([ilnc' ilnc_l5' ilnc_l6'])
  % presults = multcompare(stats)
cd(save_folder);saveas(gcf, 'contra_ILN_removeL5L6.pdf');

%% Panel E /F  : Comparison  MEAN ipsi and contra 
%ipsi
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 300, 650, 400]);set(gcf,'color','w');t=tiledlayout('horizontal');
t.TileSpacing = 'compact';t.Padding = 'compact';
nexttile
hold on;line([0 3.5],[nanmean(ilni) nanmean(ilni)],'LineStyle',':','Color','r','LineWidth',1);
for i=1:24
 pp1=plot([1 2 3],[ilni_l5(i) ilni(i) ilni_l6(i)] ,'-','Color',[0.5 0.5 0.5 0.3]);
 hold on;sc1=scatter(1,ilni_l5(i),40);sc1.MarkerFaceColor = [1 1 1];sc1.MarkerEdgeColor = [0.5 0.5 0.5];sc1.MarkerFaceAlpha=0.5;
 hold on;sc1=scatter(2,ilni(i),40);sc1.MarkerFaceColor = [0.5 0.5 0.5];sc1.MarkerEdgeColor = [1 1 1];sc1.MarkerFaceAlpha=0.5;
 hold on;sc1=scatter(3,ilni_l6(i),40);sc1.MarkerFaceColor = [0 0 0];sc1.MarkerEdgeColor = [0.5 0.5 0.5];sc1.MarkerFaceAlpha=0.5;
end

hold on;h=errorbar([0.8],[nanmean(ilni_l5)],[nanstd(ilni_l5)/sqrt(length(ilni_l5))]...
    , 'LineStyle', 'none', 'MarkerSize',8,... 
        'Color', [0.8 0.8 0.8], 'LineWidth', 1,'CapSize',0);
hold on;h=errorbar([1.8],[nanmean(ilni)],[nanstd(ilni)/sqrt(length(ilni))]...
    , 'LineStyle', 'none',  'MarkerSize',8,... 
        'Color', [0.5 0.5 0.5], 'LineWidth', 1,'CapSize',0);
hold on;h=errorbar([2.8],[nanmean(ilni_l6)],[nanstd(ilni_l6)/sqrt(length(ilni_l6))]...
    , 'LineStyle', 'none', 'MarkerSize',8,... 
        'Color', [0 0 0], 'LineWidth', 1,'CapSize',0);
hold on;sc1=scatter(0.8,nanmean(ilni_l5),60,'^');sc1.MarkerFaceColor = [1 1 1];sc1.MarkerEdgeColor = [0.5 0.5 0.5];
hold on;sc1=scatter(1.8,nanmean(ilni),60,'^');sc1.MarkerFaceColor = [0.5 0.5 0.5];sc1.MarkerEdgeColor = [1 1 1];
hold on;sc1=scatter(2.8,nanmean(ilni_l6),60,'^');sc1.MarkerFaceColor = [0 0 0];sc1.MarkerEdgeColor = [0.5 0.5 0.5];

xlim([0 4]);xticks([0:1:4]);yticks([0:0.2:1]);
box off;;xticklabels({'','Excl L5','Default','Excl L6',''});ylabel('ipsilateral fILN');set(gca,'FontSize',11);set(gca,'TickDir','out');
xtickangle(45);offsetAxes;ylim([0 1]);
%stats ipsi

%[hac pac] = adtest(ilni)

%contra

nexttile
hold on;line([0 3.5],[nanmean(ilni) nanmean(ilni)],'LineStyle',':','Color','r','LineWidth',1);
for i=1:24
 pp1=plot([1 2 3],[ilnc_l5(i) ilnc(i) ilnc_l6(i)] ,'-','Color',[0.5 0.5 0.5 0.3]);
 hold on;sc1=scatter(1,ilnc_l5(i),40);sc1.MarkerFaceColor = [1 1 1];sc1.MarkerEdgeColor = [0.5 0.5 0.5];sc1.MarkerFaceAlpha=0.5;
 hold on;sc1=scatter(2,ilnc(i),40);sc1.MarkerFaceColor = [0.5 0.5 0.5];sc1.MarkerEdgeColor = [1 1 1];sc1.MarkerFaceAlpha=0.5;
 hold on;sc1=scatter(3,ilnc_l6(i),40);sc1.MarkerFaceColor = [0 0 0];sc1.MarkerEdgeColor = [0.5 0.5 0.5];sc1.MarkerFaceAlpha=0.5;
end
hold on;h=errorbar([0.8],[nanmean(ilnc_l5)],[nanstd(ilnc_l5)/sqrt(length(ilnc_l5))]...
    , 'LineStyle', 'none', 'MarkerSize',8,... 
        'Color', [0.8 0.8 0.8], 'LineWidth', 1,'CapSize',0);
hold on;h=errorbar([1.8],[nanmean(ilnc)],[nanstd(ilnc)/sqrt(length(ilnc))]...
    , 'LineStyle', 'none',  'MarkerSize',8,... 
        'Color', [0.5 0.5 0.5], 'LineWidth', 1,'CapSize',0);
hold on;h=errorbar([2.8],[nanmean(ilnc_l6)],[nanstd(ilnc_l6)/sqrt(length(ilnc_l6))]...
    , 'LineStyle', 'none', 'MarkerSize',8,... 
        'Color', [0 0 0], 'LineWidth', 1,'CapSize',0);
hold on;sc1=scatter(0.8,nanmean(ilnc_l5),60,'^');sc1.MarkerFaceColor = [1 1 1];sc1.MarkerEdgeColor = [0.5 0.5 0.5];
hold on;sc1=scatter(1.8,nanmean(ilnc),60,'^');sc1.MarkerFaceColor = [0.5 0.5 0.5];sc1.MarkerEdgeColor = [1 1 1];
hold on;sc1=scatter(2.8,nanmean(ilnc_l6),60,'^');sc1.MarkerFaceColor = [0 0 0];sc1.MarkerEdgeColor = [0.5 0.5 0.5];
box off;xlim([0 4]);xticklabels({'','Excl L5','Default','Excl L6',''});ylabel('contralateral fILN');set(gca,'FontSize',11);set(gca,'TickDir','out');
xtickangle(45);offsetAxes;ylim([0 1]);yticks([0:0.2:1]);
cd(save_folder);saveas(gcf, 'line_paired_ILN_removeL5L6.pdf');


%% END


%% Supplementary
%% 
%% Figure 1–figure supplement 1:
%% Panel B Plot all cells ipis+contra  TOTAL numbers for supplement
temp_color=[v1_color ;s1_color; m1_color];
p1=[];p1=squeeze(nansum(nansum(i_v1aam)))+squeeze(nansum(nansum(c_v1aam)));
p2=[];p2=squeeze(nansum(nansum(i_s1aam)))+squeeze(nansum(nansum(c_s1aam)));
p3=[];p3=squeeze(nansum(nansum(i_m1aam)))+squeeze(nansum(nansum(c_m1aam)));
temp_p= [p1 p2 p3];
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 150, 270]);set(gcf,'color','w');tiledlayout("horizontal");
%title(['Total cell nr: ' num2str(nansum(nansum(temp_p)))],'FontWeight','normal')
%title(['T: 3.17*10^6' ],'FontWeight','normal','FontSize',5)
for i=1:3
    hold on;
    b1=bar(i,nanmean(temp_p(:,i)),0.7);hold on;b1.FaceColor=temp_color(i,:);  set(b1,'ShowBaseLine','off');b1.EdgeColor='none';
      r=i;rng(i);r1 = r-0.01 + (0.2)*rand(length(temp_p(:,i)),1);
sc1=scatter(r1,temp_p(:,i),12,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
    hold on;errorbar([i],[nanmean(temp_p(:,i))],[nanstd(temp_p(:,i))/sqrt(length(temp_p(:,i)))]...
    , 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);hold on;
end
xlim([0 3.5]);xticks([1:3]);hold on;box off;;ylabel('Total nr of cells (ipsi + contra)');
set(gca,'FontSize',11);set(gca,'TickDir','out');box off;xtickangle(45);%ylim([0 410000]);%xticklabels({'VISp','SSp-bfd','MOp'})
h = gca;h.XAxis.Visible = 'off'; 
t1=text(0.4,-50000,'VISp','FontSize',11);set(t1,'Rotation',45);
   t1=text(0.8,-62000,'SSp-bfd','FontSize',11);set(t1,'Rotation',45);
    t1=text(2.4,-50000,'MOp','FontSize',11);set(t1,'Rotation',45);
offsetAxes
cd(save_folder);saveas(gcf, 'Total_nr_cells.pdf');
[p,tbl,stats] = anova1([temp_p])
presults = multcompare(stats)

%% Figure 1–figure supplement 2:
%% Panel D Bolus volume vs cell numbers
% Bolus cell numbers for m2
px=[];px=nansum(nansum(i_m2aa))+nansum(nansum(c_m2aa));

fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 270, 270]);set(gcf,'color','w');
temp_color=[v1_color ;s1_color; m1_color];
p1=[];p1=squeeze(nansum(nansum(i_v1aam)))+squeeze(nansum(nansum(c_v1aam)));
p2=[];p2=squeeze(nansum(nansum(i_s1aam)))+squeeze(nansum(nansum(c_s1aam)));
p3=[];p3=squeeze(nansum(nansum(i_m1aam)))+squeeze(nansum(nansum(c_m1aam)));
sc1=scatter(cell2mat(inj_bol(:,1)),p1,50,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=temp_color(1,:);
sc1=scatter(cell2mat(inj_bol(:,2)),p2,50,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=temp_color(2,:);
sc1=scatter(cell2mat(inj_bol(:,3)),p3,50,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=temp_color(3,:);
sc2=scatter(cell2mat(inj_bol(1,4)),px,50,'k','filled');hold on;sc2.MarkerEdgeColor=[1 1 1];sc2.MarkerEdgeAlpha=0.5;sc2.MarkerFaceColor=[0.5 0.5 0.5];
hold on;box off;ylabel('Total nr of cells (ipsi + contra)');xlabel('Injection volume (mm^{3})')
set(gca,'FontSize',11);set(gca,'TickDir','out');
xlim([0.1 1]);ylim([0 420000])

t1=text(0.8,420000,'VISp','FontSize',11,'Color',temp_color(1,:));
t1=text(0.8,390000,'SSp-bfd','FontSize',11,'Color',temp_color(2,:));
t1=text(0.8,360000,'MOp','FontSize',11,'Color',temp_color(3,:));
t1=text(0.8,330000,'Excl MOp','FontSize',11,'Color',[0.5 0.5 0.5]);
offsetAxes
cd(save_folder);saveas(gcf, 'bolusvscellnumber.pdf');

%% Panel B and C: Bar plot spill over
temp_color=[v1_color ;s1_color; m1_color; [0.5 0.5 0.5]];
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 800, 270]);set(gcf,'color','w');t=tiledlayout("horizontal");
t.TileSpacing = 'compact';
t.Padding = 'compact';
p1=[];p1={cell2mat(inj_perc(:,12:15)) cell2mat(inj_perc(:,7:11)) [cell2mat(inj_perc(:,1:6)) [0; 0; 0; 0; 0; 0]] cell2mat(inj_perc(:,16:17))};
areas_com={{'VISp', 'WM', 'VISl','OpticRad'};{'SSp-bfd','SSp-un','SSp-n','WM','OpticRad'};{'MOp','SSp-ul','SSp-ll','SSs','Cingb','MOs','WM'};{'MOp','MOs'}};
for i=1:4
    nexttile
b1=bar([1:size(p1{i},2)],[nanmean(p1{i})],0.8);b1.FaceColor=temp_color(i,:);set(b1,'ShowBaseLine','off');b1.EdgeColor='none';
 hold on;errorbar([1:size(p1{i},2)],[nanmean(p1{i})],[nanstd(p1{i})/sqrt(size(p1{i},1))]...
    , 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1,'CapSize',0);hold on;
ylabel('Percentage of viral bolus (%)');
xticklabels(areas_com{i}(:));xtickangle(45);
box off;set(gca,'FontSize',10);set(gca,'TickDir','out');
xlim([0.5 7.5]);ylim([0 100])
if i==2 | i==3 | i==4
        %h = gca;h.YAxis.Visible = 'off'  ;
end
offsetAxes
end
cd(save_folder);saveas(gcf, 'Percentage bolus.pdf');




%% Figure 4–figure supplement 1:  hierarhcy HVA
%
%% Panel D
dat_all1={};dat_all1={i_v1aafm i_s1aafm i_m1aafm};
dat_all2={};dat_all2={c_v1aafm c_s1aafm c_m1aafm};
%idx_modules={frontal_idx lateral_idx somamo_idx visual_idx medial_idx aud_idx};
idx_modules={somamo_idx visual_idx aud_idx frontal_idx lateral_idx medial_idx };
%module_names={'Pfron','Lat','SoMo','Vis','Med','Aud'};
module_names={'SoMo','Vis','Aud','Pfron','Lat','Med',};
pan_title={'VISp','SSp-bfd','MOp'};
temp_color=[v1_color ;s1_color; m1_color];
temp_c1=[0.8 0.8 0.8];
temp_c2=[0.3 0.3 0.3];
ventral_areas=[];dorsal_areas=[];
ventral_areas=[27 28 29 30];
dorsal_areas=[32 33 34 35 36];

p_fig3f=[];
%plot
%fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 250, 330]);set(gcf,'color','w');t=tiledlayout('vertical');
%t.TileSpacing = 'compact';t.Padding = 'compact';

rm1=[];rm2=[];rm1=dat_all1{1};rm2=dat_all2{1};
%calculate and plot 
i_metric=[];c_metric=[];
%calculate indexes across areas
for i=1:length(rm1)
temp_l=[];temp_l2=[];
temp_l=squeeze(rm1(:,i,:));
temp_l2=squeeze(rm2(:,i,:));
[i_metric(:,i,:)] = anatomy_indexcalc(temp_l);
[c_metric(:,i,:)] = anatomy_indexcalc(temp_l2);
end
visp_iln_i=i_metric(:,:,1);visp_hindex_i=i_metric(:,:,4);
visp_iln_c=c_metric(:,:,1);
visp_iln_i_l6=i_metric(:,:,7);
visp_iln_i_l5=i_metric(:,:,8);
visp_iln_i_l56=visp_iln_i-((visp_iln_i-visp_iln_i_l6)+(visp_iln_i-visp_iln_i_l5));

vareas_names=cortex_names(27:36,2)
kk1=[];kk2=[];
[ventral_iln kk1]=sort(nanmean(visp_iln_i(:,ventral_areas)),'descend');
ventral_ilnl6=nanmean(visp_iln_i_l6(:,ventral_areas));
ventral_ilnl56=nanmean(visp_iln_i_l56(:,ventral_areas));
ventral_ilnc=nanmean(visp_iln_c(:,ventral_areas));

v_ilnl6=ventral_ilnl6(kk1);
v_ilnl56=ventral_ilnl56(kk1);
v_ilnc=ventral_ilnc(kk1);
%yerr3=[];yerr3=nanstd(visp_iln_i_l6(:,ventral_areas))./sqrt(size(visp_iln_i_l6(:,ventral_areas),1));
%yerr3=[yerr3 NaN];
yerr5=[];yerr5=nanstd(visp_iln_i_l56(:,ventral_areas))./sqrt(size(visp_iln_i_l56(:,ventral_areas),1));
yerr5=[yerr5 NaN];
%yerr4=[];yerr4=nanstd(visp_iln_i_l6(:,dorsal_areas))./sqrt(size(visp_iln_i_l6(:,dorsal_areas),1));
yerr6=[];yerr6=nanstd(visp_iln_i_l56(:,dorsal_areas))./sqrt(size(visp_iln_i_l56(:,dorsal_areas),1));
yerr7=[];yerr7=nanstd(visp_iln_c(:,ventral_areas))./sqrt(size(visp_iln_c(:,ventral_areas),1));
yerr7=[yerr7 NaN];
kk1=[kk1 5];
[dorsal_iln kk2]=sort(nanmean(visp_iln_i(:,dorsal_areas)),'ascend');
dorsal_ilnl6=nanmean(visp_iln_i_l6(:,dorsal_areas));
dorsal_ilnl56=nanmean(visp_iln_i_l56(:,dorsal_areas));
dorsal_ilnc=nanmean(visp_iln_c(:,dorsal_areas));
yerr1=[];yerr1=nanstd(visp_iln_i(:,ventral_areas))./sqrt(size(visp_iln_i(:,ventral_areas),1));
yerr1=[yerr1 NaN];

yerr8=[];yerr8=nanstd(visp_iln_c(:,dorsal_areas))./sqrt(size(visp_iln_c(:,dorsal_areas),1));
%yerr8=[yerr8 NaN];

yerr2=[];yerr2=nanstd(visp_iln_i(:,dorsal_areas))./sqrt(size(visp_iln_i(:,dorsal_areas),1));
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 300, 330]);set(gcf,'color','w')
 hold on;h=errorbar([1:10],[ventral_iln NaN dorsal_iln],[yerr1(kk1) yerr2(kk2)], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
  %hold on;h=errorbar([1:10],[v_ilnl6 NaN dorsal_ilnl6(kk2)],[yerr3(kk1) yerr4(kk2)], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
  %hold on;h=errorbar([1:10],[v_ilnl56 NaN dorsal_ilnl56(kk2)],[yerr5(kk1) yerr6(kk2)], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
hold on;h=errorbar([1:10],[v_ilnc NaN dorsal_ilnc(kk2)],[yerr7(kk1) yerr8(kk2)], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);

  hold on;sc1=scatter([1:10],[ventral_iln NaN dorsal_iln] ,60,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=1;sc1.MarkerFaceColor=[0.75 0.75 0.75];
%hold on;sc2=scatter([1:10],[v_ilnl6 NaN dorsal_ilnl6(kk2)] ,60,'k','filled');hold on;sc2.MarkerEdgeColor=[1 1 1];sc2.MarkerEdgeAlpha=1;sc2.MarkerFaceColor=[0 0 0];
%hold on;sc2=scatter([1:10],[v_ilnl56 NaN dorsal_ilnl56(kk2)] ,60,'k','filled');hold on;sc2.MarkerEdgeColor=[1 1 1];sc2.MarkerEdgeAlpha=1;sc2.MarkerFaceColor=[0.4 0.4 0.4];
hold on;sc2=scatter([1:10],[v_ilnc NaN dorsal_ilnc(kk2)] ,60,'k','filled');hold on;sc2.MarkerEdgeColor=[1 1 1];sc2.MarkerEdgeAlpha=1;sc2.MarkerFaceColor=[0.3 0.3 0.3];
xticks([1:1:10]);xticklabels(vareas_names([kk1 kk2+5]));
ylabel('fILN ipsi');set(gca,'FontSize',11);set(gca,'TickDir','out');box off;%xtickangle(45);
line([5 5],[0.4 0.9],'LineStyle','--','Color','k');
text(2,0.9,'ventral','FontSize',11);
text(7,0.9,'dorsal','FontSize',11);
ylabel('fILN');text(7,0.4,'ipsi','FontSize',11,'Color',[0.75 0.75 0.75]);text(7,0.45,'contra','FontSize',11,'Color',[0.3 0.3 0.3]);
offsetAxes;
cd(save_folder);saveas(gcf, 'ILN_HVA_sorted.pdf');

%% Panel C: heatmap
dat_all1={};dat_all1={i_v1aafm i_s1aafm i_m1aafm};
rm1=[];rm2=[];rm1=dat_all1{1};rm2=dat_all2{1};
tm_hva=nanmean(rm1(2:6,27:36,:),3);

fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 1000, 350, 250]);set(gcf,'color','w');
nexttile
h=imagesc(tm_hva(:,[kk1 kk2+5]));
set(h, 'AlphaData', ~isnan(tm_hva(:,[kk1 kk2+5])))
cmap(temp_color(1,:),100,5,5);a1=colorbar;
ylabel(a1,'Fraction of neuron','FontSize',11,'Rotation',270);
[hh oo]=find(isnan(tm_hva(:,[kk1 kk2+5])));
coord_b=[oo hh];
for u=1:length(coord_b)
    hold on
    text(coord_b(u,1)-0.25,coord_b(u,2)-0.1,'x','FontSize',14,'Color','k')
end
box off;set(gca,'TickDir','out');
hold on;yticklabels({'L2/3','L4','L5','L6a','L6b'});
xticks([1:1:10]);xticklabels(vareas_names([kk1 kk2+5]));
%cd(save_folder);saveas(gcf, 'ILNHVA_heatmap.pdf');
%offsetAxes;


%% Figure 1–figure supplement 4
%use gad_perc
%plot
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 500, 200]);set(gcf,'color','w');
nexttile
for i=1:4
    gad_mean=[];gad_sem=[];gad_mean=gad_perc{i,1};gad_sem=gad_perc{i,2};
    b1=bar(i,gad_mean,0.8);hold on;b1(1).FaceColor=[1 1 1];b1(1).EdgeColor=[0 0 0];set(b1,'ShowBaseLine','off')
    hold on;h=errorbar([i],[gad_mean],[gad_sem], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
end

xticks([1:4]);xticklabels({'L2/3','L4','L5','L6'});ylabel('Percentage GAD+ (%)')
box off;set(gca,'FontSize',11);set(gca,'TickDir','out');
offsetAxes;ylim([0 4]);

nexttile
 gad_mean=[];gad_sem=[];gad_mean=0;gad_sem=0;
    b1=bar(1,gad_mean,0.8);hold on;b1(1).FaceColor=[1 1 1];b1(1).EdgeColor=[0 0 0];set(b1,'ShowBaseLine','off')
    hold on;h=errorbar([1],[gad_mean],[gad_sem], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
    xticks([0:1]);xticklabels({'','L6'});ylabel('Percentage NTSR1+ (%)');ylim([0 1])
box off;set(gca,'FontSize',11);set(gca,'TickDir','out');offsetAxes;ylim([0 4]);
offsetAxes;

%% REVIEWS
%1)Plot total number of cells vs L6 fraction per area 

%Plot total number of cells vs L6 fraction per area 

dat_all1={};dat_all1={i_v1aafm i_s1aafm i_m1aafm;c_v1aafm c_s1aafm c_m1aafm};
dat_all3={};dat_all3={v1a_tiam s1a_tiam m1a_tiam;v1a_tcam s1a_tcam m1a_tcam};
idx=1;
 module_names={'VISp','SSp-bfd','MOp'};
temp_color=[v1_color ;s1_color; m1_color];

for m=1:2
   fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 1200, 300]);set(gcf,'color','w'); tiledlayout("horizontal");
for j=1:3
    nexttile
    temp_c=temp_color(j,:);
rm1=[];rm1=dat_all1{m,j};
rm3=[];rm3=dat_all3{m,j};

l6i=squeeze(rm1(5,:,:));
hmetric=[];temp_metric1=[];
%calculate indexes across areas
for i=1:length(rm1)
temp_l=[];
temp_l=squeeze(rm1(:,i,:));
[hmetric(:,i,:)] = anatomy_indexcalc(temp_l);
end
temp_metric1=nanmean(hmetric(:,:,idx));


% xerr=nanstd(rm3,[],2)/sqrt(size(rm3,2));
% yerr=nanstd(l6i,[],2)/sqrt(size(rm3,2));
% hold on;errorbar(nanmean(rm3,2),nanmean(l6i,2),yerr,yerr,xerr,xerr, 'LineStyle', 'none', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5,'CapSize',0);hold on;

hold on;
sc1=scatter(nanmean(rm3,2),nanmean(l6i,2),50,temp_metric1,'filled');sc1.MarkerFaceColor=temp_c;sc1.MarkerEdgeColor='w';
%sc1=scatter3(nanmean(rm3,2),nanmean(l6i,2),temp_metric1,50,'filled');sc1.MarkerFaceColor=temp_c;sc1.MarkerEdgeColor='w';
%view(45,45)
 % a1=colorbar;
 % ylabel(a1,'fILN','FontSize',11,'Rotation',270);
r_a=[];p_a=[];
for i=1:size(rm3,2)
r1=[];p1=[];r2=[];p2=[];
[r1 p1] = corr(rm3(:,i),l6i(:,i),'Type','Spearman','Rows','complete');
[r2 p2] = corr(rm3(:,i),hmetric(i,:,idx)','Type','Spearman','Rows','complete');
[r3 p3] = corr(l6i(:,i),hmetric(i,:,idx)','Type','Spearman','Rows','complete');
r_a1(i)=r1(1);
p_a1(i)=p1(1);
r_a2(i)=r2(1);
p_a2(i)=p2(1);
r_a3(i)=r3(1);
p_a3(i)=p3(1);
end
p_a
combined_p = harmonic_mean_pvalue(p_a)
hold on;%text(1000,1,...
   % ['r= ' num2str(round(nanmean(r_a),2)) ' +- ' num2str(round(nanstd(r_a)/sqrt(length(r_a)),2))],'FontSize',11,'Color',temp_c);
set(gca, 'XScale', 'log');
ylabel('L6 fraction of cells');xlabel('Total number of cells per area (log)');ylim([-0.05 1]);zlabel('fILN');
 box off;set(gca,'FontSize',11);set(gca,'TickDir','out');
 %grid on;
 % if m==1
 % title('Ipsilateral');
 % else m==2
 %  title('Contralateral');
 % end
 if m==1
 title(['ipsi ' module_names{j}],'Color',temp_color(j,:),'FontWeight','normal');
 else m==2
     title(['contra ' module_names{j}],'Color',temp_color(j,:),'FontWeight','normal');
 end
end
%axis square

 %offsetAxes;
 % if m==1
 %     title('Ipsilateral');
 % else
 %     title('Contralateral');
 % end
end
%% 



%% decrease per area 
 module_names={'VISp','SSp-bfd','MOp'};
temp_color=[v1_color ;s1_color; m1_color];
dat_all1={};dat_all1={i_v1aam i_s1aam i_m1aam};
dat_all2={};dat_all2={c_v1aam c_s1aam c_m1aam};
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 800, 300]);set(gcf,'color','w');tiledlayout("horizontal")
for j=1:3
    ipsij=squeeze(sum(dat_all1{j}));
    contraj=squeeze(sum(dat_all2{j}));
    for k=1:size(ipsij,2)
        ipa(:,k)=ipsij(:,k)./(ipsij(:,k)+contraj(:,k));
        cpa(:,k)=contraj(:,k)./(ipsij(:,k)+contraj(:,k));
    end
    tempi=nanmean(ipa,2);
    tempc=nanmean(cpa,2);
    tempis=nanstd(ipa,[],2)./sqrt(6);
    tempcs=nanstd(cpa,[],2)./sqrt(6);
    %lines between paired data points
    nexttile
for i=1:length(tempi)
     pl=plot([1,2],[tempi(i),tempc(i)],'color',[0.5 0.5 0.5]); hold on;   
     hold on;er1=errorbar(1,tempi(i),tempis(i),'ok','MarkerFaceColor',temp_color(j,:),'Markersize',8);
     hold on;er1=errorbar(2,tempc(i),tempcs(i),'ok','MarkerFaceColor',temp_color(j,:),'Markersize',8);
end


xlim([0 3]);ylim([0 1])
ylabel('Fraction of cells in each cortical area');xticks([1:1:2]);xticklabels({'ipsi','contra'})
 box off;set(gca,'FontSize',11);set(gca,'TickDir','out');
 offsetAxes;
title(module_names{j},'Color',temp_color(j,:),'FontWeight','normal');
end


%% 
lc=[2 4 5];

anim_incl={1:6 7:12 13:18};
labels_b={'L2/3','L5','L6a'};temp_color=[v1_color ;s1_color; m1_color];panel_tit={'VISp','SSp-bfd','MOp'};
 module_names={'VISp','SSp-bfd','MOp'};
temp_color=[v1_color ;s1_color; m1_color];
dat_all1={};dat_all1={i_v1aam i_s1aam i_m1aam};
dat_all2={};dat_all2={c_v1aam c_s1aam c_m1aam};
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 800, 300]);set(gcf,'color','w');tiledlayout("horizontal")
for m=1:3
    nexttile
for j=1:3
    ipsij=[];contraj=[];ipsir=[];contrar=[];
    ipsij=dat_all1{j};
    contraj=dat_all2{j};
    for i=1:6
    ipsir(:,i)=nansum(ipsij(lc(m),:,i))/(nansum(ipsij(lc(m),:,i))+nansum(contraj(lc(m),:,i)));
    contrar(:,i)=nansum(contraj(lc(m),:,i))/(nansum(ipsij(lc(m),:,i))+nansum(contraj(lc(m),:,i)));
    end
     tempi=[];tempc=[];tempis=[];tempcs=[];
    tempi=nanmean(ipsir);
    tempc=nanmean(contrar);
    tempis=nanstd(ipsir)./sqrt(6);
    tempcs=nanstd(contrar)./sqrt(6);

     hold on

     pl=plot([1,2],[tempi,tempc],'color',[0.5 0.5 0.5]); hold on;   
     hold on;er1=errorbar(1,tempi,tempis,'ok','MarkerFaceColor',temp_color(j,:),'Markersize',8);
     hold on;er1=errorbar(2,tempc,tempcs,'ok','MarkerFaceColor',temp_color(j,:),'Markersize',8);
     xlim([0 3]);ylim([0 1])
ylabel('Fraction of cells');xticks([1:1:2]);xticklabels({'ipsi','contra'})
 box off;set(gca,'FontSize',11);set(gca,'TickDir','out');
 offsetAxes;
title(labels_b{m},'Color','k','FontWeight','normal');

end

end
text(2.3,0.3,module_names{1},'Color',temp_color(1,:));hold on;
text(2.3,0.2,module_names{2},'Color',temp_color(2,:));hold on;
text(2.3,0.1,module_names{3},'Color',temp_color(3,:));

%% Preservation index 
frontal_idx=[1:8];lateral_idx=[9:16 44:45];somamo_idx=[17:26];visual_idx=[27:36];medial_idx=[37:39];aud_idx=[40:43];
dat_all1={};dat_all1={i_v1aam i_s1aam i_m1aam};
dat_all2={};dat_all2={c_v1aam c_s1aam c_m1aam};
idx_modules={frontal_idx lateral_idx somamo_idx visual_idx medial_idx aud_idx};
temp_color=[v1_color ;s1_color; m1_color];
senso_mo=[idx_modules{3} idx_modules{4} idx_modules{6}];
higher_a=[idx_modules{1} idx_modules{2} idx_modules{5}];
pan_title={'VISp','SSp-bfd','MOp'};
temp_color=[v1_color ;s1_color; m1_color];
temp_c1=[0.8 0.8 0.8];
temp_c2=[0.3 0.3 0.3];
%fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 800, 300]);set(gcf,'color','w');tiledlayout("horizontal")
for j=1:3
  %nexttile; 
  da23=[];da5=[];da6=[];
rm1=[];rm2=[];rm1=dat_all1{j};rm2=dat_all2{j};
da23=(squeeze(rm1(2,senso_mo,:))-squeeze(rm2(2,senso_mo,:)))./(squeeze(rm1(2,senso_mo,:))+squeeze(rm2(2,senso_mo,:)));
da5=(squeeze(rm1(4,senso_mo,:))-squeeze(rm2(4,senso_mo,:)))./(squeeze(rm1(4,senso_mo,:))+squeeze(rm2(4,senso_mo,:)));
da6=(squeeze(rm1(5,senso_mo,:))-squeeze(rm2(5,senso_mo,:)))./(squeeze(rm1(5,senso_mo,:))+squeeze(rm2(5,senso_mo,:)));
da23m=1-nanmean(da23,2);
da5m=1-nanmean(da5,2);
da6m=1-nanmean(da6,2);
[rr pp]=adtest(da23m)

% imagesc([da23m da5m da6m]);
% colorbar

% l23smi=squeeze(nanmean(rm1(2,senso_mo,:),2));
% l5smi=squeeze(nanmean(rm1(4,senso_mo,:),2));
% l6smi=squeeze(nanmean(rm1(5,senso_mo,:),2));
% l23smc=squeeze(nanmean(rm2(2,senso_mo,:),2));
% l5smc=squeeze(nanmean(rm2(4,senso_mo,:),2));
% l6smc=squeeze(nanmean(rm2(5,senso_mo,:),2));
g1=[];g2=[];g3=[];
p1=[];p2=[];p3=[];
p1=da23m;p2=da5m;p3=da6m;
par=[];par=[p1; p2; p3];
g1=[];g1=ones(1,length(p1));
g2=[];g2=ones(1,length(p2))*2;
g3=[];g3=ones(1,length(p3))*3;
gro=[];gro=[g1 g2 g3]';
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 200, 300]);set(gcf,'color','w');
%hold on;line([0 3.5],[nanmedian(p1) nanmedian(p1)],'LineStyle',':','Color','r','LineWidth',1);
%hold on;line([0.8 6.5],[1 1],'LineStyle',':','Color','k','LineWidth',1)
violins = violinplot(par, gro,'ViolinColor',[[0.5 0.5 0.5];[0.8 0.8 0.8];[0.2 0.2 0.2]],'ShowMean', false,'ShowMedian', true,'ViolinAlpha',0,'MedianColor',[0 0 0],'BoxColor',[0 0 0]);box off;
xlim([0.5 3.5]);ylabel('Preservation index');set(gca,'FontSize',12);
% h = gca;h.XAxis.Visible = 'off';
% set(gca,'xtick',[]);
xticklabels({'L2/3','L5','L6'});
%ylim([-1.2 1.2]);
ylim([-0.2 1.2]);
set(gca,'TickDir','out'); 
title(pan_title{j},'Color',temp_color(j,:));hold on;
 if j==1
     line([1 2],[0.9 0.9],'Color','k');
     text(1.5,0.92,'*','FontSize',18);hold on;
     line([1 3],[1 1],'Color','k');
     text(1.8,1.02,'*','FontSize',18);
 elseif j==2
       line([1 3],[1.08 1.08],'Color','k');
     text(1.8,1.1,'***','FontSize',18);
 else j==3
 end
[p,tbl,stats] = kruskalwallis(par, gro);
 presults = multcompare(stats)

end
%[dchange] = anatomy_indexhemi(layers_i,layers_c)










%% ATTACHED FUNCTIONS
%% 1 cell selecter 
function filter_out = cell_selecter(data,  varargin)

p = inputParser();
p.KeepUnmatched = false;
p.CaseSensitive = false;
p.StructExpand  = false;
p.PartialMatching = true;

% filters based on data structure
addParameter(p, 'label', nan); 
addParameter(p, 'area' , nan) % when the filter is set to nan it is not used 
addParameter(p, 'wc', nan);
addParameter(p, 'sol', nan); 
addParameter(p, 'amp', nan);
addParameter(p, 'pair', nan);
addParameter(p, 'distx', nan);
addParameter(p, 'disty', nan);
addParameter(p, 'geno', nan);
addParameter(p, 'drugs', nan);
addParameter(p, 'optovariant', nan);
addParameter(p, 'layer', nan);
addParameter(p, 'qualityinput', nan);
addParameter(p, 'virustype', nan);


parse(p, varargin{:})
cell_filters = p.Results;

% remove fields with nan
fn = fieldnames(cell_filters);
for k=1:numel(fn)
    if isnumeric(cell_filters.(fn{k})) && isnan(cell_filters.(fn{k}))
        cell_filters = rmfield(cell_filters,fn{k});
    end
end

data2 = data;
% loop through filter and create filter vectors
fn = fieldnames(cell_filters);
filters_array = [];
for k = 1:length(fn)
    if strcmp(fn{k},'amp')
        temp_filter = cellfun(@(x) all(cell_filters.amp==x),{data2.(fn{k})});
        filters_array(k,:) = temp_filter;
    elseif strcmp(fn{k},'distx')
        temp_filter = cell_filters.transduction_QC <= [data2.transduction_QC];
        filters_array(k,:) = temp_filter;
    elseif strcmp(fn{k},'disty')
        temp_filter = cellfun(@(x) ~isempty(x),{data2.scracm_red});
        filters_array(k,:) = temp_filter;
    elseif strcmp(fn{k},'label')
        temp_filter = cellfun(@(x) all(cell_filters.label==x),{data2.(fn{k})});
        filters_array(k,:) = temp_filter;
   
     elseif strcmp(fn{k},'sol')
       temp_filter = cellfun(@(x) all(cell_filters.sol==x),{data2.(fn{k})});
        filters_array(k,:) = temp_filter;
         elseif strcmp(fn{k},'geno')
       temp_filter = cellfun(@(x) all(cell_filters.geno==x),{data2.(fn{k})});
        filters_array(k,:) = temp_filter;
           elseif strcmp(fn{k},'area')
       temp_filter = cellfun(@(x) all(cell_filters.area==x),{data2.(fn{k})});
        filters_array(k,:) = temp_filter;
        elseif strcmp(fn{k},'wc')
       temp_filter = cellfun(@(x) all(cell_filters.wc==x),{data2.(fn{k})});
        filters_array(k,:) = temp_filter;
    
          elseif strcmp(fn{k},'pair')
       temp_filter = cellfun(@(x) all(cell_filters.pair==x),{data2.(fn{k})});
        filters_array(k,:) = temp_filter;
              elseif strcmp(fn{k},'drugs')
       temp_filter = cellfun(@(x) all(cell_filters.drugs==x),{data2.(fn{k})});
        filters_array(k,:) = temp_filter;
            elseif strcmp(fn{k},'optovariant')
       temp_filter = cellfun(@(x) all(cell_filters.optovariant==x),{data2.(fn{k})});
        filters_array(k,:) = temp_filter;
          elseif strcmp(fn{k},'layer')
       temp_filter = cellfun(@(x) all(cell_filters.layer==x),{data2.(fn{k})});
        filters_array(k,:) = temp_filter;
           elseif strcmp(fn{k},'qualityinput')
       temp_filter = cellfun(@(x) all(cell_filters.qualityinput==x),{data2.(fn{k})});
        filters_array(k,:) = temp_filter;
             elseif strcmp(fn{k},'virustype')
       temp_filter = cellfun(@(x) all(cell_filters.virustype==x),{data2.(fn{k})});
        filters_array(k,:) = temp_filter;
    else
        temp_filter = cellfun(@(x) all(cell_filters.(fn{k})==x),{data2.(fn{k})});
        filters_array(k,:) = temp_filter;
    end
end
if isempty(filters_array)
    filter_out = ones(1,length(data2));
else
    filter_out = all(filters_array,1);
end
disp([num2str(sum(filter_out)) ' cells/areas/animals pass all filters' ])
end


%% 2 function to read out cell numbers and fractions 
function [i_animal c_animal i_areas_animal c_areas_animal ...
    i_areas_animalf c_areas_animalf i_areas_animalcd c_areas_animalcd] = anatomy_cellnr(data,type,cortex_names)
aa=[18 25 31 41];
for i=1:length(type)
  ipsi_all=[];contra_all=[];ipsi_allcd=[];contra_allcd=[];
  %k=areas
  for k=1:length(cortex_names)
  cortex_idx=[];cortex_idx=find(contains(string(data(type(i)).area_names),cortex_names{k,3})==1);
  a_names=[];a_names=data(type(i)).area_names(cortex_idx);
  %sort from L1 to L6b in a given area
  rb=[];sort_v=[];[rb sort_v]=sort(a_names);
  %save cell numbers with the sorted indx for ipsi and contra
  ipsi_nr=[];ipsi_nr=data(type(i)).ipsi_cellnr(cortex_idx(sort_v));
  contra_nr=[];contra_nr=data(type(i)).contra_cellnr(cortex_idx(sort_v));
  ipsi_cd=[];ipsi_cd=data(type(i)).ipsi_celldens(cortex_idx(sort_v));
  contra_cd=[];contra_cd=data(type(i)).contra_celldens(cortex_idx(sort_v));
  %distinguish areas that dont have a L4 (5 layers only) from the rest
  %five layers
  if length(sort_v)==5
      ipsi_nr=[ipsi_nr(1:2); NaN ;ipsi_nr(3:end)];
      contra_nr=[contra_nr(1:2); NaN ;contra_nr(3:end)];
      ipsi_cd=[ipsi_cd(1:2); NaN ;ipsi_cd(3:end)];
      contra_cd=[contra_cd(1:2); NaN ;contra_cd(3:end)];
  %six layers
  else
      ipsi_nr=ipsi_nr;contra_nr=contra_nr;ipsi_cd=ipsi_cd;contra_cd=contra_cd;
  end
  %save across areas
  ipsi_all(:,k)=ipsi_nr;
  contra_all(:,k)=contra_nr;
  ipsi_allcd(:,k)=ipsi_cd;
  contra_allcd(:,k)=contra_cd;
  end
  %remove respective injection area on the ipsi site
  idx_inj=[];idx_inj=find(nansum(ipsi_all)==max(nansum(ipsi_all)));
  disp([num2str(cortex_names{idx_inj,3}) ' was injected'])
    if  any(idx_inj==aa)==1
        idx_inj=idx_inj;
    else
        idx_inj=idx_inj-1;
    end
  ipsi_all(:,idx_inj)=[NaN;NaN;NaN;NaN;NaN;NaN];
  ipsi_allcd(:,idx_inj)=[NaN;NaN;NaN;NaN;NaN;NaN];

  ihemi=[];chemi=[];
  %calculate fraction across layers (minus L1)
        for l=2:6
        ihemi(l)=nansum(ipsi_all(l,:))/(nansum(nansum(ipsi_all))-nansum(nansum(ipsi_all(1,:))));
        chemi(l)=nansum(contra_all(l,:))/(nansum(nansum(contra_all))-nansum(nansum(contra_all(1,:))));
        end
  %save fractions across animlas regardless of area
  i_animal(:,i)=ihemi;
  c_animal(:,i)=chemi;
    area_li=[];area_lc=[];
    for m=1:length(cortex_names)
       temp_li=[];temp_lc=[];
       for l=2:6
       temp_li(l)=sum(ipsi_all(l,m))/(nansum(ipsi_all(:,m))-nansum(ipsi_all(1,m)));
       temp_lc(l)=sum(contra_all(l,m))/(nansum(contra_all(:,m))-nansum(contra_all(1,m)));
       end
    area_li(:,m)=temp_li;
    area_lc(:,m)=temp_lc;
   end
  i_areas_animal(:,:,i)=ipsi_all;
  c_areas_animal(:,:,i)=contra_all;
  i_areas_animalf(:,:,i)=area_li;
  c_areas_animalf(:,:,i)=area_lc;
  i_areas_animalcd(:,:,i)=ipsi_allcd;
  c_areas_animalcd(:,:,i)=contra_allcd;

end

end


%% 3 anatomy_indexcalc
function [indexcalc] = anatomy_indexcalc(layers)
for i=1:size(layers,2)
%1) iln=(L5+L6a/b) / (L2/3+L5+L6a/b);
iln(i)=(layers(4,i)+layers(5,i)+layers(6,i)) / (layers(2,i)+layers(4,i)+layers(5,i)+layers(6,i));
%2) L6i=(L6a/b-L2/3) / (L6a/b+L2/3);
L6i(i)=((layers(5,i)+layers(6,i))-layers(2,i)) / ((layers(5,i)+layers(6,i))+layers(2,i));
%3) L6ai=(L6a-L2/3) / (L6a+L2/3);
L6ai(i)=(layers(5,i)-layers(2,i)) / (layers(5,i)+layers(2,i));
%4) h_index=1-((L2/3+L4) / (L2/3+L4+L5));
h_index(i)=1-((layers(2,i)+layers(3,i)) / (layers(2,i)+layers(3,i)+layers(4,i)));
%5) pure L6a
%(L6a)/(L23 + L5 +L6ab)
pL6a(i)=(layers(5,i)) / (layers(2,i)+layers(4,i)+layers(5,i)+layers(6,i));
%6) fraction layer 6a/5
pL6aL5(i)=(layers(5,i)) / (layers(4,i));
%7)
%ILN without L6
iln_n6(i)=(layers(4,i)) / (layers(2,i)+layers(4,i));
%8)
%ILN without L5
iln_n5(i)=(layers(5,i)+layers(6,i)) / (layers(2,i)+layers(5,i)+layers(6,i));
% %9)
% %ILN without L5 and L6
% iln_n5(i)=(layers(5,i)+layers(6,i)) / (layers(2,i)+layers(5,i)+layers(6,i));

end

indexcalc=[iln' L6i' L6ai' h_index' pL6a' pL6aL5' iln_n6' iln_n5'];

% for k=1:size(indexcalc,2)
%     if sum(isnan(indexcalc(:,k)))>=size(indexcalc,2)-1
%     indexcalc(find(~isnan(indexcalc(:,k))),k)=NaN;
%     else
%     end
% end
end


%% 4 zero_area_layers
function [data3_m data4_m] = zero_area_layers(data3,data4,thr,an)
data=squeeze(nansum(data3(:,:,:)));
data_m=data;

data3_m=data3;
data4_m=data4;

for i=1:size(data,2)
        data_m(find(data(:,i)<thr),i)=0;
end

for i=1:size(data,2)
        data3_m(:,find(data(:,i)<thr),i)=zeros(6,length(find(data(:,i)<thr)));  
        data4_m(:,find(data(:,i)<thr),i)=zeros(6,length(find(data(:,i)<thr)));  

end



ind_a=find(sum(data_m>0,2)<an);
for i=1:size(data_m,2)
data3_m(:,ind_a,i)=zeros(6,length(ind_a));
data4_m(:,ind_a,i)=zeros(6,length(ind_a));
end

idxnan=find(sum(isnan(data3(3,:,:)),3)==size(data,2));
for i=1:size(data,2)
data3_m(3,idxnan,i)=NaN;
data4_m(3,idxnan,i)=NaN;
end

end
%% 5 anatomy_testarea

function [p_bin p_bin_c] = anatomy_testarea(hemi_idx,ipsi_val,contra_val)

for i=1:length(hemi_idx);
    h=[];p=[];pi=[];
   
   [h p] = ttest(hemi_idx(i,:),0);
    %[h p] = ttest(contra_val(i,:),ipsi_val(i,:),"Tail","right");
    try
    %[p h] = signrank(hemi_idx(i,:),0,"Tail","right");
     [jj(i) pp] = adtest(hemi_idx(i,:));
    catch
        %p=0.9;
         jj(i)=0;
    end
    p_bin_c(i)=p<0.05/(length(hemi_idx)-1);

    [h pi] = ttest(hemi_idx(i,:),0);
    %[h pi] = ttest(contra_val(i,:),ipsi_val(i,:),"Tail","left");
   try
     %[pi h] = signrank(hemi_idx(i,:),0,"Tail","left");
   catch
        %pi=0.9;
   end
    p_bin(i)=pi<0.05/(length(hemi_idx)-1);
end
end


%% 6 fILN_visualize

function fILN_visualize(x_cord,y_cord,fILNi,fILNc,color_scheme,title_pan)


%Colormap definition 
%tt=[];
%[tt]=buildcmap('xbwrz');
minbar=0;
maxbar=1;


% fig7= figure;set(fig7, 'Name', 'fILN plot');set(fig7, 'Position', [400, 300, 400, 400]);set(gcf,'color','w');t=tiledlayout('horizontal');
 fig7= figure;set(fig7, 'Name', 'fILN plot');set(fig7, 'Position', [400, 300, 450, 350]);set(gcf,'color','w');t=tiledlayout('horizontal');
t.TileSpacing = 'compact';t.Padding = 'compact';

   
%fILN=[]; fILN=iln_v{k};
nexttile
find(isnan(fILNi))
hold on;sc1=scatter(x_cord,y_cord,150,fILNi,'filled');sc1.MarkerEdgeColor='k';
hold on;sc1=scatter(x_cord(find(isnan(fILNi))),y_cord(find(isnan(fILNi))),150,[0.5 0.5 0.5],'filled');sc1.MarkerEdgeColor='k';

%hold on;sc1=scatter(x_cord,y_cord,150,'w','filled');sc1.MarkerEdgeColor='k';
set(gca,'Ydir','reverse');axis off;

%colormap(tt);
for i=1:45
 if i<10
text(x_cord(i)-18,y_cord(i)-5,num2str(i),"FontSize",9,'Color',[1 1 1]);%8
 else
text(x_cord(i)-45,y_cord(i)-5,num2str(i),"FontSize",9,'Color',[1 1 1]);%20
 end
end
if contains(title_pan,'VISp')==1
 hold on;sc1=scatter(x_cord(31),y_cord(31),150,'w','filled');sc1.MarkerEdgeColor=color_scheme;  sc1.LineWidth=1;sc1.Marker='o';
elseif contains(title_pan,'SSp-bfd')==1
    hold on;sc1=scatter(x_cord(18),y_cord(18),150,'w','filled');sc1.MarkerEdgeColor=color_scheme; sc1.Marker='o'; sc1.LineWidth=1;
else contains(title_pan,'MOp')==1
    hold on;sc1=scatter(x_cord(25),y_cord(25),150,'w','filled');sc1.MarkerEdgeColor=color_scheme; sc1.Marker='o'; sc1.LineWidth=1;
end
title('ipsi','FontWeight','normal',"FontSize",11);
%caxis([minbar maxbar]) 
%xlim([100 1500]) 
%ylim([100 1300]) 
nexttile
hold on;sc1=scatter(x_cord,y_cord,150,fILNc,'filled');sc1.MarkerEdgeColor='k';
hold on;sc1=scatter(x_cord(find(isnan(fILNc))),y_cord(find(isnan(fILNc))),150,[0.5 0.5 0.5],'filled');sc1.MarkerEdgeColor='k';

%hold on;sc1=scatter(x_cord,y_cord,150,'w','filled');sc1.MarkerEdgeColor='k';
set(gca,'Ydir','reverse');axis off;
set(gca,'Ydir','reverse');axis off;
set(gca,'xdir','reverse');axis off;
%colormap(tt);
colormap(slanCM('rainbow'))
for i=1:45
 if i<10
text(x_cord(i)+18,y_cord(i)-5,num2str(i),"FontSize",9,'Color',[1 1 1]);
 else
text(x_cord(i)+45,y_cord(i)-5,num2str(i),"FontSize",9,'Color',[1 1 1]);
 end
end
title('contra','FontWeight','normal',"FontSize",11);
text(1400,-100,title_pan,"Color",color_scheme,"FontSize",11)
%line([1500 1500], [0 1200],'LineStyle','--','Color','k')
cc=colorbar('southoutside');set(gca,'FontSize',9)
cc.Ticks = [0 0.25 0.5 0.75 1]
ylabel(cc,'fILN','FontSize',11,'Rotation',0);
caxis([minbar maxbar]) 
%xlim([100 1500])
%ylim([100 1300]) 

end
%% 7 iln_change

function [data_ilni data_ilnc rm1_error rm2_error order_i order_c temp_metric1 temp_metric2] = iln_change(rm1,rm2,idx)

for i=1:length(rm1)
temp_l=[];temp_l2=[];
temp_l=squeeze(rm1(:,i,:));
temp_l2=squeeze(rm2(:,i,:));
[i_metric(:,i,:)] = anatomy_indexcalc(temp_l);
[c_metric(:,i,:)] = anatomy_indexcalc(temp_l2);
end 

temp_metric1=i_metric(:,:,idx);temp_metric2=c_metric(:,:,idx);
temp_rm1=[];
temp_rm2=[];
temp_rm1=nanmean(temp_metric1);
temp_rm2=nanmean(temp_metric2);
rm1_error=nanstd(temp_metric1)/sqrt(size(temp_metric1,1));
rm2_error=nanstd(temp_metric2)/sqrt(size(temp_metric2,1));
[data_xxi order_i]=sort(temp_rm1,'ascend');
[data_xxc order_c]=sort(temp_rm2,'ascend');
data_ilnc=temp_rm2;
data_ilni=temp_rm1;
end
%% 8 anatomy_correlation

%calculates correlation values and average/SEM of data points to be correlted
function [value1_all value2_all value1 value2 xerr yerr r_a p_a r_avg p_avg value3_all xerr_c] = anatomy_correlation(rm1, rm2, idx)
%input to function: 
%rm1=
%rm2=
%idx=indexes to be correlated 
%output
%r_a/p_a= individual animal/modules based correlation/pvalues
%avg p_avg= correlation on average data points
%dependencies: anatomy_indexcalc.m

i_metric=[];c_metric=[];
%calculate indexes across areas
for i=1:size(rm1,2)
temp_l=[];temp_l2=[];
temp_l=squeeze(rm1(:,i,:));
temp_l2=squeeze(rm2(:,i,:));
[i_metric(:,i,:)] = anatomy_indexcalc(temp_l);
[c_metric(:,i,:)] = anatomy_indexcalc(temp_l2);
end 
temp_metric1=i_metric(:,:,idx);temp_metric2=c_metric(:,:,idx);

%Calculate correaltions across Ns
r_a=[];p_a=[];
for i=1:size(temp_metric1,1)
r=[];p=[];
[r p] = corr(temp_metric1(i,:)',(temp_metric2(i,:)-temp_metric1(i,:))','Type','Spearman','Rows','complete');
r_a(i)=r(1);
p_a(i)=p(1);
end
%average values 

value1=nanmean(temp_metric1);
value2=nanmean(temp_metric2)-nanmean(temp_metric1);
value3_all=temp_metric2;
value1_all=temp_metric1;
value2_all=temp_metric2-temp_metric1;
%xSEM
xerr=[];
xerr=nanstd(temp_metric1);
xerr_c=nanstd(temp_metric2);
temp_div=[];
for m=1:45
    temp_z=[];
    temp_z=temp_metric1(:,m);
temp_div(m)=sqrt(length(temp_z(~isnan(temp_z))));
end
xerr=xerr./temp_div;

for m=1:45
    temp_z=[];
    temp_z=temp_metric2(:,m);
temp_div(m)=sqrt(length(temp_z(~isnan(temp_z))));
end
xerr_c=xerr_c./temp_div;
%ySEM
yerr=[];
yerr=nanstd(temp_metric2-temp_metric1);
temp_div=[];
for m=1:45
    temp_z=[];
    temp_z=temp_metric2(:,m)-temp_metric1(:,m);
temp_div(m)=sqrt(length(temp_z(~isnan(temp_z))));
end
yerr=yerr./temp_div;

%average correlation 
[r_avg p_avg] = corr(value1',value2','Type','Spearman','Rows','complete');
end
%% 9 harmonic_mean_pvalue
function combined_p = harmonic_mean_pvalue(p_values)
    % Input:
    %   p_values - Array of p-values (e.g., [0.05, 0.01, 0.2])
    % Output:
    %   combined_p - Combined p-value using the Harmonic Mean p-value method

    % Ensure p-values are valid
    if any(p_values <= 0 | p_values > 1)
        error('P-values must be in the range (0, 1].');
    end

    % Number of p-values
    m = length(p_values);

    % Compute the harmonic mean p-value
    HMP = m / sum(1 ./ p_values);

    % Adjust the HMP to account for multiple testing
    combined_p = min(HMP * m, 1); % Ensure the combined p-value is capped at 1

    % Display results
    fprintf('Harmonic Mean Combined P-Value: %f\n', combined_p);
end
%% 10 anatomy_indexhemi
function [dchange] = anatomy_indexhemi(layers_i,layers_c)
for i=1:size(layers_i,2)
%1) iln=(L5+L6a/b) / (L2/3+L5+L6a/b);
dchange(i)=1/(abs(layers_c(:,i)-layers_i(:,i)))
end

end

%% 



%% add ons

%% M2 vs M1

fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 800, 300]);set(gcf,'color','w');
%p1=plot((nansum(i_m2aa))./nansum(nansum(i_m2aa)),'-ro');hold on;
%p2=plot(nanmean(((nansum(i_m1aa))./nansum(nansum(i_m1aa))),3),'-bo');hold on;
%p1=plot((nansum(i_m2aa))./nansum(nansum(i_m2aa))./nanmean(((nansum(i_m1aa))./nansum(nansum(i_m1aa))),3),'-ko');hold on;
% p1=plot(nanmean((nansum(i_m1aa)),3)./(nansum(i_m2aa)),'-bo');hold on;
% p2=plot(nanmean((nansum(i_m1aa)),3)./(nansum(i_m2aa)),'-ro');hold on;

p1=plot((nansum(i_m2aa))./nansum(nansum(i_m2aa)),'-ro');hold on;
p2=plot(nanmean((nansum(i_m1aa)),3)./nansum(nansum(i_m2aa)),'-bo');hold on;

box off;xticks([1:1:45]);xticklabels(cortex_names(:,2));
ylabel('Fraction norm. to M2')
set(gca,'TickDir','out');set(gca,'FontSize',11);legend({'M2 norm M2','M1 norm M2'})
title('ipsi');offsetAxes;ylim([-0.01 0.35])
%% 

fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 800, 300]);set(gcf,'color','w');
p1=plot((nansum(c_m2aa))./nansum(nansum(c_m2aa)),'-ro');hold on;
p2=plot(nanmean(((nansum(c_m1aa))./nansum(nansum(c_m1aa))),3),'-bo');hold on;
box off;xticks([1:1:45]);xticklabels(cortex_names(:,2));
ylabel('Fraction of cells')
set(gca,'TickDir','out');set(gca,'FontSize',11);legend({'M2','M1'})
title('contra');ylim([-0.01 0.5]);offsetAxes
%% 
fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 800, 300]);set(gcf,'color','w');
%p1=plot((nansum(i_m2aa))./nansum(nansum(i_m2aa)),'-ro');hold on;
%p2=plot(nanmean(((nansum(i_m1aa))./nansum(nansum(i_m1aa))),3),'-bo');hold on;
%p1=plot((nansum(i_m2aa))./nansum(nansum(i_m2aa))./nanmean(((nansum(i_m1aa))./nansum(nansum(i_m1aa))),3),'-ko');hold on;
% p1=plot(nanmean((nansum(i_m1aa)),3)./(nansum(i_m2aa)),'-bo');hold on;
% p2=plot(nanmean((nansum(i_m1aa)),3)./(nansum(i_m2aa)),'-ro');hold on;

p1=plot((nansum(c_m2aa))./nansum(nansum(c_m2aa)),'-ro');hold on;
p2=plot(nanmean((nansum(c_m1aa)),3)./nansum(nansum(c_m2aa)),'-bo');hold on;

box off;xticks([1:1:45]);xticklabels(cortex_names(:,2));
ylabel('Fraction norm. to M2')
set(gca,'TickDir','out');set(gca,'FontSize',11);legend({'M2 norm M2','M1 norm M2'})
title('contra');offsetAxes;ylim([-0.01 0.35])

%% 

fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 800, 300]);set(gcf,'color','w');
p1=plot((nansum(i_m2aa(:,27:36,:)))./nansum(nansum(i_m2aa(:,27:36,:))),'-ro');hold on;
p2=plot(nanmean(((nansum(i_m1aa(:,27:36,:)))./nansum(nansum(i_m1aa(:,27:36,:)))),3),'-bo');hold on;
box off;xticks([1:1:10]);xticklabels(cortex_names(27:36,2));
ylabel('Fraction of cells')
set(gca,'TickDir','out');set(gca,'FontSize',11);legend({'M2','M1'})
title('ipsi');ylim([-0.01 0.1]);offsetAxes