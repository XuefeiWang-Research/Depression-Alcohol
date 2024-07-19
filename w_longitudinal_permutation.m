%% permu_longitutional   MDD-audit(all)
clear;clc

% Merge all 10 threshold PRS data
root_path = '/Users/gray/Desktop/loadWXF/IMA_analy/GWAS/GWAStoPRS/result_PRSv125/mdd1VSalco1/seGWAS_BETAsum/';
i=1;  id = num2str(i*5, "%03d");
file = ['clumped',id];
file = [root_path,file];
file = [file, '/outpre_mdd/outpre_mdd.profile'];
selec_SNP = readtable(file,'FileType','text');
PRS=selec_SNP(:,[2,6]);   %select the IID and SCORE column;
PRS.Properties.VariableNames{2} = ['SCORE',id];

for i = 2:10
    id = num2str(i*5, "%03d");
    file = ['clumped',id];
    file = [root_path,file];
    file = [file, '/outpre_mdd/outpre_mdd.profile'];
    selec_SNP = readtable(file,'FileType','text');
    [~,inda,indb] = intersect(selec_SNP.IID,PRS.IID);
    selec_SNP=selec_SNP(inda,:);     PRS=PRS(indb,:);  
    PRS.P2=selec_SNP{:,6};%
    PRS.Properties.VariableNames{i+1} = ['SCORE',id];
end


% R_mean
load('/Users/gray/Desktop/loadWXF/IMA_analy/Apaper_result/MR/correctPermu/IMAGEN/MDD_audit.mat');

%  4 phenos
dirOutput = dir('/Users/gray/Desktop/loadWXF/IMA_analy/Behavior/independent_file/people_sitesex/drunk_audit/auditALL_*');
folder = string({dirOutput.folder}');
file = string({dirOutput.name}');
filepath = strcat(folder, '/', file);

IMA_BL=readtable(filepath{1}); IMA_BL.audit_UseProb = IMA_BL.audit_symp + IMA_BL.audit_prob;
IMA_FU1=readtable(filepath{2}); IMA_FU1.audit_UseProb = IMA_FU1.audit_symp + IMA_FU1.audit_prob;
IMA_FU2=readtable(filepath{3}); IMA_FU2.audit_UseProb = IMA_FU2.audit_symp + IMA_FU2.audit_prob;
IMA_FU3=readtable(filepath{4}); IMA_FU3.audit_UseProb = IMA_FU3.audit_symp + IMA_FU3.audit_prob;

IDs_union = union(union(union(IMA_BL.SubID, IMA_FU1.SubID), IMA_FU2.SubID), IMA_FU3.SubID);


    
for j=[12:14,27];
    phe_all_BL=IMA_BL;
    ind1=find(isnan(table2array(phe_all_BL(:,j)))==1);
    phe_all_BL(ind1,:)=[];
    [~,inda,indb] = intersect(PRS.IID,phe_all_BL.SubID);  
    PRS_all_BL=PRS(inda,:);  phe_all_BL=phe_all_BL(indb,:);   
    cova_all_BL = table2array(phe_all_BL(:,[18:26])); % warning

    phe_all_FU1=IMA_FU1;
    ind1=find(isnan(table2array(phe_all_FU1(:,j)))==1);
    phe_all_FU1(ind1,:)=[];
    [~,inda,indb] = intersect(PRS.IID,phe_all_FU1.SubID);  
    PRS_all_FU1=PRS(inda,:);     phe_all_FU1=phe_all_FU1(indb,:);   
    cova_all_FU1 = table2array(phe_all_FU1(:,[18:26])); % warning

    phe_all_FU2=IMA_FU2;
    ind1=find(isnan(table2array(phe_all_FU2(:,j)))==1);
    phe_all_FU2(ind1,:)=[];
    [~,inda,indb] = intersect(PRS.IID,phe_all_FU2.SubID);  
    PRS_all_FU2=PRS(inda,:);     phe_all_FU2=phe_all_FU2(indb,:);   
    cova_all_FU2 = table2array(phe_all_FU2(:,[18:26])); % warning

    phe_all_FU3=IMA_FU3;
    ind1=find(isnan(table2array(phe_all_FU3(:,j)))==1);
    phe_all_FU3(ind1,:)=[];
    [~,inda,indb] = intersect(PRS.IID,phe_all_FU3.SubID);  
    PRS_all_FU3=PRS(inda,:);     phe_all_FU3=phe_all_FU3(indb,:);   
    cova_all_FU3 = table2array(phe_all_FU3(:,[18:26])); % warning


    for k=1:1000
        rank_all = randperm(length(IDs_union));
        new_IDs=IDs_union(rank_all);

        [~,ia,ib] = intersect(new_IDs,phe_all_BL.SubID,'stable');
        phe_all_BLnew=phe_all_BL(ib,:);     
        beh_all_BLnew=table2array(phe_all_BLnew(:,j));
        [r1,p1]=partialcorr(PRS_all_BL{:,2:size(PRS_all_BL,2)},beh_all_BLnew,cova_all_BL,'Tail','right');
        re_per_all_BLr(:,k)=r1;  re_per_all_BLp(:,k)=p1;  re_per_all_BLnum(:,k) = length(PRS_all_BL{:,1});

        [~,ia,ib] = intersect(new_IDs,phe_all_FU1.SubID,'stable');
        phe_all_FU1new=phe_all_FU1(ib,:);     
        beh_all_FU1new=table2array(phe_all_FU1new(:,j));
        [r1,p1]=partialcorr(PRS_all_FU1{:,2:size(PRS_all_FU1,2)},beh_all_FU1new,cova_all_FU1,'Tail','right');
        re_per_all_FU1r(:,k)=r1;  re_per_all_FU1p(:,k)=p1;  re_per_all_FU1num(:,k) = length(PRS_all_FU1{:,1});

        [~,ia,ib] = intersect(new_IDs,phe_all_FU2.SubID,'stable');
        phe_all_FU2new=phe_all_FU2(ib,:);     
        beh_all_FU2new=table2array(phe_all_FU2new(:,j));
        [r1,p1]=partialcorr(PRS_all_FU2{:,2:size(PRS_all_FU2,2)},beh_all_FU2new,cova_all_FU2,'Tail','right');
        re_per_all_FU2r(:,k)=r1;  re_per_all_FU2p(:,k)=p1;  re_per_all_FU2num(:,k) = length(PRS_all_FU2{:,1});

        [~,ia,ib] = intersect(new_IDs,phe_all_FU3.SubID,'stable');
        phe_all_FU3new=phe_all_FU3(ib,:);     
        beh_all_FU3new=table2array(phe_all_FU3new(:,j));
        [r1,p1]=partialcorr(PRS_all_FU3{:,2:size(PRS_all_FU3,2)},beh_all_FU3new,cova_all_FU3,'Tail','right');
        re_per_all_FU3r(:,k)=r1;  re_per_all_FU3p(:,k)=p1;  re_per_all_FU3num(:,k) = length(PRS_all_FU3{:,1});

      end

    
    observe=[re_ob{j,1};re_ob{j,2};re_ob{j,3};re_ob{j,4}];   
    r_mean_all=mean(observe(:,1));  % all  

    % all  P
    mean_all2=mean([re_per_all_BLr;re_per_all_FU1r;re_per_all_FU2r;re_per_all_FU3r]);    
    m = length(find(mean_all2 > r_mean_all))   % if r_mean=positive,> ; else negative <
    p = m/1000

    result_all(j,1)=p;    
end

save('re_per_all_r.mat','re_per_all_BLr','re_per_all_FU1r','re_per_all_FU2r','re_per_all_FU3r');

re_per_all_BLr=array2table(re_per_all_BLr);
re_per_all_FU1r=array2table(re_per_all_FU1r);
re_per_all_FU2r=array2table(re_per_all_FU2r);
re_per_all_FU3r=array2table(re_per_all_FU3r);



writetable(re_per_all_BLr,'re_per_all_BLr.csv');
writetable(re_per_all_FU1r,'re_per_all_FU1r.csv');
writetable(re_per_all_FU2r,'re_per_all_FU2r.csv');
writetable(re_per_all_FU3r,'re_per_all_FU3r.csv');



%% permu_longitutional   MDD-audit(male)
clear;clc

% Merge all 10 threshold PRS data
root_path = '/Users/gray/Desktop/loadWXF/IMA_analy/GWAS/GWAStoPRS/result_PRSv125/mdd1VSalco1/seGWAS_BETAsum/';
i=1;  id = num2str(i*5, "%03d");
file = ['clumped',id];
file = [root_path,file];
file = [file, '/outpre_mdd/outpre_mdd.profile'];
selec_SNP = readtable(file,'FileType','text');
PRS=selec_SNP(:,[2,6]);   %select the IID and SCORE column;
PRS.Properties.VariableNames{2} = ['SCORE',id];

for i = 2:10
    id = num2str(i*5, "%03d");
    file = ['clumped',id];
    file = [root_path,file];
    file = [file, '/outpre_mdd/outpre_mdd.profile'];
    selec_SNP = readtable(file,'FileType','text');
    [~,inda,indb] = intersect(selec_SNP.IID,PRS.IID);
    selec_SNP=selec_SNP(inda,:);     PRS=PRS(indb,:);  
    PRS.P2=selec_SNP{:,6};%
    PRS.Properties.VariableNames{i+1} = ['SCORE',id];
end


% R_mean
load('/Users/gray/Desktop/loadWXF/IMA_analy/Apaper_result/MR/correctPermu/IMAGEN/MDD_audit.mat');

%  4 phenos
dirOutput = dir('/Users/gray/Desktop/loadWXF/IMA_analy/Behavior/independent_file/people_sitesex/drunk_audit/auditALL_*');
folder = string({dirOutput.folder}');
file = string({dirOutput.name}');
filepath = strcat(folder, '/', file);

IMA_BL=readtable(filepath{1}); IMA_BL.audit_UseProb = IMA_BL.audit_symp + IMA_BL.audit_prob;
IMA_FU1=readtable(filepath{2}); IMA_FU1.audit_UseProb = IMA_FU1.audit_symp + IMA_FU1.audit_prob;
IMA_FU2=readtable(filepath{3}); IMA_FU2.audit_UseProb = IMA_FU2.audit_symp + IMA_FU2.audit_prob;
IMA_FU3=readtable(filepath{4}); IMA_FU3.audit_UseProb = IMA_FU3.audit_symp + IMA_FU3.audit_prob;


%male
ind=find(IMA_BL.Gender_Male==1);      IMA_BL=IMA_BL(ind,:);
ind=find(IMA_FU1.Gender_Male==1);     IMA_FU1=IMA_FU1(ind,:);
ind=find(IMA_FU2.Gender_Male==1);     IMA_FU2=IMA_FU2(ind,:);
ind=find(IMA_FU3.Gender_Male==1);     IMA_FU3=IMA_FU3(ind,:);


IDs_union = union(union(union(IMA_BL.SubID, IMA_FU1.SubID), IMA_FU2.SubID), IMA_FU3.SubID);

    
for j=[12:14,27];
    phe_all_BL=IMA_BL;
    ind1=find(isnan(table2array(phe_all_BL(:,j)))==1);
    phe_all_BL(ind1,:)=[];
    [~,inda,indb] = intersect(PRS.IID,phe_all_BL.SubID);  
    PRS_all_BL=PRS(inda,:);  phe_all_BL=phe_all_BL(indb,:);   
    cova_all_BL = table2array(phe_all_BL(:,[18:24,26])); % warning

    phe_all_FU1=IMA_FU1;
    ind1=find(isnan(table2array(phe_all_FU1(:,j)))==1);
    phe_all_FU1(ind1,:)=[];
    [~,inda,indb] = intersect(PRS.IID,phe_all_FU1.SubID);  
    PRS_all_FU1=PRS(inda,:);     phe_all_FU1=phe_all_FU1(indb,:);   
    cova_all_FU1 = table2array(phe_all_FU1(:,[18:24,26])); % warning

    phe_all_FU2=IMA_FU2;
    ind1=find(isnan(table2array(phe_all_FU2(:,j)))==1);
    phe_all_FU2(ind1,:)=[];
    [~,inda,indb] = intersect(PRS.IID,phe_all_FU2.SubID);  
    PRS_all_FU2=PRS(inda,:);     phe_all_FU2=phe_all_FU2(indb,:);   
    cova_all_FU2 = table2array(phe_all_FU2(:,[18:24,26])); % warning

    phe_all_FU3=IMA_FU3;
    ind1=find(isnan(table2array(phe_all_FU3(:,j)))==1);
    phe_all_FU3(ind1,:)=[];
    [~,inda,indb] = intersect(PRS.IID,phe_all_FU3.SubID);  
    PRS_all_FU3=PRS(inda,:);     phe_all_FU3=phe_all_FU3(indb,:);   
    cova_all_FU3 = table2array(phe_all_FU3(:,[18:24,26])); % warning


    for k=1:1000
        rank_all = randperm(length(IDs_union));
        new_IDs=IDs_union(rank_all);

        [~,ia,ib] = intersect(new_IDs,phe_all_BL.SubID,'stable');
        phe_all_BLnew=phe_all_BL(ib,:);     
        beh_all_BLnew=table2array(phe_all_BLnew(:,j));
        [r1,p1]=partialcorr(PRS_all_BL{:,2:size(PRS_all_BL,2)},beh_all_BLnew,cova_all_BL,'Tail','right');
        re_per_all_BLr(:,k)=r1;  re_per_all_BLp(:,k)=p1;  re_per_all_BLnum(:,k) = length(PRS_all_BL{:,1});

        [~,ia,ib] = intersect(new_IDs,phe_all_FU1.SubID,'stable');
        phe_all_FU1new=phe_all_FU1(ib,:);     
        beh_all_FU1new=table2array(phe_all_FU1new(:,j));
        [r1,p1]=partialcorr(PRS_all_FU1{:,2:size(PRS_all_FU1,2)},beh_all_FU1new,cova_all_FU1,'Tail','right');
        re_per_all_FU1r(:,k)=r1;  re_per_all_FU1p(:,k)=p1;  re_per_all_FU1num(:,k) = length(PRS_all_FU1{:,1});

        [~,ia,ib] = intersect(new_IDs,phe_all_FU2.SubID,'stable');
        phe_all_FU2new=phe_all_FU2(ib,:);     
        beh_all_FU2new=table2array(phe_all_FU2new(:,j));
        [r1,p1]=partialcorr(PRS_all_FU2{:,2:size(PRS_all_FU2,2)},beh_all_FU2new,cova_all_FU2,'Tail','right');
        re_per_all_FU2r(:,k)=r1;  re_per_all_FU2p(:,k)=p1;  re_per_all_FU2num(:,k) = length(PRS_all_FU2{:,1});

        [~,ia,ib] = intersect(new_IDs,phe_all_FU3.SubID,'stable');
        phe_all_FU3new=phe_all_FU3(ib,:);     
        beh_all_FU3new=table2array(phe_all_FU3new(:,j));
        [r1,p1]=partialcorr(PRS_all_FU3{:,2:size(PRS_all_FU3,2)},beh_all_FU3new,cova_all_FU3,'Tail','right');
        re_per_all_FU3r(:,k)=r1;  re_per_all_FU3p(:,k)=p1;  re_per_all_FU3num(:,k) = length(PRS_all_FU3{:,1});

    end


    observe=[re_ob{j,1};re_ob{j,2};re_ob{j,3};re_ob{j,4}];   
    r_mean_male=mean(observe(:,3));  % male  

    % all  P
    mean_all2=mean([re_per_all_BLr;re_per_all_FU1r;re_per_all_FU2r;re_per_all_FU3r]);    
    m = length(find(mean_all2 > r_mean_male))   % if r_mean=positive,> ; else negative <
    p = m/1000

    result_male(j,1)=p;  

end

%% permu_longitutional   MDD-audit(female)
clear;clc

% Merge all 10 threshold PRS data
root_path = '/Users/gray/Desktop/loadWXF/IMA_analy/GWAS/GWAStoPRS/result_PRSv125/mdd1VSalco1/seGWAS_BETAsum/';
i=1;  id = num2str(i*5, "%03d");
file = ['clumped',id];
file = [root_path,file];
file = [file, '/outpre_mdd/outpre_mdd.profile'];
selec_SNP = readtable(file,'FileType','text');
PRS=selec_SNP(:,[2,6]);   %select the IID and SCORE column;
PRS.Properties.VariableNames{2} = ['SCORE',id];

for i = 2:10
    id = num2str(i*5, "%03d");
    file = ['clumped',id];
    file = [root_path,file];
    file = [file, '/outpre_mdd/outpre_mdd.profile'];
    selec_SNP = readtable(file,'FileType','text');
    [~,inda,indb] = intersect(selec_SNP.IID,PRS.IID);
    selec_SNP=selec_SNP(inda,:);     PRS=PRS(indb,:);  
    PRS.P2=selec_SNP{:,6};%
    PRS.Properties.VariableNames{i+1} = ['SCORE',id];
end


% R_mean
load('/Users/gray/Desktop/loadWXF/IMA_analy/Apaper_result/MR/correctPermu/IMAGEN/MDD_audit.mat');

%  4 phenos
dirOutput = dir('/Users/gray/Desktop/loadWXF/IMA_analy/Behavior/independent_file/people_sitesex/drunk_audit/auditALL_*');
folder = string({dirOutput.folder}');
file = string({dirOutput.name}');
filepath = strcat(folder, '/', file);

IMA_BL=readtable(filepath{1}); IMA_BL.audit_UseProb = IMA_BL.audit_symp + IMA_BL.audit_prob;
IMA_FU1=readtable(filepath{2}); IMA_FU1.audit_UseProb = IMA_FU1.audit_symp + IMA_FU1.audit_prob;
IMA_FU2=readtable(filepath{3}); IMA_FU2.audit_UseProb = IMA_FU2.audit_symp + IMA_FU2.audit_prob;
IMA_FU3=readtable(filepath{4}); IMA_FU3.audit_UseProb = IMA_FU3.audit_symp + IMA_FU3.audit_prob;


%female
ind=find(IMA_BL.Gender_Male==0);      IMA_BL=IMA_BL(ind,:);
ind=find(IMA_FU1.Gender_Male==0);     IMA_FU1=IMA_FU1(ind,:);
ind=find(IMA_FU2.Gender_Male==0);     IMA_FU2=IMA_FU2(ind,:);
ind=find(IMA_FU3.Gender_Male==0);     IMA_FU3=IMA_FU3(ind,:);


IDs_union = union(union(union(IMA_BL.SubID, IMA_FU1.SubID), IMA_FU2.SubID), IMA_FU3.SubID);


    
for j=[12:14,27];
    phe_all_BL=IMA_BL;
    ind1=find(isnan(table2array(phe_all_BL(:,j)))==1);
    phe_all_BL(ind1,:)=[];
    [~,inda,indb] = intersect(PRS.IID,phe_all_BL.SubID);  
    PRS_all_BL=PRS(inda,:);  phe_all_BL=phe_all_BL(indb,:);   
    cova_all_BL = table2array(phe_all_BL(:,[18:24,26])); % warning

    phe_all_FU1=IMA_FU1;
    ind1=find(isnan(table2array(phe_all_FU1(:,j)))==1);
    phe_all_FU1(ind1,:)=[];
    [~,inda,indb] = intersect(PRS.IID,phe_all_FU1.SubID);  
    PRS_all_FU1=PRS(inda,:);     phe_all_FU1=phe_all_FU1(indb,:);   
    cova_all_FU1 = table2array(phe_all_FU1(:,[18:24,26])); % warning

    phe_all_FU2=IMA_FU2;
    ind1=find(isnan(table2array(phe_all_FU2(:,j)))==1);
    phe_all_FU2(ind1,:)=[];
    [~,inda,indb] = intersect(PRS.IID,phe_all_FU2.SubID);  
    PRS_all_FU2=PRS(inda,:);     phe_all_FU2=phe_all_FU2(indb,:);   
    cova_all_FU2 = table2array(phe_all_FU2(:,[18:24,26])); % warning

    phe_all_FU3=IMA_FU3;
    ind1=find(isnan(table2array(phe_all_FU3(:,j)))==1);
    phe_all_FU3(ind1,:)=[];
    [~,inda,indb] = intersect(PRS.IID,phe_all_FU3.SubID);  
    PRS_all_FU3=PRS(inda,:);     phe_all_FU3=phe_all_FU3(indb,:);   
    cova_all_FU3 = table2array(phe_all_FU3(:,[18:24,26])); % warning


    for k=1:1000
        rank_all = randperm(length(IDs_union));
        new_IDs=IDs_union(rank_all);

        [~,ia,ib] = intersect(new_IDs,phe_all_BL.SubID,'stable');
        phe_all_BLnew=phe_all_BL(ib,:);     
        beh_all_BLnew=table2array(phe_all_BLnew(:,j));
        [r1,p1]=partialcorr(PRS_all_BL{:,2:size(PRS_all_BL,2)},beh_all_BLnew,cova_all_BL,'Tail','right');
        re_per_all_BLr(:,k)=r1;  re_per_all_BLp(:,k)=p1;  re_per_all_BLnum(:,k) = length(PRS_all_BL{:,1});

        [~,ia,ib] = intersect(new_IDs,phe_all_FU1.SubID,'stable');
        phe_all_FU1new=phe_all_FU1(ib,:);     
        beh_all_FU1new=table2array(phe_all_FU1new(:,j));
        [r1,p1]=partialcorr(PRS_all_FU1{:,2:size(PRS_all_FU1,2)},beh_all_FU1new,cova_all_FU1,'Tail','right');
        re_per_all_FU1r(:,k)=r1;  re_per_all_FU1p(:,k)=p1;  re_per_all_FU1num(:,k) = length(PRS_all_FU1{:,1});

        [~,ia,ib] = intersect(new_IDs,phe_all_FU2.SubID,'stable');
        phe_all_FU2new=phe_all_FU2(ib,:);     
        beh_all_FU2new=table2array(phe_all_FU2new(:,j));
        [r1,p1]=partialcorr(PRS_all_FU2{:,2:size(PRS_all_FU2,2)},beh_all_FU2new,cova_all_FU2,'Tail','right');
        re_per_all_FU2r(:,k)=r1;  re_per_all_FU2p(:,k)=p1;  re_per_all_FU2num(:,k) = length(PRS_all_FU2{:,1});

        [~,ia,ib] = intersect(new_IDs,phe_all_FU3.SubID,'stable');
        phe_all_FU3new=phe_all_FU3(ib,:);     
        beh_all_FU3new=table2array(phe_all_FU3new(:,j));
        [r1,p1]=partialcorr(PRS_all_FU3{:,2:size(PRS_all_FU3,2)},beh_all_FU3new,cova_all_FU3,'Tail','right');
        re_per_all_FU3r(:,k)=r1;  re_per_all_FU3p(:,k)=p1;  re_per_all_FU3num(:,k) = length(PRS_all_FU3{:,1});

    end

    observe=[re_ob{j,1};re_ob{j,2};re_ob{j,3};re_ob{j,4}];   
    r_mean_female=mean(observe(:,5));  % female 

    % all  P
    mean_all2=mean([re_per_all_BLr;re_per_all_FU1r;re_per_all_FU2r;re_per_all_FU3r]);    
    m = length(find(mean_all2 > r_mean_female))   % if r_mean=positive,> ; else negative <
    p = m/1000

    result_female(j,1)=p;  

end

