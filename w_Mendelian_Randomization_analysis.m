%% Prepare GWAS and selectPRS %%
clear;clc
MDD1=readtable('/home1/WangXF/IMA_analy/GWAS/GWAStoPRS/BASE/MDD1.txt');
ALcohol1=readtable('/home1/WangXF/IMA_analy/GWAS/GWAStoPRS/BASE/Alcohol1.txt');
[~,id1,id2]=intersect(ALcohol1(:,2),MDD1(:,2));   %SNP
al1=ALcohol1(id1,:);
md1=MDD1(id2,:);

beta_mdd=log(table2array(md1(:,9)));           %beta =log(OR)
md1(:,9)=array2table(beta_mdd);
md1.Properties.VariableNames(9)={'BETA'};


root_path='/home1/WangXF/IMA_analy/GWAS/GWAStoPRS/result_PRSv125/mdd1VSalco1/selectGWAS/';
for i = 1:10
    id = num2str(i*5, "%03d");
    file = ['clumped',id];
    file_folder = [root_path,file];
    unix(['mkdir -p ', file_folder]);
    thr = 0.05 * i;
    
    %%%%%%%%%  remove pleiotropic SNPs %%%%%%%% %  0.05:0.05:0.5  % %%%%
    ind1=find((al1.P<0.05)&(md1.P>thr));  %Alco to mdd causal
    select_alco1=al1(ind1,:);
    output_alco = [file_folder, '/select_alco1'];
    writetable(select_alco1, output_alco, 'Delimiter',' ','WriteVariableNames',false);
    
    output_mdd1 = [file_folder, '/select_mdd1'];
    ind2=find((md1.P<0.05)&(al1.P>thr));  %mdd to Alco causal
    select_mdd1=md1(ind2,:);
    writetable(select_mdd1, output_mdd1, 'Delimiter',' ','WriteVariableNames',false);
    
    
    %%%%%%%%%%%%%%%  calculate PRS %%%%%%%%%%%%%%%%
    %mdd
    unix(['mkdir -p ', file_folder, '/outpre_mdd/']);
    snp='/home1/WangXF/IMA_analy/GWAS/GWAStoPRS/result_PRSv125/newPRS/mdd1/flipped_target';  %after CLUMPING .profile
    outpre=[file_folder, '/outpre_mdd/outpre_mdd'];  % output name
    unix(['plink1.9 --bfile ',snp,' --score ', file_folder, '/select_mdd1.txt 2 4 9 sum --out ',outpre]); % SNP,A1,BETA
    
    %alco
    unix(['mkdir -p ', file_folder, '/outpre_alco/']);
    snp='/home1/WangXF/IMA_analy/GWAS/GWAStoPRS/result_PRSv125/newPRS/Alcohol1/flipped_target';   %after CLUMPING .profile
    outpre=[file_folder, 'outpre_alco/outpre_alco'];  % output name
    unix(['plink1.9 --bfile ',snp,' --score ', file_folder, '/select_alco1.txt 2 4 7 --out ',outpre]);  % SNP,A1,BETA

end

%% MR analysis and single timepoint permutation  %%

% 10 threshold ALCO-audit_freq 

clear;clc
IMA=readtable('/Users/gray/Desktop/loadWXF/IMA_analy/Behavior/independent_file/people_sitesex/drunk_audit/auditALL_BL.csv');

%sex
a=table2array(IMA(:,25));    
ind1=find(a==1);    IMA_M = IMA(ind1,:);   % 1=male; 0=female
ind1=find(a==0);    IMA_F = IMA(ind1,:);  

root_path = '/Users/gray/Desktop/loadWXF/IMA_analy/GWAS/GWAStoPRS/result_PRSv125/mdd1VSalco1/selectGWASforPRS/';
result_observe=zeros(10,9);
result_permu_all=zeros(10,1000);
result_permu_male=zeros(10,1000);
result_permu_female=zeros(10,1000);

for i = 1:10
    i
    id = num2str(i*5, "%03d");
    file = ['clumped',id];
    file = [root_path,file];
    file = [file, '/outpre_alco/outpre_alco.profile'];
    selec_SNP = readtable(file,'FileType','text');

    for j=12 % behavior
        hcp_all=IMA;
        ind1=find(isnan(table2array(hcp_all(:,j)))==1);
        hcp_all(ind1,:)=[];
        [~,inda,indb] = intersect(selec_SNP.IID,hcp_all.SubID);  
        selec_SNP_all=selec_SNP(inda,:);     hcp_all=hcp_all(indb,:);  
        PRS_all=table2array(selec_SNP_all(:,6));   
        cova_all = table2array(hcp_all(:,[18:25])); % warning
        beh_all=table2array(hcp_all(:,j));
        [r1,p1]=partialcorr(PRS_all,beh_all,cova_all,'Tail','right');
        result_observe(i,1) = r1;
        result_observe(i,2) = p1;  % one-tailed
        result_observe(i,7) = length(PRS_all);

        hcp_male=IMA_M;
        ind1=find(isnan(table2array(hcp_male(:,j)))==1);
        hcp_male(ind1,:)=[];
        [~,inda,indb] = intersect(selec_SNP.IID,hcp_male.SubID);  
        selec_SNP_male=selec_SNP(inda,:);     hcp_male=hcp_male(indb,:);  
        PRS_male=table2array(selec_SNP_male(:,6));   
        cova_male = table2array(hcp_male(:,[18:24])); % warning
        beh_male=table2array(hcp_male(:,j));
        [r2,p2]=partialcorr(PRS_male,beh_male,cova_male,'Tail','right');
        result_observe(i,3) = r2;
        result_observe(i,4) = p2;  % one-tailed
        result_observe(i,8) = length(PRS_male);


        hcp_female=IMA_F;
        ind1=find(isnan(table2array(hcp_female(:,j)))==1);
        hcp_female(ind1,:)=[];
        [~,inda,indb] = intersect(selec_SNP.IID,hcp_female.SubID);  
        selec_SNP_female=selec_SNP(inda,:);     hcp_female=hcp_female(indb,:);  
        PRS_female=table2array(selec_SNP_female(:,6));   
        cova_female = table2array(hcp_female(:,[18:24])); % warning
        beh_female=table2array(hcp_female(:,j));
        [r3,p3]=partialcorr(PRS_female,beh_female,cova_female,'Tail','right');
        result_observe(i,5) = r3;
        result_observe(i,6) = p3;  % one-tailed
        result_observe(i,9) = length(PRS_female);



        for k=1:1000
            rank_all = randperm(length(beh_all));
            beh_rank_all=beh_all(rank_all);
            cova_rank_all=cova_all(rank_all,:);
            [r4,p4]=partialcorr(PRS_all,beh_rank_all,cova_rank_all);
            result_permu_all(i,k)=r4;
        end

        for l=1:1000
            rank_male = randperm(length(beh_male));
            beh_rank_male=beh_male(rank_male);
            cova_rank_male=cova_male(rank_male,:);
            [r5,p5]=partialcorr(PRS_male,beh_rank_male,cova_rank_male);
            result_permu_male(i,l)=r5;
        end

        for m=1:1000
            rank_female = randperm(length(beh_female));
            beh_rank_female=beh_female(rank_female);
            cova_rank_female=cova_female(rank_female,:);
            [r6,p6]=partialcorr(PRS_female,beh_rank_female,cova_rank_female);
            result_permu_female(i,m)=r6;
        end

    end

end


result_P=zeros(1,3);
% all
mean_all=mean(result_permu_all);  r_mean=mean(result_observe(:,1))  
m = length(find(mean_all > r_mean))   % if r_mean=positive,> ; else negative <
p_all = m/1000;  result_P(1,1)=p_all;
% male
mean_male=mean(result_permu_male);
r_mean=mean(result_observe(:,3))  
m = length(find(mean_male > r_mean))   % if r_mean=positive,> ; else negative <
p_male = m/1000;    result_P(1,2)=p_male;
% female
mean_female=mean(result_permu_female);
r_mean=mean(result_observe(:,5))  
m = length(find(mean_female > r_mean))   % if r_mean=positive,> ; else negative <
p_female = m/1000;    result_P(1,3)=p_female;



%%  MDD_PRS with audit Symp+Prob  %%
clear;clc
mdd=readtable('/Users/gray/Desktop/loadWXF/IMA_analy/GWAS/GWAStoPRS/result_PRSv125/newPRS/mdd1/PRSice_SCORES_AT_ALL_THRESHOLDS.txt');
depbeh=readtable('/Users/gray/Desktop/loadWXF/IMA_analy/Behavior/independent_file/people_sitesex/drunk_audit/auditALL_FU3.csv');
depbeh.audit_sum=depbeh.audit_symp + depbeh.audit_prob;
% all
[~,inda,indb] = intersect(depbeh.SubID,mdd.IID);  
depbeh_all=depbeh(inda,:);   mdd_all=mdd(indb,:);
best_all=table2array(mdd_all(:,11));  beh_all=table2array(depbeh_all(:,27));  cova_all=table2array(depbeh_all(:,[18:26]));
[r,p]=partialcorr(best_all,beh_all,cova_all,'Tail','right') %one-tailed
result=zeros(1,9);   result(1,1)=r;   result(1,2)=p;   result(1,3)=length(best_all);

a=table2array(depbeh(:,25));    
ind1=find(a==1);    depbeh_M=depbeh(ind1,:);   %male
ind2=find(a==0);    depbeh_F=depbeh(ind2,:);   %female
% male
[~,inda,indb] = intersect(depbeh_M.SubID,mdd.IID);  
depbeh_M=depbeh_M(inda,:);   mdd_M=mdd(indb,:);
best_M=table2array(mdd_M(:,11));  beh_M=table2array(depbeh_M(:,27)); cova_M=table2array(depbeh_M(:,[18:24,26]));
[r,p]=partialcorr(best_M,beh_M,cova_M,'Tail','right')
result(1,4)=r;   result(1,5)=p;   result(1,6)=length(best_M);
% female
[~,inda,indb] = intersect(depbeh_F.SubID,mdd.IID);  
depbeh_F=depbeh_F(inda,:);   mdd_F=mdd(indb,:);
best_F=table2array(mdd_F(:,11));  beh_F=table2array(depbeh_F(:,27));  cova_F=table2array(depbeh_F(:,[18:24,26]));
[r,p]=partialcorr(best_F,beh_F,cova_F,'Tail','right')
result(1,7)=r;   result(1,8)=p;   result(1,9)=length(best_F);

%% 10 threshold MDD-audit Symp+Prob  %%
clear;clc
IMA=readtable('/Users/gray/Desktop/loadWXF/IMA_analy/Behavior/independent_file/people_sitesex/drunk_audit/auditALL_FU3.csv');
IMA.audit_sum=IMA.audit_symp + IMA.audit_prob;

%sex
a=table2array(IMA(:,25));    
ind1=find(a==1);    IMA_M = IMA(ind1,:);   % 1=male; 0=female
ind1=find(a==0);    IMA_F = IMA(ind1,:);  

root_path = '/Users/gray/Desktop/loadWXF/IMA_analy/GWAS/GWAStoPRS/result_PRSv125/mdd1VSalco1/seGWAS_BETAsum/';
result_observe=zeros(10,9);
result_permu_all=zeros(10,1000);
result_permu_male=zeros(10,1000);
result_permu_female=zeros(10,1000);

for i = 1:10
    i
    id = num2str(i*5, "%03d");
    file = ['clumped',id];
    file = [root_path,file];
    file = [file, '/outpre_mdd/outpre_mdd.profile'];
    selec_SNP = readtable(file,'FileType','text');

    for j=27 % behavior
        hcp_all=IMA;
        ind1=find(isnan(table2array(hcp_all(:,j)))==1);
        hcp_all(ind1,:)=[];
        [~,inda,indb] = intersect(selec_SNP.IID,hcp_all.SubID);  
        selec_SNP_all=selec_SNP(inda,:);     hcp_all=hcp_all(indb,:);  
        PRS_all=table2array(selec_SNP_all(:,6));   
        cova_all = table2array(hcp_all(:,[18:26])); % warning
        beh_all=table2array(hcp_all(:,j));
        [r1,p1]=partialcorr(PRS_all,beh_all,cova_all,'Tail','right');
        result_observe(i,1) = r1;
        result_observe(i,2) = p1;  % one-tailed
        result_observe(i,7) = length(PRS_all);

        hcp_male=IMA_M;
        ind1=find(isnan(table2array(hcp_male(:,j)))==1);
        hcp_male(ind1,:)=[];
        [~,inda,indb] = intersect(selec_SNP.IID,hcp_male.SubID);  
        selec_SNP_male=selec_SNP(inda,:);     hcp_male=hcp_male(indb,:);  
        PRS_male=table2array(selec_SNP_male(:,6));   
        cova_male = table2array(hcp_male(:,[18:24,26])); % warning
        beh_male=table2array(hcp_male(:,j));
        [r2,p2]=partialcorr(PRS_male,beh_male,cova_male,'Tail','right');
        result_observe(i,3) = r2;
        result_observe(i,4) = p2;  % one-tailed
        result_observe(i,8) = length(PRS_male);


        hcp_female=IMA_F;
        ind1=find(isnan(table2array(hcp_female(:,j)))==1);
        hcp_female(ind1,:)=[];
        [~,inda,indb] = intersect(selec_SNP.IID,hcp_female.SubID);  
        selec_SNP_female=selec_SNP(inda,:);     hcp_female=hcp_female(indb,:);  
        PRS_female=table2array(selec_SNP_female(:,6));   
        cova_female = table2array(hcp_female(:,[18:24,26])); % warning
        beh_female=table2array(hcp_female(:,j));
        [r3,p3]=partialcorr(PRS_female,beh_female,cova_female,'Tail','right');
        result_observe(i,5) = r3;
        result_observe(i,6) = p3;  % one-tailed
        result_observe(i,9) = length(PRS_female);



        for k=1:1000
            rank_all = randperm(length(beh_all));
            beh_rank_all=beh_all(rank_all);
            cova_rank_all=cova_all(rank_all,:);
            [r4,p4]=partialcorr(PRS_all,beh_rank_all,cova_rank_all);
            result_permu_all(i,k)=r4;
        end

        for l=1:1000
            rank_male = randperm(length(beh_male));
            beh_rank_male=beh_male(rank_male);
            cova_rank_male=cova_male(rank_male,:);
            [r5,p5]=partialcorr(PRS_male,beh_rank_male,cova_rank_male);
            result_permu_male(i,l)=r5;
        end

        for m=1:1000
            rank_female = randperm(length(beh_female));
            beh_rank_female=beh_female(rank_female);
            cova_rank_female=cova_female(rank_female,:);
            [r6,p6]=partialcorr(PRS_female,beh_rank_female,cova_rank_female);
            result_permu_female(i,m)=r6;
        end

    end

end


result_P=zeros(1,3);
% all
mean_all=mean(result_permu_all);  r_mean=mean(result_observe(:,1))  
m = length(find(mean_all > r_mean))   % if r_mean=positive,> ; else negative <
p_all = m/1000;  result_P(1,1)=p_all;
% male
mean_male=mean(result_permu_male);
r_mean=mean(result_observe(:,3))  
m = length(find(mean_male > r_mean))   % if r_mean=positive,> ; else negative <
p_male = m/1000;    result_P(1,2)=p_male;
% female
mean_female=mean(result_permu_female);
r_mean=mean(result_observe(:,5))  
m = length(find(mean_female > r_mean))   % if r_mean=positive,> ; else negative <
p_female = m/1000;    result_P(1,3)=p_female;
