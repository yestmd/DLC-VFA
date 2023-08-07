#!/usr/bin/env python
# coding: utf-8

# ## To combine the VFA orientation from 1st and 2nd cohort data

# In[1]:


import glob as glob
import pandas as pd
import os
import shutil
import numpy as np
from datetime import datetime
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.stats import sem
from statsmodels.stats.anova import AnovaRM
import scikit_posthocs as sp
import statsmodels.formula.api as smf


# In[2]:


import statsmodels.formula.api as smf
import scipy.stats as stats
import statsmodels.api as sm
import researchpy as rp


# In[3]:


import pingouin as pg


# In[4]:


#sync_temp2 path
coh1_male_path="D:\\CPT_TOT_DS_Stage_3\\CPT_DS_TOT_Visual_analysis\\vfa male tot_based on 1000000 model_fixed\\dlc_results_namechanged\\results\\temp_2"
coh1_female_path="D:\\CPT_TOT_DS_Stage_3\\CPT_DS_TOT_Visual_analysis\\vfa female tot_based on 1030000 model\\results\\temp_2"
coh2_path="D:\\CPT_TOT_DS_Stage_3_2nd\\VFA_baseline\\DLC_output\\results\\temp_2"


# ### plot 45min bin blind vs oriented time

# In[5]:


coh1_male_code=pd.read_csv(r"D:\CPT_TOT_DS_Stage_3\CPT_DS_TOT_Visual_analysis\file_frame_info_male.csv", index_col=False)
coh1_female_code=pd.read_csv(r"D:\CPT_TOT_DS_Stage_3\CPT_DS_TOT_Visual_analysis\file_frame_info_female.csv", index_col=False)
coh2_tot_code=pd.read_csv(r"D:\CPT_TOT_DS_Stage_3_2nd\VFA_baseline\coh2_tot_code.csv", index_col=False)


# In[6]:


coh1_male_code.head(1)


# In[7]:


coh1_female_code.head(1)


# In[8]:


coh2_tot_code.head(1)


# In[9]:


#generate 45min bin group df for each session
path_list=[coh1_male_path, coh1_female_path, coh2_path]
coh=["coh1", "coh1", "coh2"]
delte_str=[4,12,12]
code_file_name=[coh1_male_code["file_name"].tolist(), coh1_female_code["file_name"].tolist(), coh2_tot_code["vfa_csv_2"].tolist()]
nfps=[coh1_male_code["n_fps"].tolist(), coh1_female_code["n_fps"].tolist(), coh2_tot_code["n_fps"].tolist()]

for i,m, h,j, d in zip(path_list,code_file_name, nfps, coh, delte_str) :
    os.chdir(i)
    
    save_path=i[:-7]+"\\45_oriented"
    if os.path.isdir(save_path) == False:
        os.mkdir(save_path)
        
    for n, b in zip(m, h):
        df=pd.read_csv(n, index_col=[0])
        df.dropna(subset=["index", "total_visual"], inplace=True)
        df=df.groupby(by=["45min_bin", "blind_or_not"], as_index=False).count()
        df["cohort"]=j
        df["dt_id"]=n[:-d]
        df["animal_id"]=df["dt_id"][1][11:]
        df["id_cohort"]=df["animal_id"]+ "-" + df["cohort"]                    
        df["duration"]=df['index']/b
        df.to_csv(save_path + "\\%s"%(n), index=False)


# In[10]:


#generate file list from three parts of data
path_list=[coh1_male_path, coh1_female_path, coh2_path]
df_list=[]

for i in path_list:
    file_list=glob.glob(i[:-7]+"\\45_oriented\\*.csv")
    df_list=df_list+file_list
    
#combine 45min oriented bin df from three parts

df=pd.DataFrame()
for i in df_list:
    df_=pd.read_csv(i)
    if len(df_)<8:
        print(i)
    df=pd.concat([df, df_], ignore_index=True)
df['group_id']=df['animal_id'].apply(lambda x:x[0])
df['gander']=np.where(df['group_id'].isin(["W", "X", "Y", "Z"]), "male", "female")

#plot tot baseline and state
title="male_tot_baseline_45min_bin"
t=df.groupby(by=["id_cohort", "45min_bin", "blind_or_not", "gander"], as_index=False)['duration'].mean()
t=t.rename(columns={"45min_bin": "min_45bin"})
plot=sns.pointplot(data=t[(t["min_45bin"]!="after_session") & (t["min_45bin"]!="before_session") & (t["gander"]=="male")], 
                   x="min_45bin", y="duration", hue="blind_or_not", capsize=.01, errorbar=('ci', 68))
plot.set(xlabel="")
plot.set(ylabel="total_duration")
plt.title(title)
plt.xticks(rotation = 45, ha='right')

plt.legend(bbox_to_anchor=(1.25, 1),borderaxespad=0)
   
save_path=r"D:\CPT_TOT_manuscript\data processing\fig\vfa_45"
if os.path.isdir(save_path) == False:
    os.mkdir(save_path)
    
plt.savefig(save_path+"\%s.png"%(title), dpi=800, bbox_inches="tight")

try:
    print(title)
    data=t[(t["min_45bin"]!="after_session") & (t["min_45bin"]!="before_session") & (t["gander"]=="male")]
    aovrm = AnovaRM(data, 'duration', 'id_cohort', within=['min_45bin', 'blind_or_not'])
    res = aovrm.fit()

    print(res)
except ValueError as e:
    print(str(e))


# In[11]:


title="female_tot_baseline_45min_bin"
t=df.groupby(by=["id_cohort", "45min_bin", "blind_or_not", "gander"], as_index=False)['duration'].mean()
t=t.rename(columns={"45min_bin": "min_45bin"})
plot=sns.pointplot(data=t[(t["min_45bin"]!="after_session") & (t["min_45bin"]!="before_session") & (t["gander"]=="female")], 
                   x="min_45bin", y="duration", hue="blind_or_not", capsize=.01, errorbar=('ci', 68))
plot.set(xlabel="")
plot.set(ylabel="total_duration")
plt.title(title)
plt.xticks(rotation = 45, ha='right')

plt.legend(bbox_to_anchor=(1.25, 1),borderaxespad=0)
   
save_path=r"D:\CPT_TOT_manuscript\data processing\fig\vfa_45\%s.png"%(title)
    
plt.savefig(save_path, dpi=800, bbox_inches="tight")

try:
    print(title)
    data=t[(t["min_45bin"]!="after_session") & (t["min_45bin"]!="before_session") & (t["gander"]=="female")]
    aovrm = AnovaRM(data, 'duration', 'id_cohort', within=['min_45bin', 'blind_or_not'])
    res = aovrm.fit()

    print(res)
except ValueError as e:
    print(str(e))


# ### plot 15min bin blind vs oriented time

# In[12]:


#generate 15min bin grouped file for each session
path_list=[coh1_male_path, coh1_female_path, coh2_path]
coh=["coh1", "coh1", "coh2"]
delte_str=[4,12,12]
code_file_name=[coh1_male_code["file_name"].tolist(), coh1_female_code["file_name"].tolist(), coh2_tot_code["vfa_csv_2"].tolist()]
nfps=[coh1_male_code["n_fps"].tolist(), coh1_female_code["n_fps"].tolist(), coh2_tot_code["n_fps"].tolist()]

for i,m, h,j, d in zip(path_list,code_file_name, nfps, coh, delte_str) :
    os.chdir(i)
    
    save_path=i[:-7]+"\\15_oriented"
    if os.path.isdir(save_path) == False:
        os.mkdir(save_path)
        
    for n, b in zip(m, h):
        df=pd.read_csv(n, index_col=[0]) 
        df.dropna(subset=["index", "total_visual"], inplace=True)
        df=df.groupby(by=["15min_bin", "blind_or_not"], as_index=False).count()
        df["cohort"]=j
        df["dt_id"]=n[:-d]
        df["animal_id"]=df["dt_id"][1][11:]
        df["id_cohort"]=df["animal_id"]+ "-" + df["cohort"]                    
        df["duration"]=df['id']/b
        df.to_csv(save_path + "\\%s"%(n), index=False)


# In[25]:


path_list=[coh1_male_path, coh1_female_path, coh2_path]
df_list=[]

for i in path_list:
    file_list=glob.glob(i[:-7]+"\\15_oriented\\*.csv")
    df_list=df_list+file_list
    
#combine 15min bin df from three parts

df=pd.DataFrame()
for i in df_list:
    df_=pd.read_csv(i)
    if len(df_)<8:
        print(i)
    df=pd.concat([df, df_], ignore_index=True)
df['group_id']=df['animal_id'].apply(lambda x:x[0])
df['gander']=np.where(df['group_id'].isin(["W", "X", "Y","Z"]), "male", "female")

gander=["male", "female"]
color_list=['royalblue', 'coral']
for i, p in zip(gander, color_list) :
    title="%s_tot_baseline_15min_bin"%(i)
    t=df.groupby(by=["id_cohort", "15min_bin", "blind_or_not", "gander"], as_index=False)['duration'].mean()
    t=t.rename(columns={"15min_bin": "min_15bin"})
    
    plot=sns.pointplot(data=t[(t["min_15bin"]!="after_session") & (t["min_15bin"]!="before_session") & (t["gander"]==i)], 
                   x="min_15bin", y="duration", hue="blind_or_not",  hue_order=["oriented", "blind"], palette=[p, p], 
                       linestyles=['-', '-.'], markers=["o", "x"], 
                       capsize=.01, errorbar=('ci', 68))
    plot.set(xlabel="")
    plot.set(ylabel="total_duration")
    plt.title(title)
    plt.xticks(rotation = 15, ha='right')
#     h = plt.gca().get_lines()
#     lg = plt.legend(handles=h, labels=['oriented', 'blind'], loc='best')
    plt.legend(bbox_to_anchor=(1.25, 1),borderaxespad=0)
    
    save_path=r"D:\CPT_TOT_manuscript\data processing\fig\vfa_15"
    
    if os.path.isdir(save_path ) == False:
        os.mkdir(save_path)
    
    plt.savefig(save_path+"\%s.png"%(title), dpi=800, bbox_inches="tight")
    plt.close()
    
    try:
        print(title)
        data=t[(t["min_15bin"]!="after_session") & (t["min_15bin"]!="before_session") & (t["gander"]==i)]
        aovrm = AnovaRM(data, 'duration', 'id_cohort', within=['min_15bin', 'blind_or_not'])
        res = aovrm.fit()

        print(res)
        
        print("try pingouin method")
        aov=pg.rm_anova(dv='duration', within=['min_15bin', 'blind_or_not'], subject='id_cohort', data=data)
        
        pg.print_table(aov)
        
        # Optional post-hoc tests
        print("post hoc......")
        psthc=pg.pairwise_tests(dv='duration', within=['min_15bin', 'blind_or_not'], subject='id_cohort', data=data)
        print(psthc)
        
    except ValueError as e:
        
        print(str(e))
        print("---------------------------------------------------")
        print("using mixed effect model...")
        print("title")
        print("with interation...")
        t_=t[(t["min_15bin"]!="after_session") & (t["min_15bin"]!="before_session") & (t["gander"]==i)].copy()
        model=smf.mixedlm("duration ~ C(min_15bin, Treatment('15bin_1')) + C(blind_or_not) + C(min_15bin, Treatment('15bin_1')):C(blind_or_not)", 
                          data=t_, groups=t_['id_cohort']).fit()
        display(model.summary())
    
        print("                                                   ")
        print("excluding interaction...")
        model=smf.mixedlm("duration ~ C(min_15bin, Treatment('15bin_1')) + C(blind_or_not)", 
                          data=t_, groups=t_['id_cohort']).fit()    
        display(model.summary())
    
    
        print("KDE plot generating......")
        fig = plt.figure(figsize = (16, 9))
        ax = sns.distplot(model.resid, hist = False, kde_kws = {"shade" : True, "lw": 1}, fit = stats.norm)
        ax.set_title(title + "KDE Plot of Model Residuals (Blue) and Normal Distribution (Black)")
        ax.set_xlabel("Residuals")
        save_stat=save_path+"\stat"
        if os.path.isdir(save_stat ) == False:
            os.mkdir(save_stat)
        
        plt.savefig(save_stat+"\%s_KDE.png"%title)
        plt.close()
        print("KDE plot saved")
        
        print("                                          ")
        
        print("Q-Q plot generating......")
        fig = plt.figure(figsize = (16, 9))
        ax = fig.add_subplot(111)
        sm.qqplot(model.resid, dist = stats.norm, line = 's', ax = ax)
        ax.set_title(file_name +"Q-Q Plot")
        plt.savefig(save_stat+"\%s_Q-Q_plot.png"%title)
        plt.close()
    
        print("KDE plot saved")
    
        labels = ["Statistic", "p-value"]

        norm_res = stats.shapiro(model.resid)

        for key, val in dict(zip(labels, norm_res)).items():
            print(key, val)


# In[26]:


# combine male and female into same plot, seperate oriented and blind
path_list=[coh1_male_path, coh1_female_path, coh2_path]
df_list=[]

for i in path_list:
    file_list=glob.glob(i[:-7]+"\\15_oriented\\*.csv")
    df_list=df_list+file_list
    
#combine 15min bin df from three parts

df=pd.DataFrame()
for i in df_list:
    df_=pd.read_csv(i)
    if len(df_)<8:
        print(i)
    df=pd.concat([df, df_], ignore_index=True)
df['group_id']=df['animal_id'].apply(lambda x:x[0])
df['gander']=np.where(df['group_id'].isin(["W", "X", "Y","Z"]), "male", "female")

oritation=["oriented", "blind"]
line_styles=['-', '-.']
markers=["o", "x"]
for i,j,p in zip(oritation, line_styles, markers):
    title="%s_tot_baseline_15min_bin"%(i)
    t=df.groupby(by=["id_cohort", "15min_bin", "blind_or_not", "gander"], as_index=False)['duration'].mean()
    t=t.rename(columns={"15min_bin": "min_15bin"})
    plot=sns.pointplot(data=t[(t["min_15bin"]!="after_session") & (t["min_15bin"]!="before_session") & (t["blind_or_not"]==i)], 
                   x="min_15bin", y="duration", hue="gander", linestyles=[j, j], markers=[p,p], hue_order=['male', 'female'], palette=['royalblue', 'coral'],
                       capsize=.01, errorbar=('ci', 68))
    plot.set(xlabel="")
    plot.set(ylabel="total_duration")
    plt.title(title)
    plt.xticks(rotation = 15, ha='right')
    
    plt.legend(bbox_to_anchor=(1.25, 1),borderaxespad=0)
    
    save_path=r"D:\CPT_TOT_manuscript\data processing\fig\vfa_15"
    
    if os.path.isdir(save_path ) == False:
        os.mkdir(save_path)
    
    plt.savefig(save_path+"\%s.png"%(title), dpi=800, bbox_inches="tight")
    plt.close()
    
    try:
        print(title)
        data=t[(t["min_15bin"]!="after_session") & (t["min_15bin"]!="before_session") & (t["blind_or_not"]==i)]
        aovrm = AnovaRM(data, 'duration', 'id_cohort', within=['min_15bin', 'gander'])
        res = aovrm.fit()
        print(res)
        
    except ValueError as e:
        print(str(e))
        print("---------------------------------------------------")
        print("using mixed effect model...")
        print("title")
        print("with interation...")
        t_=t[(t["min_15bin"]!="after_session") & (t["min_15bin"]!="before_session") & (t["blind_or_not"]==i)].copy()
        model=smf.mixedlm("duration ~ C(min_15bin, Treatment('15bin_1')) + C(gander) + C(min_15bin, Treatment('15bin_1')):C(gander)", 
                          data = t_, groups=t_['id_cohort']).fit()
        display(model.summary())

        print("                                                   ")
        print("excluding interaction...")

        model=smf.mixedlm("duration ~ C(min_15bin, Treatment('15bin_1')) + C(gander)", 
                          data =t_, groups=t_['id_cohort']).fit()    
        display(model.summary())


        print("KDE plot generating......")
        fig = plt.figure(figsize = (16, 9))
        ax =  sns.distplot(model.resid, kde=True, kde_kws = {"shade" : True, "lw": 1}, fit = stats.norm)

        ax.set_title(title + "KDE Plot of Model Residuals (Blue) and Normal Distribution (Black)")
        ax.set_xlabel("Residuals")
        save_stat=save_path+"\stat"
        if os.path.isdir(save_stat ) == False:
            os.mkdir(save_stat)

        plt.savefig(save_stat+"\%s_KDE.png"%title)
        plt.close()
        print("KDE plot saved")

        print("                                          ")

        print("Q-Q plot generating......")
        fig = plt.figure(figsize = (16, 9))
        ax = fig.add_subplot(111)
        sm.qqplot(model.resid, dist = stats.norm, line = 's', ax = ax)
        ax.set_title(title +"Q-Q Plot")
        plt.savefig(save_stat+"\%s_Q-Q_plot.png"%title)
        plt.close()

        print("KDE plot saved")

        labels = ["Statistic", "p-value"]

        norm_res = stats.shapiro(model.resid)

        for key, val in dict(zip(labels, norm_res)).items():
            print(key, val)


# In[ ]:




