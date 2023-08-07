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


import pingouin as pg


# In[3]:


#sync_temp2 path
coh1_male_path="D:\\CPT_TOT_DS_Stage_3\\CPT_DS_TOT_Visual_analysis\\vfa male tot_based on 1000000 model_fixed\\dlc_results_namechanged\\results\\temp_2"
coh1_female_path="D:\\CPT_TOT_DS_Stage_3\\CPT_DS_TOT_Visual_analysis\\vfa female tot_based on 1030000 model\\results\\temp_2"
coh2_path="D:\\CPT_TOT_DS_Stage_3_2nd\\VFA_baseline\\DLC_output\\results\\temp_2"


# ### plot 45min bin blind vs oriented time

# In[4]:


coh1_male_code=pd.read_csv(r"D:\CPT_TOT_DS_Stage_3\CPT_DS_TOT_Visual_analysis\file_frame_info_male.csv", index_col=False)
coh1_female_code=pd.read_csv(r"D:\CPT_TOT_DS_Stage_3\CPT_DS_TOT_Visual_analysis\file_frame_info_female.csv", index_col=False)
coh2_tot_code=pd.read_csv(r"D:\CPT_TOT_DS_Stage_3_2nd\VFA_baseline\coh2_tot_code.csv", index_col=False)


# In[5]:


#generate 45min and 15min bin group df for each session
path_list=[coh1_male_path, coh1_female_path, coh2_path]
coh=["coh1", "coh1", "coh2"]
delte_str=[4,12,12]
code_file_name=[coh1_male_code["file_name"].tolist(), coh1_female_code["file_name"].tolist(), coh2_tot_code["vfa_csv_2"].tolist()]
nfps=[coh1_male_code["n_fps"].tolist(), coh1_female_code["n_fps"].tolist(), coh2_tot_code["n_fps"].tolist()]

for i,m, h,j, d in zip(path_list,code_file_name, nfps, coh, delte_str) :
    os.chdir(i)
    
    save_path=i[:-7]+"\\orient_detail"
    if os.path.isdir(save_path) == False:
        os.mkdir(save_path)
        
    for n, b in zip(m, h):
        df=pd.read_csv(n, index_col=[0])
        df.dropna(subset=["index", "total_visual"], inplace=True)
        df_=df[["blindright","lateralLeftright","lateralRightright", "frontalright", "total_visual", "45min_bin", "15min_bin", "cohort"]].copy()        
        
        df_45=df_.groupby(by=["45min_bin", "cohort"], as_index=False).mean()
        df_45["cohort"]=j
        df_45["dt_id"]=n[:-d]
        df_45["animal_id"]=df_45["dt_id"][1][11:]
        df_45["id_cohort"]=df_45["animal_id"]+ "-" + df_45["cohort"]                    
        df_45.to_csv(save_path + "\\%s_45min.csv"%(n), index=False)

        df_15=df_.groupby(by=["15min_bin", "cohort"], as_index=False).mean()
        df_15["cohort"]=j
        df_15["dt_id"]=n[:-d]
        df_15["animal_id"]=df_15["dt_id"][1][11:]
        df_15["id_cohort"]=df_15["animal_id"]+ "-" + df_15["cohort"]                    
        df_15.to_csv(save_path + "\\%s_15min.csv"%(n), index=False)


# ### plot 45min orentation detail

# In[6]:


#generate file list from three parts of data
path_list=[coh1_male_path, coh1_female_path, coh2_path]
df_list=[]

for i in path_list:
    file_list=glob.glob(i[:-7]+"\\orient_detail\\*45min*.csv")
    df_list=df_list+file_list
    
#combine 45min oriented bin df from three parts

df=pd.DataFrame()
for i in df_list:
    df_=pd.read_csv(i)
#     if len(df_)<8:
#         print(i)
    df=pd.concat([df, df_], ignore_index=True)
df['group_id']=df['animal_id'].apply(lambda x:x[0])
df['gander']=np.where(df['group_id'].isin(["W", "X", "Y", "Z"]), "male", "female")


# In[7]:


#prepare df for plot
t=df.groupby(by=["id_cohort", "45min_bin", "gander"], as_index=False).mean()
t=t.rename(columns={"45min_bin": "min_45bin"})
t_melt=t.melt(id_vars=["id_cohort", "min_45bin", "gander"], var_name="visual_direction", value_name="screen_ratio")


# ### "lateralLeftright","lateralRightright", "frontalright" comparesion

# In[8]:


#plot tot baseline and state
#split all direction into "blindright","lateralLeftright","lateralRightright", "frontalright"
#leave "total_visual" alone
for i in ["male","female"]:    
    df=t_melt[(t_melt["min_45bin"]!="after_session") & (t_melt["min_45bin"]!="before_session") 
              & (t_melt["gander"]==i) & (t_melt["visual_direction"]!="total_visual")
              & (t_melt["visual_direction"]!="blindright")].copy() 
              #]
    
    plot=sns.lineplot(data=df, x="min_45bin", y="screen_ratio", hue="visual_direction", style="visual_direction")
    title=i+"_dominate_oriented_field_in_TOT_45minbin"
    plot.set(xlabel="")
    plot.set(ylabel="screen_ratio_in_visual_field")
    plt.title(title)
    plt.xticks(rotation = 45, ha='right')
    plt.legend(bbox_to_anchor=(1.25, 1),borderaxespad=0)
    
    #create save path
    save_path=r"D:\CPT_TOT_manuscript\data processing\fig\orient_detail"       
    if os.path.isdir(save_path ) == False:
        os.mkdir(save_path)
    
    plt.savefig(save_path+"\%s.png"%(title), dpi=800, bbox_inches="tight")
    plt.close()
    
    #state
    try:
        print(title)
       
        aovrm = AnovaRM(df, 'screen_ratio', 'id_cohort', within=['min_45bin', 'visual_direction'])
        res = aovrm.fit()
        print(res)
        
        print("try pingouin method")
        aov=pg.rm_anova(dv='screen_ratio', within=['min_45bin','visual_direction'], subject='id_cohort', data=df)
        
        pg.print_table(aov)
        
        # Optional post-hoc tests
        psthc=pg.pairwise_tests(dv='screen_ratio', within=['min_45bin','visual_direction'], subject='id_cohort', data=df)
        print(psthc)
        #pg.print_table(psthc, tablefmt="pretty")
    except ValueError as e:
        print(str(e))


# ### "blindright","total_visual" comparesion

# In[9]:


#plot tot baseline and state
#split all direction into "blindright","lateralLeftright","lateralRightright", "frontalright"
#leave "total_visual" alone
for i in ["male","female"]:    
    df=t_melt[(t_melt["min_45bin"]!="after_session") & (t_melt["min_45bin"]!="before_session") 
              & (t_melt["gander"]==i) & (t_melt["visual_direction"]!="lateralLeftright") 
              & (t_melt["visual_direction"]!="lateralRightright")
              & (t_melt["visual_direction"]!="frontalright")].copy()
    
    plot=sns.lineplot(data=df, x="min_45bin", y="screen_ratio", hue="visual_direction", style="visual_direction")
    title=i+"_dominate_visual_field_in_TOT_45min_bin(blind vs total visual)"
    plot.set(xlabel="")
    plot.set(ylabel="screen_ratio_in_visual_field")
    plt.title(title)
    plt.xticks(rotation = 45, ha='right')
    plt.legend(bbox_to_anchor=(1.25, 1),borderaxespad=0)
    
    #create save path
    save_path=r"D:\CPT_TOT_manuscript\data processing\fig\orient_detail"       
    if os.path.isdir(save_path ) == False:
        os.mkdir(save_path)
    
    plt.savefig(save_path+"\%s.png"%(title), dpi=800, bbox_inches="tight")
    plt.close()
    
    #state
    try:
        print(title)
       
        aovrm = AnovaRM(df, 'screen_ratio', 'id_cohort', within=['min_45bin', 'visual_direction'])
        res = aovrm.fit()
        print(res)
        
        print("try pingouin method")
        aov=pg.rm_anova(dv='screen_ratio', within=['min_45bin','visual_direction'], subject='id_cohort', data=df)
        
        pg.print_table(aov)
        
        # Optional post-hoc tests
        psthc=pg.pairwise_tests(dv='screen_ratio', within=['min_45bin','visual_direction'], subject='id_cohort', data=df)
        print(psthc)
        #pg.print_table(psthc, tablefmt="pretty")
    except ValueError as e:
        print(str(e))


# ### plot 15min orentation detail

# In[10]:


#generate file list from three parts of data
path_list=[coh1_male_path, coh1_female_path, coh2_path]
df_list=[]

for i in path_list:
    file_list=glob.glob(i[:-7]+"\\orient_detail\\*15min*.csv")
    df_list=df_list+file_list
    
#combine 15min oriented bin df from three parts

df=pd.DataFrame()
for i in df_list:
    df_=pd.read_csv(i)
#     if len(df_)<8:
#         print(i)
    df=pd.concat([df, df_], ignore_index=True)
df['group_id']=df['animal_id'].apply(lambda x:x[0])
df['gander']=np.where(df['group_id'].isin(["W", "X", "Y", "Z"]), "male", "female")


# In[11]:


#prepare df for plot
t=df.groupby(by=["id_cohort", "15min_bin", "gander"], as_index=False).mean()
t=t.rename(columns={"15min_bin": "min_15bin"})
t_melt=t.melt(id_vars=["id_cohort", "min_15bin", "gander"], var_name="visual_direction", value_name="screen_ratio")


# In[12]:


#plot tot baseline and state
#split all direction into "blindright","lateralLeftright","lateralRightright", "frontalright"
#leave "total_visual" alone
for i in ["male","female"]:    
    df=t_melt[(t_melt["min_15bin"]!="after_session") & (t_melt["min_15bin"]!="before_session") 
              & (t_melt["gander"]==i) & (t_melt["visual_direction"]!="total_visual") 
              & (t_melt["visual_direction"]!="blindright")].copy() 

    plot=sns.lineplot(data=df, x="min_15bin", y="screen_ratio", hue="visual_direction", style="visual_direction")
    title=i+"_dominate_oriented_field_in_TOT_15minbin"
    plot.set(xlabel="")
    plot.set(ylabel="screen_ratio_in_visual_field")
    plt.title(title)
    plt.xticks(rotation = 15, ha='right')
    plt.legend(bbox_to_anchor=(1.25, 1),borderaxespad=0)
    
    #create save path
    save_path=r"D:\CPT_TOT_manuscript\data processing\fig\orient_detail"       
    if os.path.isdir(save_path ) == False:
        os.mkdir(save_path)
    
    plt.savefig(save_path+"\%s.png"%(title), dpi=800, bbox_inches="tight")
    plt.close()
    
    #state
    try:
        print(title)
       
        aovrm = AnovaRM(df, 'screen_ratio', 'id_cohort', within=['min_15bin', 'visual_direction'])
        res = aovrm.fit()
        print(res)
        
        print("try pingouin method")
        aov=pg.rm_anova(dv='screen_ratio', within=['min_15bin','visual_direction'], subject='id_cohort', data=df)
        
        pg.print_table(aov)
        
        # Optional post-hoc tests
        psthc=pg.pairwise_tests(dv='screen_ratio', within=['min_15bin','visual_direction'], subject='id_cohort', data=df)
        print(psthc)
        #pg.print_table(psthc, tablefmt="pretty")
    except ValueError as e:
        print(str(e))


# ### blind vs total visual

# In[13]:


#plot tot baseline and state
#split all direction into "blindright","lateralLeftright","lateralRightright", "frontalright"
#leave "total_visual" alone
for i in ["male","female"]:    
    df=t_melt[(t_melt["min_15bin"]!="after_session") & (t_melt["min_15bin"]!="before_session") 
              & (t_melt["gander"]==i) & (t_melt["visual_direction"]!="lateralLeftright") 
              & (t_melt["visual_direction"]!="lateralRightright")
              & (t_melt["visual_direction"]!="frontalright")].copy()
    
    plot=sns.lineplot(data=df, x="min_15bin", y="screen_ratio", hue="visual_direction", style="visual_direction")
    title=i+"_dominate_visual_field_in_TOT_15min_bin(blind vs total visual)"
    plot.set(xlabel="")
    plot.set(ylabel="screen_ratio_in_visual_field")
    plt.title(title)
    plt.xticks(rotation = 15, ha='right')
    plt.legend(bbox_to_anchor=(1.25, 1),borderaxespad=0)
    
    #create save path
    save_path=r"D:\CPT_TOT_manuscript\data processing\fig\orient_detail"       
    if os.path.isdir(save_path ) == False:
        os.mkdir(save_path)
    
    plt.savefig(save_path+"\%s.png"%(title), dpi=800, bbox_inches="tight")
    plt.close()
    
    #state
    try:
        print(title)
       
        aovrm = AnovaRM(df, 'screen_ratio', 'id_cohort', within=['min_15bin', 'visual_direction'])
        res = aovrm.fit()
        print(res)
        
        print("try pingouin method")
        aov=pg.rm_anova(dv='screen_ratio', within=['min_15bin','visual_direction'], subject='id_cohort', data=df)
        
        pg.print_table(aov)
        
        # Optional post-hoc tests
        psthc=pg.pairwise_tests(dv='screen_ratio', within=['min_15bin','visual_direction'], subject='id_cohort', data=df)
        print(psthc)
        #pg.print_table(psthc, tablefmt="pretty")
    except ValueError as e:
        print(str(e))


# In[ ]:




