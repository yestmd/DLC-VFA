#!/usr/bin/env python
# coding: utf-8

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


# In[37]:


from scipy import stats


# In[39]:


import statsmodels.api as sm


# In[2]:


import pingouin as pg


# In[3]:


#sync_temp2 path
coh1_male_path="D:\\CPT_TOT_DS_Stage_3\\CPT_DS_TOT_Visual_analysis\\vfa male tot_based on 1000000 model_fixed\\dlc_results_namechanged\\results\\temp_2"
coh1_female_path="D:\\CPT_TOT_DS_Stage_3\\CPT_DS_TOT_Visual_analysis\\vfa female tot_based on 1030000 model\\results\\temp_2"
coh2_path="D:\\CPT_TOT_DS_Stage_3_2nd\\VFA_baseline\\DLC_output\\results\\temp_2"


# ### plot 45min bin vfa_event

# In[4]:


coh1_male_code=pd.read_csv(r"D:\CPT_TOT_DS_Stage_3\CPT_DS_TOT_Visual_analysis\file_frame_info_male.csv", index_col=False)
coh1_female_code=pd.read_csv(r"D:\CPT_TOT_DS_Stage_3\CPT_DS_TOT_Visual_analysis\file_frame_info_female.csv", index_col=False)
coh2_tot_code=pd.read_csv(r"D:\CPT_TOT_DS_Stage_3_2nd\VFA_baseline\coh2_tot_code.csv", index_col=False)


# In[5]:


coh1_male_code.head(2)


# In[6]:


coh1_female_code.head(2)


# In[7]:


coh2_tot_code.head(2)


# ### calculate 45min bin

# In[8]:


path_list=[coh1_male_path, coh1_female_path, coh2_path]
coh=["coh1", "coh1", "coh2"]
delte_str=[4,12,12]
code_file_name=[coh1_male_code["file_name"].tolist(), coh1_female_code["file_name"].tolist(), coh2_tot_code["vfa_csv_2"].tolist()]
nfps=[coh1_male_code["n_fps"].tolist(), coh1_female_code["n_fps"].tolist(), coh2_tot_code["n_fps"].tolist()]

for i,m, h,j, d in zip(path_list,code_file_name, nfps, coh, delte_str) :
    os.chdir(i)
    
    save_path=i[:-7]+"\\45_oriented_event"
    if os.path.isdir(save_path) == False:
        os.mkdir(save_path)
        
    for n, b in zip(m, h):
        df=pd.read_csv(n, index_col=[0]) 
        df=df[df["count"]==1]
        df=df.groupby(by=["45min_bin", "event", "ori"], as_index=False).count()
        df["cohort"]=j
        df["dt_id"]=n[:-d]
        df["animal_id"]=df["dt_id"][1][11:]
        df["id_cohort"]=df["animal_id"]+ "-" + df["cohort"]                    
        df.to_csv(save_path + "\\%s"%(n), index=False)


# In[9]:


#combine 45min bin df from three parts

path_list=[coh1_male_path, coh1_female_path, coh2_path]
df_list=[]

for i in path_list:
    file_list=glob.glob(i[:-7]+"\\45_oriented_event\\*.csv")
    df_list=df_list+file_list

df45=pd.DataFrame()
for i in df_list:
    df_=pd.read_csv(i)
#     if len(df_)<8:
#         print(i)
    df45=pd.concat([df45, df_], ignore_index=True)
    
df45['group_id']=df45['animal_id'].apply(lambda x:x[0])
df45['gander']=np.where(df45['group_id'].isin(["W", "X", "Y","Z"]), "male", "female")


# In[10]:


err=df45[(df45["event"]=="Hits - #") & (df45["ori"]=="blind")].copy()
err=err[["45min_bin", "event", "ori", "dt_id", "id_cohort", "count"]].sort_values(by=["count", "id_cohort"])
display(err.shape)
err.sort_values(by='dt_id')


# #### check undeterminded frame

# In[11]:


err_undete=df45[(df45["event"]=="Mistakes - #") & (df45["ori"]=="undetermined")].copy()
err_undete=err_undete[["45min_bin", "event", "ori", "dt_id", "id_cohort", "count"]].sort_values(by=["count", "id_cohort"])
display(err_undete.shape)
err_undete.sort_values(by='dt_id')


# In[12]:


err_undete.to_csv(r"D:\CPT_TOT_manuscript\err_undeter_after.csv")


# ### calculate d'

# In[8]:


def hr_calcu(hit, miss):
    hr_rate = hit/(hit+miss)
    return hr_rate

def fr_calcu(mistake, correct_rejection):
    fr_rate = mistake/(mistake+correct_rejection)
    return fr_rate

def d_calcu(hit, miss, mistake, correct_rejection):
    d=norm.ppf(hr_calcu(hit, miss))-norm.ppf(fr_calcu(mistake, correct_rejection))
    return d  

def si_calcu(hit, miss, mistake, correct_rejection):
    si=(hr_calcu(hit, miss)-fr_calcu(mistake, correct_rejection))/(2*(hr_calcu(hit, miss)+fr_calcu(mistake, correct_rejection))-pow((hr_calcu(hit, miss)+fr_calcu(mistake, correct_rejection)), 2))
    return si

def c_calcu(hit, miss, mistake, correct_rejection):
    c=-(norm.ppf(hr_calcu(hit, miss))+norm.ppf(fr_calcu(mistake, correct_rejection)))/2
    return c

def ri_calcu(hit, miss, mistake, correct_rejection):
    ri=(hr_calcu(hit, miss)+fr_calcu(mistake, correct_rejection)-1)/(1-pow((hr_calcu(hit, miss)-fr_calcu(mistake, correct_rejection)), 2))
    return ri


# In[14]:


h=pd.pivot_table(df45, values='count', columns=['event'], index=['45min_bin', 'ori', 'id_cohort', 'gander'], aggfunc=np.mean)
h=h.reset_index()


# In[15]:


h[h["ori"]=='blind']


# In[16]:


df=h
try:
    df.loc[:, 'Hit_Rate'] =hr_calcu(df.loc[:, 'Hits - #'],df.loc[:, 'Misses - #'])
    df.loc[:, 'False_Alarm_Rate'] =fr_calcu(df.loc[:, 'Mistakes - #'], df.loc[:, 'Correct Rejections - #'])
    df.loc[:, 'c'] =c_calcu(df.loc[:, 'Hits - #'], df.loc[:, 'Misses - #'], df.loc[:, 'Mistakes - #'], df.loc[:, 'Correct Rejections - #'])
    df.loc[:, 'd'] =d_calcu(df.loc[:, 'Hits - #'], df.loc[:, 'Misses - #'], df.loc[:, 'Mistakes - #'], df.loc[:, 'Correct Rejections - #'])
    df.loc[:, 'si'] =si_calcu(df.loc[:, 'Hits - #'], df.loc[:, 'Misses - #'], df.loc[:, 'Mistakes - #'], df.loc[:, 'Correct Rejections - #'])
    df.loc[:, 'ri'] =ri_calcu(df.loc[:, 'Hits - #'], df.loc[:, 'Misses - #'], df.loc[:, 'Mistakes - #'], df.loc[:, 'Correct Rejections - #'])
    #df.loc[:, 'rp']=df.loc[:, 'Center Screen Touch - Time']/df.loc[:, 'Stim Onset - End ITI']
except ValueError as e:
    print(str(e))


# In[17]:


df45=df.melt(id_vars=["id_cohort", "45min_bin", "ori",'gander'], var_name="event", value_name="value")


# In[18]:


np.isinf(df45['value']).values.sum()


# In[19]:


df45[df45['value'].isna()]


# In[20]:


#err.sort_values(by='dt_id').to_csv(r"D:\CPT_TOT_manuscript\err_dismatch_after.csv")


# In[21]:


df45=df45.groupby(by=["id_cohort", "45min_bin", "event", "ori",'gander'], as_index=False)["value"].mean()


# In[22]:


df45["event"].unique()


# In[23]:


df45[df45["ori"]=='blind']


# ### plot 45min ori 

# In[24]:


#generate oriented vs blind plot
gander=['male', 'female']
event=['Correct Rejections - #', 'Hits - #','Misses - #', 'Mistakes - #', 'Hit_Rate', 'False_Alarm_Rate', 'c', 'd', 'ri', 'si']

for i in gander:
    for j in event:
        df=df45[(df45["event"]==j) & (df45["gander"]==i) & (df45["ori"]!="undetermined")]
        df=df.rename(columns={"45min_bin": "min_45bin"})
        plot=sns.barplot(data=df, x="ori", y="value", hue='min_45bin',
            order=["oriented", "blind"], capsize=.01, errorbar=('ci', 68))
        title=j+"-"+i
        plt.title(title)
        plot.set(xlabel="")
        plot.set(ylabel="")
        plt.legend(bbox_to_anchor=(1.26, 1),borderaxespad=0)
        save_path=r"D:\CPT_TOT_manuscript\data processing\fig\vfa_45_event_ori"
    
        if os.path.isdir(save_path ) == False:
            os.mkdir(save_path)
    
        plt.savefig(save_path+"\%s.png"%(title), dpi=800, bbox_inches="tight")
        plt.close()
        
        try:
            print(title)
            print("oriented")
            aovrm = AnovaRM(df[df["ori"]=="oriented"], 'value', 'id_cohort', within=['min_45bin','ori'])
            res = aovrm.fit()
            print(res)
            
            print("blind")
            aovrm = AnovaRM(df[df["ori"]=="blind"], 'value', 'id_cohort', within=['min_45bin','ori'])
            res = aovrm.fit()
            print(res)
            
            psthc=pg.pairwise_tests(dv='value', within=['min_45bin','ori'], subject='id_cohort', data=df[df["ori"]=="oriented"])
            print(psthc)
        #pg.print_table(psthc, tablefmt="pretty")
        except ValueError as e:
            print(str(e))


# ### calculate 15min bin

# In[5]:


path_list=[coh1_male_path, coh1_female_path, coh2_path]
coh=["coh1", "coh1", "coh2"]
delte_str=[4,12,12]
code_file_name=[coh1_male_code["file_name"].tolist(), coh1_female_code["file_name"].tolist(), coh2_tot_code["vfa_csv_2"].tolist()]
nfps=[coh1_male_code["n_fps"].tolist(), coh1_female_code["n_fps"].tolist(), coh2_tot_code["n_fps"].tolist()]

for i,m, h,j, d in zip(path_list,code_file_name, nfps, coh, delte_str) :
    os.chdir(i)
    
    save_path=i[:-7]+"\\15_oriented_event"
    if os.path.isdir(save_path) == False:
        os.mkdir(save_path)
        
    for n, b in zip(m, h):
        df=pd.read_csv(n, index_col=[0]) 
        df=df[df["count"]==1]
        df=df.groupby(by=["15min_bin", "event", "ori"], as_index=False).count()
        df["cohort"]=j
        df["dt_id"]=n[:-d]
        df["animal_id"]=df["dt_id"][1][11:]
        df["id_cohort"]=df["animal_id"]+ "-" + df["cohort"]                    
        df.to_csv(save_path + "\\%s"%(n), index=False)


# In[6]:


#combine 45min bin df from three parts
path_list=[coh1_male_path, coh1_female_path, coh2_path]
df_list=[]

for i in path_list:
    file_list=glob.glob(i[:-7]+"\\15_oriented_event\\*.csv")
    df_list=df_list+file_list
    

df15=pd.DataFrame()
for i in df_list:
    df_=pd.read_csv(i)
#     if len(df_)<8:
#         print(i)
    df15=pd.concat([df15, df_], ignore_index=True)
df15['group_id']=df15['animal_id'].apply(lambda x:x[0])
df15['gander']=np.where(df15['group_id'].isin(["W", "X", "Y","Z"]), "male", "female")


# In[9]:


h=pd.pivot_table(df15, values='count', columns=['event'], index=['15min_bin', 'ori', 'id_cohort', 'gander'], aggfunc=np.mean)
h=h.reset_index()

df=h
try:
    df.loc[:, 'Hit_Rate'] =hr_calcu(df.loc[:, 'Hits - #'],df.loc[:, 'Misses - #'])
    df.loc[:, 'False_Alarm_Rate'] =fr_calcu(df.loc[:, 'Mistakes - #'], df.loc[:, 'Correct Rejections - #'])
    df.loc[:, 'c'] =c_calcu(df.loc[:, 'Hits - #'], df.loc[:, 'Misses - #'], df.loc[:, 'Mistakes - #'], df.loc[:, 'Correct Rejections - #'])
    df.loc[:, 'd'] =d_calcu(df.loc[:, 'Hits - #'], df.loc[:, 'Misses - #'], df.loc[:, 'Mistakes - #'], df.loc[:, 'Correct Rejections - #'])
    df.loc[:, 'si'] =si_calcu(df.loc[:, 'Hits - #'], df.loc[:, 'Misses - #'], df.loc[:, 'Mistakes - #'], df.loc[:, 'Correct Rejections - #'])
    df.loc[:, 'ri'] =ri_calcu(df.loc[:, 'Hits - #'], df.loc[:, 'Misses - #'], df.loc[:, 'Mistakes - #'], df.loc[:, 'Correct Rejections - #'])
    #df.loc[:, 'rp']=df.loc[:, 'Center Screen Touch - Time']/df.loc[:, 'Stim Onset - End ITI']
except ValueError as e:
    print(str(e))

df15=df.melt(id_vars=["id_cohort", "15min_bin", "ori",'gander'], var_name="event", value_name="value")


# In[10]:


h


# In[11]:


df15=df15.groupby(by=["id_cohort", "15min_bin", "event", "ori",'gander'], as_index=False)["value"].mean()


# In[30]:


gander=['male', 'female']
event=['Correct Rejections - #', 'Hits - #','Misses - #', 'Mistakes - #', 'Hit_Rate', 'False_Alarm_Rate', 'c', 'd', 'ri', 'si']

for i in gander:
    for j in event:
        df=df15[(df15["event"]==j) & (df15["gander"]==i) & (df15["ori"]!="undetermined")]
        df=df.rename(columns={"15min_bin": "min_15bin"})
        plot=sns.barplot(data=df, x="ori", y="value", hue='min_15bin',
            order=["oriented", "blind"], capsize=.01, errorbar=('ci', 68))
        title=j+"-"+i
        plt.title(title)
        plot.set(xlabel="")
        plot.set(ylabel="")
        plt.legend(bbox_to_anchor=(1.26, 1),borderaxespad=0)
        save_path=r"D:\CPT_TOT_manuscript\data processing\fig\vfa_15_event_ori"
    
        if os.path.isdir(save_path ) == False:
            os.mkdir(save_path)
    
        plt.savefig(save_path+"\%s.png"%(title), dpi=800, bbox_inches="tight")
        plt.close()
        
        try:
            print(title)
            print("oriented")
            aovrm = AnovaRM(df[df["ori"]=="oriented"], 'value', 'id_cohort', within=['min_15bin','ori'])
            res = aovrm.fit()
            print(res)
            
            print("blind")
            aovrm = AnovaRM(df[df["ori"]=="blind"], 'value', 'id_cohort', within=['min_15bin','ori'])
            res = aovrm.fit()
            print(res)
            
            psthc=pg.pairwise_tests(dv='value', within=['min_15bin','ori'], subject='id_cohort', data=df[df["ori"]=="oriented"])
            print(psthc)
        #pg.print_table(psthc, tablefmt="pretty")
        except ValueError as e:
            print(str(e))


# ### plot total trial orientated-45min bin

# In[31]:


total45=df45[df45["event"].isin(['Correct Rejections - #', 'Hits - #','Misses - #', 'Mistakes - #'])]
total45=total45.groupby(by=["gander", "45min_bin", "id_cohort", "ori"], as_index=False).sum()
total45=total45.rename(columns={"45min_bin": "bin_45min"})

gander=['male', 'female']
for i in gander:
    df=total45[(total45["gander"]==i) &  (total45["ori"]!="undetermined")]
    plot=sns.barplot(data=df, x="ori", y="value", hue='bin_45min',
            order=["oriented", "blind"], capsize=.01, ci=68)
    
    plot=sns.stripplot(data=df,x="ori", y="value", hue='bin_45min', hue_order=["first_45min","last_45min"],
                         order=["oriented", "blind"], dodge=True, jitter=True, color="gray", size=3.5, legend=False)
    
    title=i+ "-total-trial-orientation-45min"
    plt.title(title)
    plot.set(xlabel="")
    plot.set(ylabel="")
    plt.legend(bbox_to_anchor=(1.26, 1),borderaxespad=0)
    save_path=r"D:\CPT_TOT_manuscript\data processing\fig\vfa_45_event_ori\total_ori"
    
    if os.path.isdir(save_path ) == False:
        os.mkdir(save_path)
    
    plt.savefig(save_path+"\%s.png"%(title), dpi=800, bbox_inches="tight")
    plt.close()
    
    try:
        aovrm = AnovaRM(df, 'value', 'id_cohort', within=['bin_45min', 'ori'])
        res = aovrm.fit()

        print(res)
    except ValueError as e:
        print(str(e))
        


# In[32]:


total45


# ### plot total trial orientated-15min bin

# In[53]:


total15=df15[df15["event"].isin(['Correct Rejections - #', 'Hits - #','Misses - #', 'Mistakes - #'])]
total15=total15.groupby(by=["gander", "15min_bin", "id_cohort", "ori"], as_index=False).sum()
total15=total15.rename(columns={"15min_bin": "bin_15min"})

gander=['male', 'female']
color_list=['royalblue', 'coral']
for i, p in zip(gander,color_list):
    df=total15[(total15["gander"]==i) &  (total15["ori"]!="undetermined")]
    plot=sns.pointplot(data=df, x="bin_15min", y="value", hue="ori",  hue_order=["oriented", "blind"], palette=[p, p], 
                       linestyles=['-', '-.'], markers=["o", "x"], 
                       capsize=.01, errorbar=('ci', 68))
    title=i+ "-total-trial-orientation-15min"
    plt.title(title)
    plot.set(xlabel="")
    plot.set(ylabel="")
    plt.legend(bbox_to_anchor=(1.26, 1),borderaxespad=0)
    save_path=r"D:\CPT_TOT_manuscript\data processing\fig\vfa_15_event_ori\total_ori"
    
    if os.path.isdir(save_path ) == False:
        os.mkdir(save_path)
    
    plt.savefig(save_path+"\%s.png"%(title), dpi=800, bbox_inches="tight")
    plt.close()
    
    try:
        print(title)
        aovrm = AnovaRM(df, 'value', 'id_cohort', within=['bin_15min', 'ori'])
        res = aovrm.fit()
        print(res)
        
        print("try pingouin method")
        aov=pg.rm_anova(dv='value', within=['bin_15min', 'ori'], subject='id_cohort', data=df)
        
        pg.print_table(aov)
        
        # Optional post-hoc tests
        print("post hoc......")
        psthc=pg.pairwise_tests(dv='value', within=['bin_15min', 'ori'], subject='id_cohort', data=df)
        print(psthc)
        
    except ValueError as e:
        print(str(e))
        print("---------------------------------------------------")
        print("using mixed effect model...")
        print(title)
        print("with interation...")
        
        model=smf.mixedlm("value ~ C(bin_15min, Treatment('15bin_1')) + C(ori) + C(bin_15min, Treatment('15bin_1')):C(ori)", 
                          data = df, groups=df['id_cohort']).fit()
        display(model.summary())
    
        print("                                                   ")
        print("excluding interaction...")
        model=smf.mixedlm("value ~ C(bin_15min, Treatment('15bin_1')) + C(ori)",
                          data =df, groups=df['id_cohort']).fit()    
        display(model.summary())
    
    
        print("KDE plot generating......")
        fig = plt.figure(figsize = (16, 9))
        ax = sns.distplot(model.resid, hist = False, kde_kws = {"shade" : True, "lw": 1}, fit = stats.norm)

        ax.set_title(title + "KDE Plot of Model Residuals (Blue) and Normal Distribution (Black)")
        ax.set_xlabel("Residuals")
        save_path_2=r"D:\CPT_TOT_manuscript\data processing\fig\vfa_15_event_ori\total_ori\stat"
        
        if os.path.isdir(save_path_2 ) == False:
            os.mkdir(save_path_2)
        plt.savefig(save_path_2+"\%s_KDE.png"%title)
        plt.close()
        print("KDE plot saved")
        
        print("                                          ")
        
        print("Q-Q plot generating......")
        fig = plt.figure(figsize = (16, 9))
        ax = fig.add_subplot(111)
        sm.qqplot(model.resid, dist = stats.norm, line = 's', ax = ax)
        ax.set_title(file_name +"Q-Q Plot")
        plt.savefig(save_path_2+"\%s_Q-Q_plot.png"%title)
        plt.close()
    
        print("KDE plot saved")
    
        labels = ["Statistic", "p-value"]

        norm_res = stats.shapiro(model.resid)

        for key, val in dict(zip(labels, norm_res)).items():
            print(key, val)
        


# In[29]:


total15.head(2)


# ### plot oriented ratio

# In[45]:


h=total15.pivot(index=["gander", "bin_15min", "id_cohort"], columns='ori', values='value').reset_index().copy()
h["oriented trials ratio"]=h["oriented"]/(h["oriented"]+h["blind"])
total15=h.melt(id_vars=["bin_15min", "id_cohort",'gander'], var_name="ori", value_name="value")
h=total15[total15["ori"]=="oriented trials ratio"]

plot=sns.lineplot(data=h, x="bin_15min", y="value", hue="gander", style="gander")
title="trial-orientation-ratio-15min"
plt.title(title)
plot.set(xlabel="")
plot.set(ylabel="")
plt.legend(bbox_to_anchor=(1.26, 1),borderaxespad=0)
save_path=r"D:\CPT_TOT_manuscript\data processing\fig\vfa_15_event_ori\total_ori"

if os.path.isdir(save_path ) == False:
    os.mkdir(save_path)

plt.savefig(save_path+"\%s.png"%(title), dpi=800, bbox_inches="tight")


print("male"+"-"+title)
h_=h[h["gander"]=="male"]
aovrm = AnovaRM(h_, 'value', 'id_cohort', within=['bin_15min'])
res = aovrm.fit()
print(res)

print("female"+"-"+title)
h_=h[h["gander"]=="female"]
aovrm = AnovaRM(h_, 'value', 'id_cohort', within=['bin_15min'])
res = aovrm.fit()
print(res)


print("---------------------------------------------------")
print("using mixed effect model...")
print(title)
print("with interation...")

model=smf.mixedlm("value ~ C(bin_15min, Treatment('15bin_1')) + C(gander) + C(bin_15min, Treatment('15bin_1')):C(gander)", 
                  data = h, groups=h['id_cohort']).fit()
display(model.summary())

print("                                                   ")
print("excluding interaction...")
model=smf.mixedlm("value ~ C(bin_15min, Treatment('15bin_1')) + C(gander)",
                  data =h, groups=h['id_cohort']).fit()    
display(model.summary())


print("KDE plot generating......")
fig = plt.figure(figsize = (16, 9))
ax = sns.distplot(model.resid, hist = False, kde_kws = {"shade" : True, "lw": 1}, fit = stats.norm)

ax.set_title(title + "KDE Plot of Model Residuals (Blue) and Normal Distribution (Black)")
ax.set_xlabel("Residuals")
save_path_2=r"D:\CPT_TOT_manuscript\data processing\fig\vfa_15_event_ori\total_ori\stat"

if os.path.isdir(save_path_2) == False:
    os.mkdir(save_path_2)
plt.savefig(save_path_2+"\%s_KDE.png"%title)
plt.close()
print("KDE plot saved")

print("                                          ")

print("Q-Q plot generating......")
fig = plt.figure(figsize = (16, 9))
ax = fig.add_subplot(111)
sm.qqplot(model.resid, dist = stats.norm, line = 's', ax = ax)
ax.set_title(title +"Q-Q Plot")
plt.savefig(save_path_2+"\%s_Q-Q_plot.png"%title)
plt.close()

print("KDE plot saved")

labels = ["Statistic", "p-value"]

norm_res = stats.shapiro(model.resid)

for key, val in dict(zip(labels, norm_res)).items():
    print(key, val)


# In[ ]:





# ### plot total orientation-45min bin

# In[34]:


total=total45.groupby(by=["gander", "id_cohort", "ori"], as_index=False).sum()


# In[35]:


gander=['male', 'female']
for i in gander:
    plot=sns.barplot(data=total[(total["gander"]==i) & (total["ori"]!="undetermined")], x="ori", y="value",
            order=["oriented", "blind"], capsize=.01, ci=68)
    plot=sns.stripplot(data=total[(total["gander"]==i) & (total["ori"]!="undetermined")],x="ori", y="value", 
                         order=["oriented", "blind"], dodge=True, color="gray", size=3.5)
    title=i+ "-total-trial-orientation"
    plt.title(title)
    plot.set(xlabel="")
    plot.set(ylabel="")
    plt.legend(bbox_to_anchor=(1.26, 1),borderaxespad=0)
    save_path=r"D:\CPT_TOT_manuscript\data processing\fig\vfa_45_event_ori\total_ori"
    
    if os.path.isdir(save_path ) == False:
        os.mkdir(save_path)
    
    plt.savefig(save_path+"\%s.png"%(title), dpi=800, bbox_inches="tight")
    plt.close()
    


# In[ ]:




