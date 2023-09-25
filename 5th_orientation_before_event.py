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

import statsmodels.formula.api as smf
import scipy.stats as stats
import statsmodels.api as sm
import researchpy as rp

import pingouin as pg

#sync_temp2 path
coh1_male_path="D:\\CPT_TOT_DS_Stage_3\\CPT_DS_TOT_Visual_analysis\\vfa male tot_based on 1000000 model_fixed\\dlc_results_namechanged\\results\\temp_2"
coh1_female_path="D:\\CPT_TOT_DS_Stage_3\\CPT_DS_TOT_Visual_analysis\\vfa female tot_based on 1030000 model\\results\\temp_2"
coh2_path="D:\\CPT_TOT_DS_Stage_3_2nd\\VFA_baseline\\DLC_output\\results\\temp_2"

coh1_male_code=pd.read_csv(r"D:\CPT_TOT_DS_Stage_3\CPT_DS_TOT_Visual_analysis\file_frame_info_male.csv", index_col=False)
coh1_female_code=pd.read_csv(r"D:\CPT_TOT_DS_Stage_3\CPT_DS_TOT_Visual_analysis\file_frame_info_female.csv", index_col=False)
coh2_tot_code=pd.read_csv(r"D:\CPT_TOT_DS_Stage_3_2nd\VFA_baseline\coh2_tot_code.csv", index_col=False)

path_list=[coh1_male_path, coh1_female_path, coh2_path]
coh=["coh1", "coh1", "coh2"]
delte_str=[4,12,12]
code_file_name=[coh1_male_code["file_name"].tolist(), coh1_female_code["file_name"].tolist(), coh2_tot_code["vfa_csv_2"].tolist()]
nfps=[coh1_male_code["n_fps"].tolist(), coh1_female_code["n_fps"].tolist(), coh2_tot_code["n_fps"].tolist()]

event_list=["Hits - #", "Misses - #" , "Correct Rejections - #", "Mistakes - #" ]
bin_list=["15bin_1", "15bin_2","15bin_3","15bin_4","15bin_5","15bin_6"]

for i,m, h,j, d in zip(path_list,code_file_name, nfps, coh, delte_str) :
    os.chdir(i)
    
    save_path=i[:-7]+"\\ori_LH"
    if os.path.isdir(save_path) == False:
        os.mkdir(save_path)
    for n, b in zip(m, h):    
        file=pd.read_csv(n)
        
        df_=file[file["event"].isin(["Stim Onset - End ITI", "Hits - #", "Misses - #" , "Correct Rejections - #", "Mistakes - #" ])].copy()

        hit_index=[]
        miss_index=[]
        mistake_index=[]
        corr_reject_index=[]
        for e in range(len(df_)-1):
            if df_.iloc[e]["event"]=="Stim Onset - End ITI" and df_.iloc[e+1]["event"]=="Hits - #":
                hit_index.append(df_.iloc[[e,e+1]]["sec_cal"].to_list())
            if df_.iloc[e]["event"]=="Stim Onset - End ITI" and df_.iloc[e+1]["event"]=="Misses - #":
                miss_index.append(df_.iloc[[e,e+1]]["sec_cal"].to_list())
            if df_.iloc[e]["event"]=="Stim Onset - End ITI" and df_.iloc[e+1]["event"]=="Mistakes - #":
                mistake_index.append(df_.iloc[[e,e+1]]["sec_cal"].to_list())
            if df_.iloc[e]["event"]=="Stim Onset - End ITI" and df_.iloc[e+1]["event"]=="Correct Rejections - #":
                corr_reject_index.append(df_.iloc[[e,e+1]]["sec_cal"].to_list())

        hit=[]
        miss=[]
        mistake=[]
        corr_reject=[]        
        for l, c in zip([hit_index,miss_index,mistake_index,corr_reject_index],
                        ["hit", "miss", "mistake", "corr_reject"]):
            for f, a in zip(l, range(1, len(l)+1)):
                b=file[(file["sec_cal"] >= f[0]) & (file["sec_cal"] <= f[1])].copy()
                b["trial_number"]=a
                b["event_type"]=c
                
                if c=="hit":
                    hit.append(b)            
                elif c=="miss":
                    miss.append(b)
                elif c=="mistake":
                    mistake.append(b)
                else:
                    corr_reject.append(b)            
        hit=pd.concat(hit, ignore_index=True)
        miss=pd.concat(miss, ignore_index=True)
        mistake=pd.concat(mistake, ignore_index=True)
        corr_reject=pd.concat(corr_reject, ignore_index=True)
        
        #add each bin's trial number into df according to event type
        event_list=["Hits - #", "Misses - #" ,"Mistakes - #", "Correct Rejections - #"]
        df=[hit, miss, mistake, corr_reject]
        for x, y in zip(event_list, df): 
            bin_n=[]
            bin_condition=[(y["15min_bin"]=="15bin_1"), (y["15min_bin"]=="15bin_2"),
                       (y["15min_bin"]=="15bin_3"),(y["15min_bin"]=="15bin_4"),
                       (y["15min_bin"]=="15bin_5"),(y["15min_bin"]=="15bin_6")]
            for z in bin_condition:
                bin_n.append(len(y[(y["event"]==x) & z]))
            y["bin_trial_n"]=np.select(bin_condition, bin_n)
        #combine four different event into individual file
        df=pd.concat([hit, miss, mistake, corr_reject], ignore_index=True)
        df.to_csv(save_path + "\\%s"%(n), index=False)

path_list=[coh1_male_path, coh1_female_path, coh2_path]
coh=["coh1", "coh1", "coh2"]
delte_str=[4,12,12]
code_file_name=[coh1_male_code["file_name"].tolist(), coh1_female_code["file_name"].tolist(), coh2_tot_code["vfa_csv_2"].tolist()]
nfps=[coh1_male_code["n_fps"].tolist(), coh1_female_code["n_fps"].tolist(), coh2_tot_code["n_fps"].tolist()]

for i,m, h, j, d in zip(path_list,code_file_name, nfps, coh, delte_str) :
    os.chdir(i[:-7]+"\\ori_LH")
    
    save_path=i[:-7]+"\\15_ori_LH"
    if os.path.isdir(save_path) == False:
        os.mkdir(save_path)
        
    for n, b in zip(m, h):
        df=pd.read_csv(n)
        df=df.groupby(by=["15min_bin", "event_type", "blind_or_not"], as_index=False).agg(
            bin_frame_sum=('index', 'count'),
            bin_trial_n_sum=('bin_trial_n', 'mean'))
        df["cohort"]=j
        df["dt_id"]=n[:-d]
        df["animal_id"]=df["dt_id"][1][11:]
        df["id_cohort"]=df["animal_id"]+ "-" + df["cohort"]       
        df["duration"]=df['bin_frame_sum']/b
        df["ori_mean"]=df["duration"]/df["bin_trial_n_sum"]
        df.to_csv(save_path + "\\%s"%(n), index=False)

path_list=[coh1_male_path, coh1_female_path, coh2_path]
df_list=[]

for i in path_list:
    file_list=glob.glob(i[:-7]+"\\15_ori_LH\\*.csv")
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


event=["hit", "miss", "mistake", "corr_reject"]
gander=["male", "female"]
color_list=['royalblue', 'coral']
for i, p in zip(gander, color_list) :
    for j in event:        
        title="%s_%s_ori_mean_before_event_15min_bin"%(i, j)
        t=df[(df["gander"]==i) & (df["event_type"]==j)].copy()
        t=t.rename(columns={"15min_bin": "min_15bin"})
        dt=t[(t["min_15bin"]!="after_session") & (t["min_15bin"]!="before_session") & (t["gander"]==i)]
        plot=sns.pointplot(data=dt, 
                       x="min_15bin", y="ori_mean", hue="blind_or_not",  hue_order=["oriented", "blind"], palette=[p, p], 
                           linestyles=['-', '-.'], markers=["o", "x"], capsize=.01, errorbar=('ci', 68))
        plot.set(xlabel="")
        plot.set(ylabel="duration_mean")
        plt.title(title)
        plt.xticks(rotation = 15, ha='right')
    #     h = plt.gca().get_lines()
    #     lg = plt.legend(handles=h, labels=['oriented', 'blind'], loc='best')
        plt.legend(bbox_to_anchor=(1.25, 1),borderaxespad=0)

        save_path=r"D:\CPT_TOT_manuscript\data processing\fig\15_ori_LH"
        
        if os.path.isdir(save_path ) == False:
            os.mkdir(save_path)
        dt.to_csv(save_path+"\%s.csv"%(title), index=False)
        plt.savefig(save_path+"\%s.png"%(title), dpi=800, bbox_inches="tight")
        plt.close()

        try:
            print(title)
            data=t[(t["min_15bin"]!="after_session") & (t["min_15bin"]!="before_session") & (t["gander"]==i)]
            aovrm = AnovaRM(data, 'ori_mean', 'id_cohort', within=['min_15bin', 'blind_or_not'])
            res = aovrm.fit()

            print(res)

            print("try pingouin method")
            aov=pg.rm_anova(dv='ori_mean', within=['min_15bin', 'blind_or_not'], subject='id_cohort', data=data)

            pg.print_table(aov)

            # Optional post-hoc tests
            print("post hoc......")
            psthc=pg.pairwise_tests(dv='ori_mean', within=['min_15bin', 'blind_or_not'], subject='id_cohort', data=data)
            print(psthc)

        except ValueError as e:

            print(str(e))
            print("---------------------------------------------------")
            print("using mixed effect model...")
            print("title")
            print("with interation...")
            t_=t[(t["min_15bin"]!="after_session") & (t["min_15bin"]!="before_session") & (t["gander"]==i)].copy()
            model=smf.mixedlm("ori_mean ~ C(min_15bin, Treatment('15bin_1')) + C(blind_or_not) + C(min_15bin, Treatment('15bin_1')):C(blind_or_not)", 
                              data=t_, groups=t_['id_cohort']).fit()
            display(model.summary())

            print("                                                   ")
            print("excluding interaction...")
            model=smf.mixedlm("ori_mean ~ C(min_15bin, Treatment('15bin_1')) + C(blind_or_not)", 
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
            ax.set_title(title +"Q-Q Plot")
            plt.savefig(save_stat+"\%s_Q-Q_plot.png"%title)
            plt.close()

            print("KDE plot saved")

            labels = ["Statistic", "p-value"]

            norm_res = stats.shapiro(model.resid)

            for key, val in dict(zip(labels, norm_res)).items():
                print(key, val)

path_list=[coh1_male_path, coh1_female_path, coh2_path]
coh=["coh1", "coh1", "coh2"]
delte_str=[4,12,12]
code_file_name=[coh1_male_code["file_name"].tolist(), coh1_female_code["file_name"].tolist(), coh2_tot_code["vfa_csv_2"].tolist()]
nfps=[coh1_male_code["n_fps"].tolist(), coh1_female_code["n_fps"].tolist(), coh2_tot_code["n_fps"].tolist()]

for i,m, h, j, d in zip(path_list,code_file_name, nfps, coh, delte_str) :
    os.chdir(i[:-7]+"\\ori_LH")
#     
    save_path=i[:-7]+"\\15_ori_ratio_LH"
    if os.path.isdir(save_path) == False:
        os.mkdir(save_path)
        
    for n, b in zip(m, h):
        df=pd.read_csv(n)
        df=df.groupby(by=["15min_bin", "trial_number", "event_type", "blind_or_not"], as_index=False).agg(
            bin_frame_sum=('index', 'count'),bin_trial_n_sum=('bin_trial_n', 'mean'),
            rewardside=('LEFTCLOSE', 'sum'), screenside=('RIGHTCLOSE', 'sum'), 
            distance_move_trial=('distanceMoved', 'sum'))        
        df["cohort"]=j
        df["dt_id"]=n[:-d]
        df["animal_id"]=df["dt_id"][1][11:]
        df["id_cohort"]=df["animal_id"]+ "-" + df["cohort"]       
        df["duration"]=df['bin_frame_sum']/b
        df["ori_mean"]=df["duration"]/df["bin_trial_n_sum"]
        
        df.to_csv(save_path + "\\%s"%(n), index=False)

path_list=[coh1_male_path, coh1_female_path, coh2_path]
coh=["coh1", "coh1", "coh2"]
delte_str=[4,12,12]
code_file_name=[coh1_male_code["file_name"].tolist(), coh1_female_code["file_name"].tolist(), coh2_tot_code["vfa_csv_2"].tolist()]
nfps=[coh1_male_code["n_fps"].tolist(), coh1_female_code["n_fps"].tolist(), coh2_tot_code["n_fps"].tolist()]

for i,m, h, j, d in zip(path_list,code_file_name, nfps, coh, delte_str) :
    os.chdir(i[:-7]+"\\15_ori_ratio_LH")
    
    save_path=i[:-7]+"\\15_loca_ori_LH"
    if os.path.isdir(save_path) == False:
        os.mkdir(save_path)
        
    for n, b in zip(m, h):
        df=pd.read_csv(n)
        df_=pd.pivot_table(df, index=["15min_bin", "event_type", "trial_number", "dt_id", "animal_id", "id_cohort"], 
                          columns=["blind_or_not"], values=["rewardside", "screenside"]).reset_index().copy()
        df_.columns=['_'.join(col) for col in df_.columns.values]
        df_["rewardside_blind"] = df_["rewardside_blind"].fillna(0)
        df_["rewardside_oriented"] = df_["rewardside_oriented"].fillna(0)
        df_["screenside_blind"] = df_["screenside_blind"].fillna(0)
        df_["screenside_oriented"] = df_["screenside_oriented"].fillna(0)
        df_=df_.melt(id_vars=["15min_bin_", "event_type_", "trial_number_", "dt_id_", "animal_id_", "id_cohort_"], var_name="loca_ori", ignore_index=True)
        df_["loca_ori_duration"]=df_["value"]/b
        
        df_.to_csv(save_path + "\\%s"%(n), index=False)



df_

def normal(mean, std, histmax=False, color="black"):
    x = np.linspace(mean-4*std, mean+4*std, 200)
    p = stats.norm.pdf(x, mean, std)
    if histmax:
        p = p*histmax/max(p)
    z = plt.plot(x, p, color, linewidth=2)

path_list=[coh1_male_path, coh1_female_path, coh2_path]
df_list=[]

for i in path_list:
    file_list=glob.glob(i[:-7]+"\\15_loca_ori_LH\\*.csv")
    df_list=df_list+file_list
    
#combine 15min bin df from three parts

df=pd.DataFrame()
for i in df_list:
    df_=pd.read_csv(i)
    if len(df_)<8:
        print(i)
    df=pd.concat([df, df_], ignore_index=True)
df['group_id']=df['animal_id_'].apply(lambda x:x[0])
df['gander']=np.where(df['group_id'].isin(["W", "X", "Y","Z"]), "male", "female")
df=df.groupby(by=["15min_bin_", "event_type_", "animal_id_", "group_id", "gander", "id_cohort_","loca_ori"], as_index=False).mean()

event=["hit", "miss", "mistake", "corr_reject"]
gander=["male", "female"]
color_list=['royalblue', 'coral']
for i, p in zip(gander, color_list) :
    for j in event:        
        title="%s_%s_location_orientation_in_LH_15min_bin"%(i, j)
        t=df[(df["gander"]==i) & (df["event_type_"]==j) & (df["15min_bin_"]!="after_session") 
             & (df["15min_bin_"]!="before_session")].copy()
        t=t.rename(columns={"15min_bin_": "min_15bin"})
        
        plot=sns.pointplot(data=t, x="min_15bin", y="loca_ori_duration", hue="loca_ori",  
                           hue_order=["screenside_oriented", "screenside_blind", 
                                      "rewardside_oriented", "rewardside_blind"], palette=[p, p, p, p], 
                           linestyles=['-', '-.', ':', '--'], markers=["o", "x","D", "*"], capsize=.01, errorbar=('ci', 68))
        plot.set(xlabel="")
        plot.set(ylabel="duration(s)")
        plt.title(title)
        plt.xticks(rotation = 15, ha='right')
    #     h = plt.gca().get_lines()
    #     lg = plt.legend(handles=h, labels=['oriented', 'blind'], loc='best')
        plt.legend(bbox_to_anchor=(1.25, 1),borderaxespad=0)

        save_path=r"D:\CPT_TOT_manuscript\data processing\fig\15_loca_ori_LH"

        if os.path.isdir(save_path ) == False:
            os.mkdir(save_path)
        t.to_csv(save_path+"\%s.png"%(title), index=False)
        plt.savefig(save_path+"\%s.png"%(title), dpi=800, bbox_inches="tight")
        plt.close()

        try:
            print(title)
            data=t
            aovrm = AnovaRM(data, 'loca_ori_duration', 'id_cohort_', within=['min_15bin', 'loca_ori'])
            res = aovrm.fit()

            print(res)

            print("try pingouin method")
            aov=pg.rm_anova(dv='loca_ori_duration', within=['min_15bin', 'loca_ori'], subject='id_cohort_', data=data)

            pg.print_table(aov)

            # Optional post-hoc tests
            print("post hoc......")
            psthc=pg.pairwise_tests(dv='loca_ori_duration', within=['min_15bin', 'loca_ori'], subject='id_cohort_', data=data)
            print(psthc)

        except ValueError as e:

            print(str(e))
            print("---------------------------------------------------")
            print("using mixed effect model...")
            print("title")
            print("with interation...")
            model=smf.mixedlm("loca_ori_duration ~ C(min_15bin, Treatment('15bin_1')) + C(loca_ori) + C(min_15bin, Treatment('15bin_1')):C(loca_ori)", 
                              data=t, groups=t['id_cohort_']).fit()
            display(model.summary())

            print("                                                   ")
            print("excluding interaction...")
            model=smf.mixedlm("loca_ori_duration ~ C(min_15bin, Treatment('15bin_1')) + C(loca_ori)", 
                              data=t, groups=t['id_cohort_']).fit()    
            display(model.summary())


            print("KDE plot generating......")
            fig = plt.figure(figsize = (16, 9))
            ax=sns.kdeplot(model.resid, fill=True, linewidth=1)
            normal(model.resid.mean(), model.resid.std(), histmax=ax.get_ylim()[1])
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

path_list=[coh1_male_path, coh1_female_path, coh2_path]
coh=["coh1", "coh1", "coh2"]
delte_str=[4,12,12]
code_file_name=[coh1_male_code["file_name"].tolist(), coh1_female_code["file_name"].tolist(), coh2_tot_code["vfa_csv_2"].tolist()]
nfps=[coh1_male_code["n_fps"].tolist(), coh1_female_code["n_fps"].tolist(), coh2_tot_code["n_fps"].tolist()]
event_list=["hit", "miss", "mistake", "corr_reject"]

for i,m, h, j, d in zip(path_list,code_file_name, nfps, coh, delte_str) :
    os.chdir(i[:-7]+"\\ori_LH")
    
    save_path=i[:-7]+"\\15_ori_index_LH"
    if os.path.isdir(save_path) == False:
        os.mkdir(save_path)
        
    save_total=i[:-7]+"\\15_ori_index_LH\\total"
    if os.path.isdir(save_total) == False:
        os.mkdir(save_total)    
    
    for n, b in zip(m, h):
        df=pd.read_csv(n)
        #for c in ["hit", "miss", "mistake", "corr_reject"]:
        df=df.groupby(by=["15min_bin", "event_type", "trial_number", "blind_or_not"], as_index=False).count()
        df["cohort"]=j
        df["dt_id"]=n[:-d]
        df["animal_id"]=df["dt_id"][1][11:]
        df["id_cohort"]=df["animal_id"]+ "-" + df["cohort"]
        df["count_"]=1
        display(df)
        df_=pd.pivot_table(df, index=["15min_bin", "event_type", "trial_number", "dt_id", "animal_id", "id_cohort"], 
                          columns=["blind_or_not"], values=["count_"]).reset_index().copy()#, aggfunc="count"
        
        df_.columns=[''.join(col) for col in df_.columns.values]
        df_["count_blind"] = df_["count_blind"].fillna(0)
        df_["count_oriented"] = df_["count_oriented"].fillna(0)
        display(df_)
        
        ##### in column "blind_minus_ori_diff", if value ==  -1.0, meaning that trial is blind trial: if value=0 or 1, meaning that trial is oriented trial. 
        ##### the orientation index is calculated by oriented trial divided by trial number in each bin. 
        df_["blind_minus_ori_diff"]=df_["count_oriented"]-df_["count_blind"]         
        df_["bin_trial_n"]=0
        for x in event_list:    
            bin_condition=[(df_["15min_bin"]=="15bin_1") & (df_["event_type"]==x), (df_["15min_bin"]=="15bin_2")  & (df_["event_type"]==x),
                       (df_["15min_bin"]=="15bin_3") & (df_["event_type"]==x), (df_["15min_bin"]=="15bin_4") & (df_["event_type"]==x),
                       (df_["15min_bin"]=="15bin_5") & (df_["event_type"]==x), (df_["15min_bin"]=="15bin_6") & (df_["event_type"]==x)]
            bin_n=[]            
            for z in bin_condition:
                bin_n.append(len(df_[z]))
            df_["bin_trial_n"]=np.select(bin_condition, bin_n, default=df_["bin_trial_n"])            
        df_["onetozero"]=df_["blind_minus_ori_diff"].replace(1, 0) 
        #so in the "bind minus ori diff"col, only has value -1 or 0, -1 means blind, 0 means ori
        display(df_)
        
        df_event=df_.groupby(by=["15min_bin", "event_type", "dt_id", "animal_id", "id_cohort"], as_index=False).agg(
            only_blind_trial_sum=('onetozero', 'sum'), bin_trial_n_sum=('bin_trial_n', 'mean'))        
        df_event["ori_trial"]=df_event["bin_trial_n_sum"]+df_event["only_blind_trial_sum"]
        df_event["ori_index_event"]=df_event["ori_trial"]/df_event["bin_trial_n_sum"]
        display(df_event)
        df_event.to_csv(save_path + "\\%s"%(n), index=False)
        
        df_total=df_.groupby(by=["15min_bin", "dt_id", "animal_id", "id_cohort"], as_index=False).agg(
            only_blind_trial_sum=('onetozero', 'sum'), bin_trial_n_sum=('bin_trial_n', 'mean'))  
        df_total["ori_trial"]=df_total["bin_trial_n_sum"]+df_total["only_blind_trial_sum"]
        df_total["ori_index_total"]=df_total["ori_trial"]/df_total["bin_trial_n_sum"]        
        df_total.to_csv(save_total+"\\%s"%(n), index=False)
        display(df_total)

path_list=[coh1_male_path, coh1_female_path, coh2_path]
df_list=[]

for i in path_list:
    file_list=glob.glob(i[:-7]+"\\15_ori_index_LH\\*.csv")
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
df=df.groupby(by=["15min_bin", "event_type", "animal_id", "group_id", "gander", "id_cohort"], as_index=False).mean()

event=["hit", "miss", "mistake", "corr_reject"]
gander=["male", "female"]
color_list=['royalblue', 'coral']
for i, p in zip(gander, color_list) :
    for j in event:        
        title="%s_%s_ori_index_before_event_15min_bin"%(i, j)
        t=df[(df["gander"]==i) & (df["event_type"]==j)].copy()
        t=t.rename(columns={"15min_bin": "min_15bin"})
        
        plot=sns.pointplot(data=t, x="min_15bin", y="ori_index_event", capsize=.01, errorbar=('ci', 68))
        plot.set(xlabel="")
        plot.set(ylabel="orientation index")
        plt.title(title)
        plt.ylim(0.75, 1.05)
        plt.xticks(rotation = 15, ha='right')
    #     h = plt.gca().get_lines()
    #     lg = plt.legend(handles=h, labels=['oriented', 'blind'], loc='best')
        plt.legend(bbox_to_anchor=(1.25, 1),borderaxespad=0)

        save_path=r"D:\CPT_TOT_manuscript\data processing\fig\15_ori_index_LH"

        if os.path.isdir(save_path ) == False:
            os.mkdir(save_path)
        t.to_csv(save_path+"\%s.csv"%(title), index=False)
        plt.savefig(save_path+"\%s.png"%(title), dpi=800, bbox_inches="tight")
        plt.close()

        try:
            print(title)
            data=t[(t["min_15bin"]!="after_session") & (t["min_15bin"]!="before_session") & (t["gander"]==i)]
            aovrm = AnovaRM(data, 'ori_index_event', 'id_cohort', within=['min_15bin', 'event_type'])
            res = aovrm.fit()

            print(res)

            print("try pingouin method")
            aov=pg.rm_anova(dv='ori_index_event', within=['min_15bin', 'event_type'], subject='id_cohort', data=data)

            pg.print_table(aov)

            # Optional post-hoc tests
            print("post hoc......")
            psthc=pg.pairwise_tests(dv='ori_index_event', within=['min_15bin', 'event_type'], subject='id_cohort', data=data)
            print(psthc)

        except ValueError as e:

            print(str(e))
            print("---------------------------------------------------")
            print("using mixed effect model...")
            print("title")
            print("with interation...")
            model=smf.mixedlm("ori_index_event ~ C(min_15bin, Treatment('15bin_1')) + C(event_type) + C(min_15bin, Treatment('15bin_1')):C(event_type)", 
                              data=t, groups=t['id_cohort']).fit()
            display(model.summary())

            print("                                                   ")
            print("excluding interaction...")
            model=smf.mixedlm("ori_index_event ~ C(min_15bin, Treatment('15bin_1')) + C(event_type)", 
                              data=t, groups=t['id_cohort']).fit()    
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
            ax.set_title(title +"Q-Q Plot")
            plt.savefig(save_stat+"\%s_Q-Q_plot.png"%title)
            plt.close()

            print("KDE plot saved")

            labels = ["Statistic", "p-value"]

            norm_res = stats.shapiro(model.resid)

            for key, val in dict(zip(labels, norm_res)).items():
                print(key, val)

sns.set(style="ticks",rc={"lines.linewidth": 0.7})

path_list=[coh1_male_path, coh1_female_path, coh2_path]
df_list=[]

for i in path_list:
    file_list=glob.glob(i[:-7]+"\\15_ori_index_LH\\*.csv")
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
df=df.groupby(by=["15min_bin", "event_type", "animal_id", "group_id", "gander", "id_cohort"], as_index=False).mean()

event=["hit", "miss", "mistake", "corr_reject"]
gander=["male", "female"]
color_list=['royalblue', 'coral']
for i, p in zip(gander, color_list) :   
        title="%s_ori_index_before_event_15min_bin"%i
        t=df[(df["gander"]==i)].copy()
        t=t.rename(columns={"15min_bin": "min_15bin"})

        plot=sns.pointplot(data=t, x="min_15bin", y="ori_index_event", hue="event_type",  hue_order=["hit", "miss", "mistake", "corr_reject"], palette=[p, p, p, p], 
                            markers=["o", "x","D", "*"], capsize=.01, errorbar=('ci', 68)) #linestyles=['-', '-.', ':', '--'],
        plt.setp(plot.collections, sizes=[70])
        plot.set(xlabel="")
        plot.set(ylabel="orientation index")
        plt.title(title)
        plt.xticks(rotation = 15, ha='right')
        plot.yaxis.set_ticks(np.arange(0.80, 1.05, 0.05))
    #     h = plt.gca().get_lines()
    #     lg = plt.legend(handles=h, labels=['oriented', 'blind'], loc='best')
        plt.legend(bbox_to_anchor=(1.30, 1),borderaxespad=0)

        save_path=r"D:\CPT_TOT_manuscript\data processing\fig\15_ori_index_LH"

        if os.path.isdir(save_path ) == False:
            os.mkdir(save_path)
        t.to_csv(save_path+"\%s.csv"%(title), index=False)
        plt.savefig(save_path+"\%s.png"%(title), dpi=800, bbox_inches="tight")
        plt.close()

        try:
            print(title)
            data=t
            aovrm = AnovaRM(data, 'ori_index_event', 'id_cohort', within=['min_15bin', 'event_type'])
            res = aovrm.fit()

            print(res)

            print("try pingouin method")
            aov=pg.rm_anova(dv='ori_index_event', within=['min_15bin', 'event_type'], subject='id_cohort', data=data)

            pg.print_table(aov)

            # Optional post-hoc tests
            print("post hoc......")
            psthc=pg.pairwise_tests(dv='ori_index_event', within=['min_15bin', 'event_type'], subject='id_cohort', data=data)
            print(psthc)

        except ValueError as e:

            print(str(e))
            print("---------------------------------------------------")
            print("using mixed effect model...")
            print("with interation...")
            model=smf.mixedlm("ori_index_event ~ C(min_15bin, Treatment('15bin_1')) + C(event_type) + C(min_15bin, Treatment('15bin_1')):C(event_type)", 
                              data=t, groups=t['id_cohort']).fit()
            display(model.summary())

            print("                                                   ")
            print("excluding interaction...")
            model=smf.mixedlm("ori_index_event ~ C(min_15bin, Treatment('15bin_1')) + C(event_type)", 
                              data=t, groups=t['id_cohort']).fit()    
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
            ax.set_title(title +"Q-Q Plot")
            plt.savefig(save_stat+"\%s_Q-Q_plot.png"%title)
            plt.close()

            print("KDE plot saved")

            labels = ["Statistic", "p-value"]

            norm_res = stats.shapiro(model.resid)

            for key, val in dict(zip(labels, norm_res)).items():
                print(key, val)

sns.set(style="ticks",rc={"lines.linewidth": 0.7})

path_list=[coh1_male_path, coh1_female_path, coh2_path]
df_list=[]

for i in path_list:
    file_list=glob.glob(i[:-7]+"\\15_ori_index_LH\\total\\*.csv")
    df_list=df_list+file_list
    
#combine 15min bin df from three parts

df=pd.DataFrame()
for i in df_list:
    df_=pd.read_csv(i)
#     if len(df_)<8:
#         print(i)
    df=pd.concat([df, df_], ignore_index=True)
df['group_id']=df['animal_id'].apply(lambda x:x[0])
df['gander']=np.where(df['group_id'].isin(["W", "X", "Y","Z"]), "male", "female")
df=df.groupby(by=["15min_bin","animal_id", "group_id", "gander", "id_cohort"], as_index=False).mean()


gander=["male", "female"]
color_list=['royalblue', 'coral']
for i, p in zip(gander, color_list) :   
    title="%s_total_ori_index_in_LH_15min_bin"%i
    t=df[(df["gander"]==i)].copy()
    t=t.rename(columns={"15min_bin": "min_15bin"})

    plot=sns.pointplot(data=t, x="min_15bin", y="ori_index_total", capsize=.01, errorbar=('ci', 68), color=p) #linestyles=['-', '-.', ':', '--'],
    plt.setp(plot.collections, sizes=[70])
    plot.set(xlabel="")
    plot.set(ylabel="orientation index total")
    plt.title(title)
    plt.xticks(rotation = 15, ha='right')
    plot.yaxis.set_ticks(np.arange(0.80, 1.05, 0.05))
#     h = plt.gca().get_lines()
#     lg = plt.legend(handles=h, labels=['oriented', 'blind'], loc='best')
    plt.legend(bbox_to_anchor=(1.30, 1),borderaxespad=0)

    save_path=r"D:\CPT_TOT_manuscript\data processing\fig\15_ori_index_LH_total"

    if os.path.isdir(save_path ) == False:
        os.mkdir(save_path)
    t.to_csv(save_path+"\%s.csv"%(title), index=False)
    plt.savefig(save_path+"\%s.png"%(title), dpi=800, bbox_inches="tight")
    plt.close()
    t.to_csv(save_path+"\%s.csv"%(title), index=False)

#         try:
#             print(title)
#             data=t
#             aovrm = AnovaRM(data, 'ori_index', 'id_cohort', within=['min_15bin', 'event_type'])
#             res = aovrm.fit()

#             print(res)

#             print("try pingouin method")
#             aov=pg.rm_anova(dv='ori_index', within=['min_15bin', 'event_type'], subject='id_cohort', data=data)

#             pg.print_table(aov)

#             # Optional post-hoc tests
#             print("post hoc......")
#             psthc=pg.pairwise_tests(dv='ori_index', within=['min_15bin', 'event_type'], subject='id_cohort', data=data)
#             print(psthc)

#         except ValueError as e:

#             print(str(e))
#             print("---------------------------------------------------")
#             print("using mixed effect model...")
#             print("with interation...")
#             model=smf.mixedlm("ori_index ~ C(min_15bin, Treatment('15bin_1')) + C(event_type) + C(min_15bin, Treatment('15bin_1')):C(event_type)", 
#                               data=t, groups=t['id_cohort']).fit()
#             display(model.summary())

#             print("                                                   ")
#             print("excluding interaction...")
#             model=smf.mixedlm("ori_index ~ C(min_15bin, Treatment('15bin_1')) + C(event_type)", 
#                               data=t, groups=t['id_cohort']).fit()    
#             display(model.summary())


#             print("KDE plot generating......")
#             fig = plt.figure(figsize = (16, 9))
#             ax = sns.distplot(model.resid, hist = False, kde_kws = {"shade" : True, "lw": 1}, fit = stats.norm)
#             ax.set_title(title + "KDE Plot of Model Residuals (Blue) and Normal Distribution (Black)")
#             ax.set_xlabel("Residuals")
#             save_stat=save_path+"\stat"
#             if os.path.isdir(save_stat ) == False:
#                 os.mkdir(save_stat)

#             plt.savefig(save_stat+"\%s_KDE.png"%title)
#             plt.close()
#             print("KDE plot saved")

#             print("                                          ")

#             print("Q-Q plot generating......")
#             fig = plt.figure(figsize = (16, 9))
#             ax = fig.add_subplot(111)
#             sm.qqplot(model.resid, dist = stats.norm, line = 's', ax = ax)
#             ax.set_title(title +"Q-Q Plot")
#             plt.savefig(save_stat+"\%s_Q-Q_plot.png"%title)
#             plt.close()

#             print("KDE plot saved")

#             labels = ["Statistic", "p-value"]

#             norm_res = stats.shapiro(model.resid)

#             for key, val in dict(zip(labels, norm_res)).items():
#                 print(key, val)


