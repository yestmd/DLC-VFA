#!/usr/bin/env python
# coding: utf-8

# ## To combine the VFA orientation from 1st and 2nd cohort data

# In[1]:


import glob as glob


# In[2]:


import pandas as pd


# In[3]:


import os
import shutil


# In[4]:


import numpy as np


# In[5]:


from datetime import datetime


# ## cohort 1st vfa tot baseline results and info

# ### handle male mice

# In[6]:


# #editing abet ii raw data for syncing:change file name and add feature in file
# origin="D:\\CPT_TOT_DS_Stage_3\\CPT_DS_TOT_Visual_analysis\\male tot abet raw data" 

# target=origin + "\\edit_raw_file"

# files=[fileName for fileName in os.listdir(origin) if fileName.endswith(".csv")]

# os.chdir(origin)
# for name in files:
#     if os.path.isdir(target) == False:
#         os.mkdir(target)

#     df=pd.read_csv(name)
#     #df.loc[:, "Schedule run date short"]=pd.to_datetime(df["Date/Time"]).dt.date
#     df["Schedule run date short"]=pd.to_datetime(df["Date/Time"], format="%m/%d/%Y %I:%M:%S %p")
#     df["Schedule run date short"]=df["Schedule run date short"].apply(lambda x: datetime.strftime(x, '%m-%d-%Y'))
#     df["file_name"]=df["Schedule run date short"] + "-" + df["Group ID"] + "-" + df["Animal ID"]
    
#     df=df.drop(['Machine Name', 'Date/Time', 'Version', 'Version Name', 'Application_Version', 
#                 'Experiment', 'Max_Number_Trials', 'Schedule_Description', 'Environment', 
#                 'Analysis Name', 'Schedule Run ID', 'Start ITI - Start'], axis=1)
#     #drop "Start ITI - Start" column, since it impact the orientation determination in next step
#     df=df.melt(id_vars=['Database','Schedule Name', 'Schedule_Start_Time', 'Schedule run date short', 'Animal ID', 'Group ID', 'file_name', 'Max_Schedule_Time'], var_name="event", value_name="sec_cal")
        
#     n=df["file_name"][1]
#     df.to_csv(target + "\\%s.csv"%n)
#     print(df["file_name"][1] + ":Files's names are updated")


# In[7]:


#load vfa results
#coh1_male=glob.glob(r"D:\CPT_TOT_DS_Stage_3\CPT_DS_TOT_Visual_analysis\vfa male tot_based on 1000000 model_fixed\dlc_results_namechanged\results\*.csv")

#load file name code and started frame info
coh1_male_code=pd.read_csv(r"D:\CPT_TOT_DS_Stage_3\CPT_DS_TOT_Visual_analysis\file_frame_info_male.csv", index_col=False)


# In[8]:


# #update code file
# # coh1_male_code=coh1_male_code.rename(columns={'fps':'n_fps', 'fm':'frame'})
# coh1_male_code.loc[:,'fm_started']=(coh1_male_code['min']*60 + coh1_male_code['sec'])*coh1_male_code['n_fps'] + coh1_male_code['frame']
# # coh1_male_code["abet_raw"]=coh1_male_code["name_2"] + '.csv'
# # coh1_male_code["abet_raw"]=coh1_male_code["abet_raw"].apply(lambda x: x.replace('-TOT', ''))
# #coh1_male_code['file_name']=coh1_male_code["name_2"]+".csv"
# coh1_male_code.to_csv(r"D:\\CPT_TOT_DS_Stage_3\\CPT_DS_TOT_Visual_analysis\\file_frame_info_male.csv", index=False)


# In[9]:


coh1_male_code.head(2)


# In[10]:


#add feature for vfa result
origin="D:\\CPT_TOT_DS_Stage_3\\CPT_DS_TOT_Visual_analysis\\vfa male tot_based on 1000000 model_fixed\\dlc_results_namechanged\\results"

target=origin+ "\\temp"

# change frames_cal column into each vfa results file 

#add temp folder into target path
if os.path.isdir(target) == False:
    os.mkdir(target)
    
#files=[fileName for fileName in os.listdir(origin) if fileName.endswith(".csv")]

os.chdir(origin)

#add sec_cal and direction feature to each vfa result file    
for i, j, h, n in zip(coh1_male_code['new_csv'], coh1_male_code["fm_started"],
                      coh1_male_code["n_fps"], coh1_male_code['abet_raw']):
    df=pd.read_csv(i)
    df['frames_cal'] = df['index']-j+1 
    df['sec_cal']=(df["frames_cal"]-1)*1/h+0.001
    df['n_fps']=h
    df["total_visual"]=df["lateralRightright"] + df["lateralLeftright"] + df["frontalright"]
    df['blind_or_not']=np.where(df['blindright']>=0.9, 'blind', 'oriented')
    
    #save editted vfa result file into temp with recovered name
    df.to_csv(target+ "\\%s"%(n), index=False)


# In[11]:


# check file name match
raw_list=glob.glob("D:\\CPT_TOT_DS_Stage_3\\CPT_DS_TOT_Visual_analysis\\male tot abet raw data\\edit_raw_file\\\*.csv")
raw_list.sort()

target="D:\\CPT_TOT_DS_Stage_3\\CPT_DS_TOT_Visual_analysis\\vfa male tot_based on 1000000 model_fixed\\dlc_results_namechanged\\results\\temp"
os.chdir(target)
vfa_list=glob.glob(r"*.csv")
vfa_list.sort()

d={'raw_list':raw_list, 'vfa_list': vfa_list}
file_list=pd.DataFrame(data=d)
file_list['match']=file_list.apply(lambda x : x.vfa_list[:-4] in x.raw_list, axis=1)
file_list['match'].unique()


# In[12]:


#sync abet raw data to vfa result
# merge cpt raw data and vfa data to sync event
#save processed file into the pathway below:

pathfile="D:\\CPT_TOT_DS_Stage_3\\CPT_DS_TOT_Visual_analysis\\vfa male tot_based on 1000000 model_fixed\\dlc_results_namechanged\\results\\temp_2"
if os.path.isdir(pathfile) == False:
    os.mkdir(pathfile)

target="D:\\CPT_TOT_DS_Stage_3\\CPT_DS_TOT_Visual_analysis\\vfa male tot_based on 1000000 model_fixed\\dlc_results_namechanged\\results\\temp"
os.chdir(target)

for i, j in zip(raw_list, vfa_list):
    df=pd.read_csv(i, index_col=None, header=0)
    evnt=df[["event", "sec_cal"]].copy()
    evnt.dropna(subset=['sec_cal'], inplace=True)
    evnt["count"]=1
    #load vfa data and reform the format of the dataframe    
    df=pd.read_csv(j)
    
    #seperating handling Center Screen Touch - Time event, since it comes with Hit, 
    #and makes the orientated hit tobe undeterminded
   
    #merge cpt raw data, not includes center screen touch event and Stim Onset - End ITI, and vfa data
    combine_1=pd.concat([evnt[(evnt["event"]!="Center Screen Touch - Time") & (evnt["event"]!="Stim Onset - End ITI")], df], 
                        ignore_index=False).sort_values(by='sec_cal').reset_index()

    #check each row, add ori_before and ori_after
    combine_1['ori_before']=""
    combine_1["ori_after"]=""

    for h in range(len(combine_1)):
        if combine_1["count"][h]==1:
            combine_1.loc[:, "ori_before"][h]=combine_1['blind_or_not'][h-1]
            combine_1.loc[:, "ori_after"][h]=combine_1['blind_or_not'][h+1]
                    
    combine_1['ori'] = np.where(combine_1["ori_before"] == combine_1["ori_after"], combine_1['ori_after'], "undetermined")
    
    
    #merge cpt raw data only includes center screen touch event and vfa data
    combine_2=pd.concat([evnt[evnt["event"]=="Center Screen Touch - Time"], df], ignore_index=False).sort_values(by='sec_cal').reset_index()

    #check each row, add ori_before and ori_after
    combine_2['ori_before']=""
    combine_2["ori_after"]=""

    for h in range(len(combine_2)):
        if combine_2["count"][h]==1:
            combine_2.loc[:, "ori_before"][h]=combine_2['blind_or_not'][h-1]
            combine_2.loc[:, "ori_after"][h]=combine_2['blind_or_not'][h+1]
                    
    combine_2['ori'] = np.where(combine_2["ori_before"] == combine_2["ori_after"], combine_2['ori_after'], "undetermined")

    
    #merge cpt raw data only Stim Onset - End ITI and vfa data
    combine_3=pd.concat([evnt[evnt["event"]=="Stim Onset - End ITI"], df], ignore_index=False).sort_values(by='sec_cal').reset_index()

    #check each row, add ori_before and ori_after
    combine_3['ori_before']=""
    combine_3["ori_after"]=""

    for h in range(len(combine_2)):
        if combine_3["count"][h]==1:
            combine_3.loc[:, "ori_before"][h]=combine_3['blind_or_not'][h-1]
            combine_3.loc[:, "ori_after"][h]=combine_3['blind_or_not'][h+1]
                    
    combine_3['ori'] = np.where(combine_3["ori_before"] == combine_3["ori_after"], combine_3['ori_after'], "undetermined")   
    
    #combine the two part results
    combine=pd.concat([combine_1, combine_2[combine_2["event"]=="Center Screen Touch - Time"], 
                       combine_3[combine_3["event"]=="Stim Onset - End ITI"]], ignore_index=False).sort_values(by='sec_cal')
    
    
    combine=combine.drop(['level_0'], axis=1)
    
    #add 45min bin and 15min bin feature

    conditions_15=[(combine["sec_cal"]<0.000),
            (combine["sec_cal"]>=0.000) & (combine["sec_cal"]<=900.000),
            (combine["sec_cal"]>900.000) & (combine["sec_cal"]<=1800.000), 
            (combine["sec_cal"]>1800.000) & (combine["sec_cal"]<=2700.000),
            (combine["sec_cal"]>2700.000) & (combine["sec_cal"]<=3600.000),
            (combine["sec_cal"]>3600.000) & (combine["sec_cal"]<=4500.000),
            (combine["sec_cal"]>4500.000) & (combine["sec_cal"]<=5400.000),
                  (combine["sec_cal"]>5400.000)]

    values_15=['before_session','15bin_1', '15bin_2', '15bin_3', '15bin_4', '15bin_5', '15bin_6', 'after_session'] 

    conditions_45=[(combine["sec_cal"]<0.000),
                   (combine["sec_cal"]>=0.000) & (combine["sec_cal"]<=2700.000),
            (combine["sec_cal"]>2700.000) & (combine["sec_cal"]<=5400.000), (combine["sec_cal"]>5400.000)]

    values_45=['before_session',"first_45min", "last_45min", "after_session"] 
    
       
    combine['45min_bin']=np.select(conditions_45, values_45)
    
    combine['15min_bin']=np.select(conditions_15, values_15)
    
    combine=combine.round(3)
    combine['cohort']='coh1'
    
    combine.to_csv(pathfile + "\\%s"%j)


# ### handle female mice

# In[13]:


# #editing abet ii raw data for syncing:change file name and add feature in file
# origin="D:\\CPT_TOT_DS_Stage_3\\CPT_DS_TOT_Visual_analysis\\female tot abet raw data\\" 

# target=origin + "edit_raw_file\\"

# files=[fileName for fileName in os.listdir(origin) if fileName.endswith(".csv")]

# os.chdir(origin)
# for name in files:
#     if os.path.isdir(target) == False:
#         os.mkdir(target)

#     df=pd.read_csv(name)
#     #df.loc[:, "Schedule run date short"]=pd.to_datetime(df["Date/Time"]).dt.date
#     df["Schedule run date short"]=pd.to_datetime(df["Date/Time"], format="%m/%d/%Y %I:%M:%S %p")
#     df["Schedule run date short"]=df["Schedule run date short"].apply(lambda x: datetime.strftime(x, '%m-%d-%Y'))
#     df["file_name"]=df["Schedule run date short"] + "-" + df["Group ID"] + "-" + df["Animal ID"]
    
#     df=df.drop(['Machine Name', 'Date/Time', 'Version', 'Version Name', 'Application_Version', 
#                 'Experiment', 'Max_Number_Trials', 'Schedule_Description', 'Environment', 
#                 'Analysis Name', 'Schedule Run ID', 'Start ITI - Start'], axis=1)
    
#     df=df.melt(id_vars=['Database','Schedule Name', 'Schedule_Start_Time', 'Schedule run date short', 'Animal ID', 'Group ID', 'file_name', 'Max_Schedule_Time'], var_name="event", value_name="sec_cal")
    
#     n=df["file_name"][1]
#     df.to_csv(target + "%s.csv"%n)
#     print(df["file_name"][1] + ":Files's names are updated")


# In[14]:


#load female file frame info
coh1_female_code=pd.read_csv(r"D:\CPT_TOT_DS_Stage_3\CPT_DS_TOT_Visual_analysis\file_frame_info_female.csv", index_col=False)
coh1_female_code.head(2)


# In[15]:


# #update code file
# # coh1_female_code=coh1_female_code.rename(columns={'fps':'n_fps', 'fm':'frame'})
# coh1_female_code.loc[:,'fm_started']=(coh1_female_code['min']*60 + coh1_female_code['sec'])*coh1_female_code['n_fps'] + coh1_female_code['frame']
# # # coh1_female_code["name_2"]=coh1_female_code["name_2"].str.replace('-tot', '')
# # # coh1_female_code['file_name']=coh1_female_code["name_2"]+".csv"
# # coh1_female_code.head(2)
# coh1_female_code.to_csv(r"D:\CPT_TOT_DS_Stage_3\CPT_DS_TOT_Visual_analysis\file_frame_info_female.csv", index=False)


# In[16]:


#change vfa result's name and add feature for vfa result
origin="D:\\CPT_TOT_DS_Stage_3\\CPT_DS_TOT_Visual_analysis\\vfa female tot_based on 1030000 model\\results"

# change frames_cal column into each vfa results file 

#add temp folder into target path
pathfile=origin+"\\temp"
if os.path.isdir(pathfile) == False:
    os.mkdir(pathfile)
    
#files=[fileName for fileName in os.listdir(origin) if fileName.endswith(".csv")]

os.chdir(origin)

#add sec_cal and direction feature to each vfa result file    
for i, j, h, n in zip(coh1_female_code['new_csv'], coh1_female_code["fm_started"], 
                      coh1_female_code["n_fps"], coh1_female_code['file_name']):
    df=pd.read_csv(i)
    df['frames_cal'] = df['index']-j+1 
    df['sec_cal']=(df["frames_cal"]-1)*1/h+0.001
    df['n_fps']=h
    df["total_visual"]=df["lateralRightright"] + df["lateralLeftright"] + df["frontalright"]
    df['blind_or_not']=np.where(df['blindright']>=0.9, 'blind', 'oriented')
    
    df.to_csv(pathfile+ "\\%s"%(n), index=False)


# In[17]:


# check file name match
raw_list=glob.glob("D:\\CPT_TOT_DS_Stage_3\\CPT_DS_TOT_Visual_analysis\\female tot abet raw data\\edit_raw_file\\*.csv")
raw_list.sort()

os.chdir(origin+"\\temp\\")
vfa_list=glob.glob(r"*.csv")
vfa_list.sort()

d={'raw_list':raw_list, 'vfa_list': vfa_list}
file_list=pd.DataFrame(data=d)
file_list['match']=file_list.apply(lambda x : x.vfa_list[:-12] in x.raw_list, axis=1)
file_list['match'].unique()


# In[18]:


len(vfa_list)


# In[19]:


#sync abet raw data to vfa result
# merge cpt raw data and vfa data to sync event
from scipy.stats import norm
import math

pathfile=origin+"\\temp_2\\"
if os.path.isdir(pathfile) == False:
    os.mkdir(pathfile)

for i, j in zip(raw_list, vfa_list):
    df=pd.read_csv(i, index_col=None, header=0)
    evnt=df[["event", "sec_cal"]].copy()
    evnt.dropna(subset=['sec_cal'], inplace=True)
    evnt["count"]=1
    #load vfa data and reform the format of the dataframe    
    df=pd.read_csv(j)
    #seperating handling Center Screen Touch - Time event, since it comes with Hit, 
    #and makes the orientated hit tobe undeterminded
   
    #merge cpt raw data, not includes center screen touch event and Stim Onset - End ITI, and vfa data
    combine_1=pd.concat([evnt[(evnt["event"]!="Center Screen Touch - Time") & (evnt["event"]!="Stim Onset - End ITI")], df], 
                        ignore_index=False).sort_values(by='sec_cal').reset_index()

    #check each row, add ori_before and ori_after
    combine_1['ori_before']=""
    combine_1["ori_after"]=""

    for h in range(len(combine_1)):
        if combine_1["count"][h]==1:
            combine_1.loc[:, "ori_before"][h]=combine_1['blind_or_not'][h-1]
            combine_1.loc[:, "ori_after"][h]=combine_1['blind_or_not'][h+1]
                    
    combine_1['ori'] = np.where(combine_1["ori_before"] == combine_1["ori_after"], combine_1['ori_after'], "undetermined")
    
    
    #merge cpt raw data only includes center screen touch event and vfa data
    combine_2=pd.concat([evnt[evnt["event"]=="Center Screen Touch - Time"], df], ignore_index=False).sort_values(by='sec_cal').reset_index()

    #check each row, add ori_before and ori_after
    combine_2['ori_before']=""
    combine_2["ori_after"]=""

    for h in range(len(combine_2)):
        if combine_2["count"][h]==1:
            combine_2.loc[:, "ori_before"][h]=combine_2['blind_or_not'][h-1]
            combine_2.loc[:, "ori_after"][h]=combine_2['blind_or_not'][h+1]
                    
    combine_2['ori'] = np.where(combine_2["ori_before"] == combine_2["ori_after"], combine_2['ori_after'], "undetermined")

    
    #merge cpt raw data only Stim Onset - End ITI and vfa data
    combine_3=pd.concat([evnt[evnt["event"]=="Stim Onset - End ITI"], df], ignore_index=False).sort_values(by='sec_cal').reset_index()

    #check each row, add ori_before and ori_after
    combine_3['ori_before']=""
    combine_3["ori_after"]=""

    for h in range(len(combine_2)):
        if combine_3["count"][h]==1:
            combine_3.loc[:, "ori_before"][h]=combine_3['blind_or_not'][h-1]
            combine_3.loc[:, "ori_after"][h]=combine_3['blind_or_not'][h+1]
                    
    combine_3['ori'] = np.where(combine_3["ori_before"] == combine_3["ori_after"], combine_3['ori_after'], "undetermined")   
    
    #combine the two part results
    combine=pd.concat([combine_1, combine_2[combine_2["event"]=="Center Screen Touch - Time"], 
                       combine_3[combine_3["event"]=="Stim Onset - End ITI"]], ignore_index=False).sort_values(by='sec_cal')
    
    #combine[(combine['ori_before'] != combine['ori_after']), ['ori']] = "undetermined"
    combine=combine.drop(['level_0'], axis=1)
    
    #add 45min bin and 15min bin feature

    conditions_15=[(combine["sec_cal"]<0.000),
            (combine["sec_cal"]>=0.000) & (combine["sec_cal"]<=900.000),
            (combine["sec_cal"]>900.000) & (combine["sec_cal"]<=1800.000), 
            (combine["sec_cal"]>1800.000) & (combine["sec_cal"]<=2700.000),
            (combine["sec_cal"]>2700.000) & (combine["sec_cal"]<=3600.000),
            (combine["sec_cal"]>3600.000) & (combine["sec_cal"]<=4500.000),
            (combine["sec_cal"]>4500.000) & (combine["sec_cal"]<=5400.000),
                  (combine["sec_cal"]>5400.000)]

    values_15=['before_session','15bin_1', '15bin_2', '15bin_3', '15bin_4', '15bin_5', '15bin_6', 'after_session'] 

    conditions_45=[(combine["sec_cal"]<0.000),
                   (combine["sec_cal"]>=0.000) & (combine["sec_cal"]<=2700.000),
            (combine["sec_cal"]>2700.000) & (combine["sec_cal"]<=5400.000), (combine["sec_cal"]>5400.000)]

    values_45=['before_session',"first_45min", "last_45min", "after_session"] 
    
       
    combine['45min_bin']=np.select(conditions_45, values_45)
    
    combine['15min_bin']=np.select(conditions_15, values_15)
    
    combine=combine.round(3)
    combine['cohort']='coh1'
    
    combine.to_csv(origin+"\\temp_2\\%s"%j)


# In[ ]:





# # get 2nd cohort vfa tot baseline results and info

# ### recover file name in result folder

# In[20]:


# #editing abet ii raw data for syncing:change file name and add feature in file
# origin="D:\\CPT_TOT_DS_Stage_3_2nd\\VFA_baseline\\ABET II raw data" 

# target=origin + "\\edit_raw_file"

# files=[fileName for fileName in os.listdir(origin) if fileName.endswith(".csv")]

# os.chdir(origin)
# for name in files:
#     if os.path.isdir(target) == False:
#         os.mkdir(target)

#     df=pd.read_csv(name)
#     #df.loc[:, "Schedule run date short"]=pd.to_datetime(df["Date/Time"]).dt.date
#     df["Schedule run date short"]=pd.to_datetime(df["Date/Time"], format="%m/%d/%Y %I:%M:%S %p")
#     df["Schedule run date short"]=df["Schedule run date short"].apply(lambda x: datetime.strftime(x, '%m-%d-%Y'))
#     df["file_name"]=df["Schedule run date short"] + "-" + df["Group ID"] + "-" + df["Animal ID"]
    
#     df=df.drop(['Machine Name', 'Date/Time', 'Version', 'Version Name', 'Application_Version', 
#                 'Experiment', 'Max_Number_Trials', 'Schedule_Description', 'Environment', 
#                 'Analysis Name', 'Schedule Run ID', 'Start ITI - Start'], axis=1)
    
#     df=df.melt(id_vars=['Database','Schedule Name', 'Schedule_Start_Time', 'Schedule run date short', 'Animal ID', 'Group ID', 'file_name', 'Max_Schedule_Time'], var_name="event", value_name="sec_cal")
        
#     n=df["file_name"][1]
#     df.to_csv(target + "\\%s.csv"%n)
#     print(df["file_name"][1] + ":Files's names are updated")


# In[21]:


# #for update name_code
# name_code=pd.read_csv(r"D:\CPT_TOT_DS_Stage_3_2nd\VFA_baseline\name_code.csv")
# started_frame=pd.read_excel(r"D:\CPT_TOT_DS_Stage_3_2nd\VFA_baseline\started frame_.xlsx",sheet_name='Sheet2', index_col=None, header=0)

# coh2_tot_code=started_frame.merge(name_code, on="file_name", how="outer")
# coh2_tot_code=coh2_tot_code.rename(columns={'vfa_id':'vfa_csv_1', 'csv_file_name':'vfa_csv_2'})
coh2_tot_code=pd.read_csv(r"D:\CPT_TOT_DS_Stage_3_2nd\VFA_baseline\coh2_tot_code.csv", index_col=False)
# coh2_tot_code.loc[:,'fm_started']=(coh2_tot_code['min']*60 + coh2_tot_code['sec'])*coh2_tot_code['n_fps'] + coh2_tot_code['frame']
# # coh2_tot_code["abet_raw"]=coh2_tot_code["vfa_csv_2"].apply(lambda x: x.replace('-TOT-1st', ''))
# # coh2_tot_code["abet_raw"]=coh2_tot_code["abet_raw"].apply(lambda x: x.replace('-TOT-2nd', ''))
# # #coh2_tot_code['file_name']=coh2_tot_code["name_2"]+".csv"

# coh2_tot_code.to_csv(r"D:\CPT_TOT_DS_Stage_3_2nd\VFA_baseline\coh2_tot_code.csv")


# In[22]:


#add feature for vfa result
origin="D:\\CPT_TOT_DS_Stage_3_2nd\\VFA_baseline\\DLC_output\\results"

save_path=origin+"\\temp"
# change frames_cal column into each vfa results file 

#add temp folder into target path
if os.path.isdir(save_path) == False:
    os.mkdir(save_path)
    
#files=[fileName for fileName in os.listdir(origin) if fileName.endswith(".csv")]

os.chdir(origin)

#add sec_cal and direction feature to each vfa result file    
for i, j, h, n in zip(coh2_tot_code['vfa_csv_1'], coh2_tot_code["fm_started"], 
                      coh2_tot_code["n_fps"], coh2_tot_code['vfa_csv_2']):
    df=pd.read_csv(i)
    df['frames_cal'] = df['index']-j+1 
    df['sec_cal']=(df["frames_cal"]-1)*1/h+0.001
    df['n_fps']=h
    df["total_visual"]=df["lateralRightright"] + df["lateralLeftright"] + df["frontalright"]
    df['blind_or_not']=np.where(df['blindright']>=0.9, 'blind', 'oriented')
    
    df.to_csv(save_path+ "\\%s"%(n), index=False)


# In[23]:


# check file name match
raw_list=glob.glob("D:\\CPT_TOT_DS_Stage_3_2nd\\VFA_baseline\\ABET II raw data\\edit_raw_file\\*.csv")
raw_list.sort()

save_path="D:\\CPT_TOT_DS_Stage_3_2nd\\VFA_baseline\\DLC_output\\results\\temp"
os.chdir(save_path)
vfa_list=glob.glob(r"*.csv")
# 04-07-2022-Z-R-16-TOT-2nd.csv was not recorded in ABET 
#since ABET program error, 
#remove it from list
#vfa_list.remove("04-07-2022-Z-R-16-TOT-2nd.csv")
#vfa_list.sort()

d={'raw_list':raw_list, 'vfa_list': vfa_list}
file_list=pd.DataFrame(data=d)
file_list['match']=file_list.apply(lambda x : x.vfa_list[:-12] in x.raw_list, axis=1)
file_list['match'].unique()


# In[24]:


file_list


# In[25]:


#sync abet raw data to vfa result
# merge cpt raw data and vfa data to sync event
from scipy.stats import norm
import math

#save the processed file into the path below
pathfile="D:\\CPT_TOT_DS_Stage_3_2nd\\VFA_baseline\\DLC_output\\results\\temp_2"
if os.path.isdir(pathfile) == False:
    os.mkdir(pathfile)

save_path="D:\\CPT_TOT_DS_Stage_3_2nd\\VFA_baseline\\DLC_output\\results\\temp"
os.chdir(save_path)

for i, j in zip(raw_list, vfa_list):
    df=pd.read_csv(i, index_col=None, header=0)
    evnt=df[["event", "sec_cal"]].copy()
    evnt.dropna(subset=['sec_cal'], inplace=True)
    evnt["count"]=1
    #load vfa data and reform the format of the dataframe    
    df=pd.read_csv(j)
    #seperating handling Center Screen Touch - Time event, since it comes with Hit, 
    #and makes the orientated hit tobe undeterminded
   
    #merge cpt raw data, not includes center screen touch event and Stim Onset - End ITI, and vfa data
    combine_1=pd.concat([evnt[(evnt["event"]!="Center Screen Touch - Time") & (evnt["event"]!="Stim Onset - End ITI")], df], 
                        ignore_index=False).sort_values(by='sec_cal').reset_index()

    #check each row, add ori_before and ori_after
    combine_1['ori_before']=""
    combine_1["ori_after"]=""

    for h in range(len(combine_1)):
        if combine_1["count"][h]==1:
            combine_1.loc[:, "ori_before"][h]=combine_1['blind_or_not'][h-1]
            combine_1.loc[:, "ori_after"][h]=combine_1['blind_or_not'][h+1]
                    
    combine_1['ori'] = np.where(combine_1["ori_before"] == combine_1["ori_after"], combine_1['ori_after'], "undetermined")
    
    
    #merge cpt raw data only includes center screen touch event and vfa data
    combine_2=pd.concat([evnt[evnt["event"]=="Center Screen Touch - Time"], df], ignore_index=False).sort_values(by='sec_cal').reset_index()

    #check each row, add ori_before and ori_after
    combine_2['ori_before']=""
    combine_2["ori_after"]=""

    for h in range(len(combine_2)):
        if combine_2["count"][h]==1:
            combine_2.loc[:, "ori_before"][h]=combine_2['blind_or_not'][h-1]
            combine_2.loc[:, "ori_after"][h]=combine_2['blind_or_not'][h+1]
                    
    combine_2['ori'] = np.where(combine_2["ori_before"] == combine_2["ori_after"], combine_2['ori_after'], "undetermined")

    
    #merge cpt raw data only Stim Onset - End ITI and vfa data
    combine_3=pd.concat([evnt[evnt["event"]=="Stim Onset - End ITI"], df], ignore_index=False).sort_values(by='sec_cal').reset_index()

    #check each row, add ori_before and ori_after
    combine_3['ori_before']=""
    combine_3["ori_after"]=""

    for h in range(len(combine_2)):
        if combine_3["count"][h]==1:
            combine_3.loc[:, "ori_before"][h]=combine_3['blind_or_not'][h-1]
            combine_3.loc[:, "ori_after"][h]=combine_3['blind_or_not'][h+1]
                    
    combine_3['ori'] = np.where(combine_3["ori_before"] == combine_3["ori_after"], combine_3['ori_after'], "undetermined")   
    
    #combine the two part results
    combine=pd.concat([combine_1, combine_2[combine_2["event"]=="Center Screen Touch - Time"], 
                       combine_3[combine_3["event"]=="Stim Onset - End ITI"]], ignore_index=False).sort_values(by='sec_cal')
    
    #combine[(combine['ori_before'] != combine['ori_after']), ['ori']] = "undetermined"
    combine=combine.drop(['level_0'], axis=1)
    
    #add 45min bin and 15min bin feature

    conditions_15=[(combine["sec_cal"]<0.000),
            (combine["sec_cal"]>=0.000) & (combine["sec_cal"]<=900.000),
            (combine["sec_cal"]>900.000) & (combine["sec_cal"]<=1800.000), 
            (combine["sec_cal"]>1800.000) & (combine["sec_cal"]<=2700.000),
            (combine["sec_cal"]>2700.000) & (combine["sec_cal"]<=3600.000),
            (combine["sec_cal"]>3600.000) & (combine["sec_cal"]<=4500.000),
            (combine["sec_cal"]>4500.000) & (combine["sec_cal"]<=5400.000),
                  (combine["sec_cal"]>5400.000)]

    values_15=['before_session','15bin_1', '15bin_2', '15bin_3', '15bin_4', '15bin_5', '15bin_6', 'after_session'] 

    conditions_45=[(combine["sec_cal"]<0.000),
                   (combine["sec_cal"]>=0.000) & (combine["sec_cal"]<=2700.000),
            (combine["sec_cal"]>2700.000) & (combine["sec_cal"]<=5400.000), (combine["sec_cal"]>5400.000)]

    values_45=['before_session',"first_45min", "last_45min", "after_session"] 
    
       
    combine['45min_bin']=np.select(conditions_45, values_45)
    
    combine['15min_bin']=np.select(conditions_15, values_15)
    
    combine=combine.round(3)
    combine['cohort']='coh2'
    
    combine.to_csv(pathfile + "\\%s"%j)


# In[ ]:




