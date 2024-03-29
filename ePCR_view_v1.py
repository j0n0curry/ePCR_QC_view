import streamlit as st
import os
from io import BytesIO
from io import StringIO
from io import TextIOWrapper
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import math

version = 'v3 - alpha'

    
#Set up main page of application / Header info / data collection / file selection / remove files / Reset

#Main landing page greeting / info




st.set_page_config(layout="wide")

st.title('ePCR analysis tool ' +str(version) + ' for Araya 1 RFL - 100 100 100 calibration')

st.write('Developed by: jonathan.curry@lgcgroup.com - alpha version - for bugs, which are likely - please email files used for debugging purposes')

st.subheader("Upload Araya csv file directly for processing - either drag and drop or click to select files on your local machine")


###function parser - parse araya files - instatiated with ArayaManager - functions can be accessed with . notation
class WellDataManager:
    """Parent class for importing data from any instrument"""
    def __init__(self,
                 files=None,
                 dfs=None,
                 run_name_split_loc=1,
                 date_time_split_loc=3,
                 group_name="",
                 replicate_name=""):

        # project attributes
        self.file_names = files
        self.df_list = dfs
        self.group_name = group_name
        self.replicate_name = replicate_name
        self.run_name = ""
        self.date_time = ""
        self.split_char_loc = run_name_split_loc
        self.run_df = pd.DataFrame()
        self.group_df = pd.DataFrame()

        # csv read attributes
        self.tables = 1
        self.index_column = [0]
        self.header_row = [1]
        self.row_count = [8]
        self.column_names = ['Row_ID', 'Col_ID', 'Value']

    def concatenate_dataframes(self, file_names):
        if file_names is not None:
            self.file_names = file_names
            self.build_dataframes(self.file_names)
            self.group_df = pd.concat([self.group_df, self.run_df], ignore_index=True)
        else:
            for each_file in self.file_names:
                self.file_names = file_names
                self.build_dataframes(each_file)
                self.group_df = pd.concat([self.group_df, self.run_df], ignore_index=True)
        return self.group_df

    def build_dataframes(self, each_file):
        self.read_file(each_file)
        self.coerce_numeric_values()
        self.run_df['Group_ID'] = self.group_name
        #self.run_df['File_root'] = each_file
        self.run_df['Run_ID'] = self.run_name
        self.run_df['date_time'] = self.date_time
        # print(self.run_df)

    def coerce_numeric_values(self):
        # may be used used to force values to numeric.  Not applicable to all instruments
        # may be used used to force values to numeric.  Not applicable to all instruments
        pass

    def read_file(self, file_name):
        """Reads Initial Data from CSV file"""
        df = pd.read_csv(file_name, header=self.header_row, nrows=self.row_count, index_col=self.index_column, dtype = 'Int32', thousands=',')
        df = df.stack()
        self.run_df = df.reset_index()
        self.run_df.columns = self.column_names

    def get_run_name(self, file_name):
        """Splits string to get run name from file name."""
        self.run_name = file_name[:self.split_char_loc]
    def get_date_time(self, file_name):
        self.date_time = file_name[-(self.split_char_loc+1):-3]
        #print(self.run_name)





#Select main function instantiate ArayaManager before for loop otherwise reads and closes files



class ArayaManager(WellDataManager):
    """Class that handles Well Data Data"""
    def __init__(self,
                 files=None,
                 dfs=None,
                 run_name_split_loc=6,
                 date_time_split_loc=1,
                 group_name="",
                 replicate_name="",
                 dyes=None,
                 separate_column=True):
        super().__init__(files,
                         dfs,
                         run_name_split_loc,
                         group_name,
                         replicate_name)

        if dyes is None:
            dyes = ['FAM', 'VIC', 'ROX']

        # Ayara-specific items
        self.separate_column_per_dye = separate_column
        self.channel_df = pd.DataFrame()
        self.dyes = dyes
        self.channels = ['CH1', 'CH2', 'CH3']
        self.length = 14

        # csv read attributes
        self.tables = 3
        self.index_column = ["<>", "<>", "<>"]
        self.header_row = [5, 23, 41]
        self.row_count = [16, 16, 16]

        if self.separate_column_per_dye:
            # TODO: generalize for other dye names
            self.column_names = ['Row_ID', 'Col_ID', 'FAM_RFU', 'VIC_RFU', 'ROX_RFU']
        else:
            self.column_names = ['Row_ID', 'Col_ID', 'RFU', 'Channel', 'Dye']

    def read_each_channel(self, df, ch):
        """Reads Individual Channel Data from CSV file"""

        # Need to shift to get rid of annoying '<>'.  Otherwise won't parse correctly.

        if '<>' in df:
            df = df.shift(periods=1, axis='columns')
            df.drop('<>', axis=1, inplace=True)

        # Stack df for various dyes and add additional columns
        df = df.stack()
        self.channel_df = df.reset_index()

        # For separate columns for each dye, rename RFU columns.  pd.concat() method does the rest!
        if self.separate_column_per_dye:
            self.channel_df.columns = self.column_names[0:3]
            self.channel_df.rename(columns={'FAM_RFU': f'{self.dyes[ch]}_RFU'},
                                   inplace=True)

        # case to stack all dyes into common RFU and Dye channels.
        else:
            self.channel_df['Channel'] = self.channels[ch]
            self.channel_df['Dye'] = self.dyes[ch]
            self.channel_df.columns = self.column_names


    def read_file(self, file_name):
        """Reads Each Channel Data from CSV file"""

        # loops through the 3 channel tables in the csv output files.
        self.run_df = pd.DataFrame()

        df = pd.read_csv(file_name,
                         header=self.header_row[0],
                         na_values="<>", thousands=',')
        dyes = ['FAM', 'VIC', 'ROX']
        dfs = {dye: pd.DataFrame() for dye in dyes}
        ranges = [(0, 16), (18, 34), (36, 52)]
        for i, dye in zip(ranges, dyes):
            dfs[dye] = df[i[0]:i[1]]

        for ch in range(self.tables):
            dye = dyes[ch]
            self.read_each_channel(dfs[dye], ch)

            # case to have separate columns for each dye
            if self.separate_column_per_dye:
                self.channel_df = self.channel_df[self.channel_df.columns.difference(self.run_df.columns)]
                self.run_df = pd.concat([self.run_df, self.channel_df], axis=1)

            # case to stack all dyes into common RFU and Dye channels.
            else:
                self.run_df = pd.concat([self.run_df, self.channel_df], ignore_index=True)

        # Force columns to correct order.  Fixed bug with concat of separate dye columns.
        self.run_df = self.run_df[self.column_names]

    def get_run_name(self, file_name):
        """Splits string to get run name from file name."""
        self.run_name = file_name[-(self.split_char_loc+4):-4]
        """Splits string of file name to collect date time"""
    def get_date_time(self, file_name):
        self.date_time = file_name[-(self.split_char_loc+50):-38]
        #add in statement to deal with appended num in duplicated upload files.
        if len(self.date_time) > self.length:
            self.date_time = file_name[-(self.split_char_loc+50):-40]
            
            
            

    

#Button to clear cache files and reselect - clearing all files

data_manager = ArayaManager()

#upload files - user select files to upload to server from local system / or other system
# files are bytes like therefore you read and open and close the file during running.
#unexpected behaviour like loss of data / loss of functions - especially int and string replacements. 
#uploaded files is a list of file like objects - not diretly readable by Pandas. 
with st.form("upload", clear_on_submit=True):
    uploaded_files = st.file_uploader("Choose a CSV file", type=['csv'], accept_multiple_files=True, key = 'key', help = 'select files by highlighting them then press upload file to start process' )
    submitted = st.form_submit_button("UPLOAD!",  help = 'upload files to start analysis!')



    
#loop through uploaded files to break open the list and pass through data IO stream 
for uploaded_file in uploaded_files:
    if submitted and uploaded_file is not None: # stops the error of 'no files and start function if there are uploaed files
        file_names = (uploaded_file.name)# to add to Araya Manager.
        
        data_manager.get_run_name(file_names) # get the Run_ID. access ArayaManager function by passing file name
        data_manager.get_date_time(file_names)
        data_manager.concatenate_dataframes(uploaded_file) # concat dataframes.
        
        
comp = data_manager.group_df #assign variable - contactenated dataframes - csv loading direct from path doesn't require this. 
print(comp)



#"""Too note - behaviour for Streamlit on webserver vs direct on python environment has a pecularity whereby the first lines of code for
#FAM_RFU / VIC_RFU / ROX_RFU will need to have .astype('int32) added and a change to stats_ROX potentially to stop integers passing through the loop
#- not sure why this behvaiour but leads to bug which indicates a problem with numpy attribute error - 'IntergerArray' object has no
#'reshape' streamlit"""

# start onwards with processing only if dataframe us greater than 0 
if len(comp) > 0:
    #conversion of string to float in ArayaManager needed - test adding in a coersion function - as per below. 
    #coerce mixed float / int nummbers from somewhere. Add to Araya Manager a method coerce_numeric.
    comp['FAM_RFU'] = comp['FAM_RFU'].astype('float').abs()
    comp['VIC_RFU'] = comp['VIC_RFU'].astype('float').abs()
    comp['ROX_RFU'] = comp['ROX_RFU'].astype('float').abs()
    #will remove this function - to the bottom - allow files to be processed and user defined and downloaed as whole set and analysed sets. 
    
    #process file attributes in to parameters for QC. Essential information. 
    comp['Well'] = comp['Row_ID']+comp['Col_ID']
    
    comp.sort_values(['Run_ID','Row_ID', 'Run_ID'], inplace = True)
    
    comp['norm_RNaseP'] =  comp['VIC_RFU'].abs() / comp['ROX_RFU']
    comp['norm_N_Cov'] =  comp["FAM_RFU"]  / comp['ROX_RFU']
    comp.index.names = ['order']
    comp.reset_index(inplace = True)
    comp['date_time'] = pd.to_datetime(comp['date_time'], format='%Y%m%d%H%M%S')
    comp[['date', 'time']] = comp['date_time'].astype(str).str.split(' ', 1, expand=True)
    
    e = len(comp)
    num_arrays = comp.Run_ID.nunique()
    date_min = comp.date.min()
    date_max = comp.date.max()
    controls = {'P19': 'A1500', 'O19' : 'A1500', 'O20': 'A1500',
                    'P21': 'NEG', 'O21': 'NEG', 'O22': 'NEG',
                    'O23':'S06', 'P23':'S06', 'O24': 'S06'}
       
    comp['control'] = comp['Well'].map(controls).fillna('patient')
else:
    st.warning('Please upload Araya files')
    st.stop()
    
def valid_array(df):
    vdf = df.groupby('Run_ID')['ROX_RFU'].mean().reset_index()
    vdf = vdf[vdf['ROX_RFU'] <= 1000] 
    vdf =  vdf['Run_ID'].to_list()
    df = df[~df.Run_ID.isin(vdf)]
    st.write('Arrays Exluded', str(vdf))
    return(df)    
    
 


def scoring(row):

    if row['norm_N_Cov'] < 3.0 and row['norm_RNaseP'] > 1.6:
        return('Negative Patient')
    elif row['norm_N_Cov'] > 3.0 and row['norm_N_Cov'] <= 9.0 and row['norm_RNaseP'] >1.1:
        return('PLOD')
    elif row['norm_N_Cov'] > 9.0 and row['norm_RNaseP'] >=1.0:
        return('N_Cov Patient Positive')
    elif row['norm_N_Cov'] > 9.0 and row['norm_RNaseP']<= 1.0:
        return('Control_N_Cov')
    elif row['norm_N_Cov'] <= 3.0 and row['norm_RNaseP'] <=1.59:
        return('No_Call')
    elif row['norm_N_Cov'] > 3.0 and row['norm_N_Cov'] <= 9.0 and row['norm_RNaseP'] <1.0:
        return'CTRL_PLOD'
    else:
        return('missing')

def void_detect_neg(row):
    if row['norm_N_Cov'] > 9.0:
        return('Detected')
    elif row['norm_N_Cov'] > 4.0 and row['norm_N_Cov'] <= 9.0:
        return('PLOD')
    else:
        return('negative')
def logit(row):
    if row['Detection'] == 'Detected':
        return(1)
    else:
        return(0)    
     


comp['Result'] = comp.apply(lambda row: scoring(row), axis = 1)   
comp['Detection'] = comp.apply(lambda row: void_detect_neg(row), axis = 1) 
comp['logit_r'] = comp.apply(lambda row: logit(row),axis = 1) 

mem = comp.groupby(['Result'])['Result'].agg(['count']).reset_index()
mem['Percentage'] = round(mem['count'] / e *100, 2)

detect = comp.groupby(['Run_ID', 'Detection'])['Detection'].count()
detect = detect.transpose()

st.subheader('Prevalence')
st.write('Number of array uploaded: ', str(num_arrays))
st.write('Number wells challenged: ', str(e))
st.write('date range min :', str(date_min), ' to date range max: ', str(date_max))
st.dataframe(mem)


####need a function here to strip out empty arrays from the data_stream - not a great idea
#####stripping data out- function should perhaps sit in loop and check the mean of ROX - if less than 1000
####delete the whole array from the dataframe - might be a point during ArayaManger to run mean of ROX and add colution to 
####filter out and delete before changing index to odrder
    
# these need to be wrapped in functions - currently outside functions for testing. Cache this and have it clear cache with
#above clear all files button.

#scoring - see if it is worth putting selection box for current cutoffs





#comp.to_csv('file_uploader_check.csv')

####start producing plots - use st.plotly_chart() container to view in app - they can be downloaded to html as well - read docs on use######## 

####use this pands view to check changes to the datafrme are correct.


#Display view of data for normalised N1N2 - order of processing over time with cutoff lines. 

def fam_pro_qc(comp):
    figN1 = px.scatter(comp, x= 'order', y = 'norm_N_Cov' ,color = 'Result', title = 'N1 N2 Calls - normalised processinng view')

    figN1.add_trace(go.Scatter(
        y=[10, 10],
        x=[comp.order.min(), comp.order.max()],
        mode="lines+markers+text",
        name="Lower_10_Positive_Boundary",
        text=["10"],
        textposition="top center",
        line=dict(color="red", dash = 'dash')
        ))

    figN1.update_traces(marker_size=6)

    figN1.add_trace(go.Scatter(
        y=[9, 9],
        x=[comp.order.min(), comp.order.max()],
        mode="lines+markers+text",
        name="Lower_9_Positive_Boundary",
        text=["9"],
        textposition="top center",
        line=dict(color="red", dash = 'dash')))



    figN1.update_xaxes(showgrid = True, gridwidth = 0.0002, gridcolor = 'grey')
    figN1.update_yaxes(range=[0, 20],showgrid = True, gridwidth = 0.0002, gridcolor = 'grey')
    st.plotly_chart(figN1, use_container_width=True)


def RP_pro_QC(comp):
    fig1bbnbb = px.scatter(comp, x= 'order', y = 'norm_RNaseP', color = 'Result', title = 'normalised RNaseP processing view' )
    fig1bbnbb.update_traces(marker_size=6)
    fig1bbnbb.update_yaxes(range=[0, 6])
    fig1bbnbb.add_trace(go.Scatter(
        y=[1.5, 1.5],
        x=[comp.order.min(), comp.order.max()],
        mode="lines+markers+text",
        name="Lower_1.5_RP Detected_Boundary",
        text=["1.5"],
        textposition="top center",
        line=dict(color="red", dash = 'dash')))
    fig1bbnbb.add_trace(go.Scatter(
        y=[2, 2],
        x=[comp.order.min(), comp.order.max()],
        mode="lines+markers+text",
        name="Lower_2_RP Detected_Boundary",
        text=["2"],
        textposition="top center",
        line=dict(color="blue", dash = 'dash')))
    st.plotly_chart(fig1bbnbb, use_container_width=True)





def plot_roxCV(comp):
    ROX_mean = round(comp.ROX_RFU.mean())
    ROX_std = round(comp.ROX_RFU.std())
    CV = round(((ROX_std/ROX_mean)*100),1)
    CT = ("CV% "+ str(CV), "Mean "+ str(ROX_mean), "standard deviation "+str(ROX_std))
    print(CT)

    fig3a = px.scatter(comp, x= comp.order, y = comp.ROX_RFU ,color = comp.Result)
    fig3a.update_traces(mode='markers', marker_line_width=0.01, marker_size=2)
    fig3a.add_trace(go.Scatter(
        y=[ROX_mean, ROX_mean],
        x=[comp.order.min(), comp.order.max()],    
        mode="lines+text",
        name="Mean average",
        text=["Mean"],
        textposition="top center",
        line=dict(color="blue", dash = 'dash')))
    fig3a.add_trace(go.Scatter(
        y=[ROX_mean + (ROX_std * 3), ROX_mean + (ROX_std * 3)],
        x=[comp.order.min(), comp.order.max()],    mode="lines+markers+text",
        name="UCL - Upper 3 SD Cutoff",
        text=["UCL"],
        textposition="top center",
        line=dict(color="Red", dash = 'dash')))

    fig3a.add_trace(go.Scatter(
        y=[1600, 1600],
        x=[comp.order.min(), comp.order.max()],    
        mode="lines+text",
        name="1600 ROX RFU LCL",
        text=["LCL"],
        textposition="top center",
        line=dict(color="Red", dash = 'dash')))
    fig3a.update_layout(title = 'ROX dispense processing ' + str(CT))
    fig3a.update_traces(marker_size=6)
    fig3a.update_yaxes(range=[0, 6000])
    
    st.plotly_chart(fig3a, use_container_width=True)
    
  

#Display ROX vs FAM plot for over chemical performace vs dispense
def roxfam(comp):
    figroxfam = px.scatter(comp, x= 'ROX_RFU', y = 'FAM_RFU' ,color = 'Result', title = 'N1 N2 detection Performance vs dispense')
    figroxfam.update_traces(marker_size=6)
    figroxfam.update_xaxes(range=[1000, 6000])
    figroxfam.update_yaxes(range=[0, 50000])


    figroxfam.add_trace(go.Scatter(
        x=[1600, 1600],
        y=[50000, -100],
        mode="lines",
        name="1600  RFU Lower Cutoff Limit",
        #text=["LCL"],
        text=["ROX 1600 lower cutoff"],
        textposition="top center",
        line=dict(color="Red", dash = 'dash')
        ))
    st.plotly_chart(figroxfam, use_container_width=True)

def cluster(comp):
    fig2b = px.scatter(comp, x= 'norm_RNaseP', y = 'norm_N_Cov',color = 'Result', title = 'Cluster Processing view')
    fig2b.update_xaxes(range=[-0.5, 4.5])
    fig2b.update_yaxes(range=[-0.5, 18])
    fig2b.update_traces(marker_size=6)
    st.plotly_chart(fig2b, use_container_width=True)






def ctrl_view(testdf, m, l , u):
    
    figdt = px.scatter(testdf, x='date_time', y='norm_N_Cov', color = 'Result')
    figdt.update_yaxes(range=[0, 20])
    figdt.update_traces(marker_size=6)
    figdt.add_trace(go.Scatter(
        y=[m, m],
        x=[testdf.date_time.min(), testdf.date_time.max()],
        mode="lines+markers+text",
        name="Val mean",
        text=["Mean"],
        textposition="top center",
        line=dict(color="red", dash = 'dash')))
    figdt.add_trace(go.Scatter(
        y=[l, l],
        x=[testdf.date_time.min(), testdf.date_time.max()],
        mode="lines+markers+text",
        name="3SD low",
        text=["-3SD"],
        textposition="top center",
        line=dict(color="yellow", dash = 'dash')))
    figdt.add_trace(go.Scatter(
        y=[u, u],
        x=[testdf.date_time.min(), testdf.date_time.max()],
        mode="lines+markers+text",
        name="3SD high",
        text=["+3SD"],
        textposition="top center",
        line=dict(color="yellow", dash = 'dash')))
    st.plotly_chart(figdt, use_container_width=True)
    
 
def ctrl_sig(testdf, sig, m, l , u):
    
    if sig == 'FAM_RFU':
        range = [0,60000]
    elif sig == 'VIC_RFU':
        range = [0,18000]
    
    figdt = px.scatter(testdf, x='date_time', y= sig, color = 'Result')
    figdt.update_yaxes(range=range)
    figdt.update_traces(marker_size=6)
    figdt.add_trace(go.Scatter(
        y=[m, m],
        x=[testdf.date_time.min(), testdf.date_time.max()],
        mode="lines+markers+text",
        name="Val mean",
        text=["Mean"],
        textposition="top center",
        line=dict(color="red", dash = 'dash')))
    figdt.add_trace(go.Scatter(
        y=[l, l],
        x=[testdf.date_time.min(), testdf.date_time.max()],
        mode="lines+markers+text",
        name="3SD low",
        text=["-3SD"],
        textposition="top center",
        line=dict(color="yellow", dash = 'dash')))
    figdt.add_trace(go.Scatter(
        y=[u, u],
        x=[testdf.date_time.min(), testdf.date_time.max()],
        mode="lines+markers+text",
        name="3SD high",
        text=["+3SD"],
        textposition="top center",
        line=dict(color="yellow", dash = 'dash')))
    st.plotly_chart(figdt, use_container_width=True)


#app layout - charts are already produced above - this allows arangment of charts in order to 
#order the chart layout as required. Place any further charts /  table genertion above this line.
#Place the st containers below this line to arrange them as required. leave tops headeders above all 
#scripts so that the buttons etc.. all continue to work with the flow of Streamlit. 

st.subheader('Nexar processing view - Trend Data')

fam_pro_qc(comp)

RP_pro_QC(comp)

plot_roxCV(comp)

st.subheader('Nexar processing - Clustering data view')


col1, col2 = st.columns(2)

with col1:
    roxfam(comp)
    
with col2:
    cluster(comp)
    

st.subheader('QC control data - mean and 3 +/- SD markers - taken from week 27 2021 post validation')


testac = comp[(comp.control == 'A1500')]

testso = comp[(comp.control == 'S06')]

testneg = comp[(comp.control == 'NEG')]

st.subheader('Accuplex nFAM')

ctrl_view(testac, 13.8, 12.8, 14.8)

st.subheader('Qnostics nFAM')

ctrl_view(testso, 12.9, 11.7, 14.1)

st.subheader('Accuplex FAM RFU')
 
ctrl_sig(testac, 'FAM_RFU', 37250, 21038, 53462) 

st.subheader('Qnostics FAM RFU')

ctrl_sig(testso, 'FAM_RFU', 35629, 24286, 46972)

 


#percentiles
def Q25(x):
    return x.quantile(0.25)

def Q50(x):
    return x.quantile(0.5)

def Q75(x):
    return x.quantile(0.75)

def ROXCV(df1):
    stats_ROX = df1.groupby(['Run_ID'])[['ROX_RFU']].agg(['count', 'mean','std','min',Q25, Q50, Q75, 'max'])

    CI95_hi_ROX = []
    CI95_lo_ROX = []
    CV_run_ROX = []
    for i in stats_ROX.index:
        c,m,s,t,u,q1,q2,v =(stats_ROX.loc[i])
        CI95_hi_ROX.append(m + 1.95*s/math.sqrt(c))
        CI95_lo_ROX.append(m - 1.95*s/math.sqrt(c))
        CV_run_ROX.append(s/m*100)
    
    stats_ROX['CI95% low ROX'] = CI95_lo_ROX
    stats_ROX['CI95% hi ROX'] = CI95_hi_ROX
    stats_ROX['ROX CV%'] = CV_run_ROX
    #stats_ROX = stats_ROX.reset_index()
    stats_ROX = stats_ROX.astype('int32')
    return(stats_ROX)

stats_ROX = ROXCV(comp)

    

def stats_FAM(df):
    
    df['FAM_RFU'] = df['FAM_RFU'].abs()
    stats_FAM = df.groupby(['Run_ID', 'Result'])['FAM_RFU'].agg(['count', 'mean','min', 'std', 'max']).fillna(0)
    

    
    CI95_hi_FAM = []
    CI95_lo_FAM = []
    CV_run_FAM = []


    for i in stats_FAM.index:
        c,m,s,t,v =(stats_FAM.loc[i])
        CI95_hi_FAM.append(m + 1.96*s/math.sqrt(c))
        CI95_lo_FAM.append(m - 1.96*s/math.sqrt(c))
        CV_run_FAM.append(100 - (s/m*100))

    stats_FAM['CI 95% low FAM'] = CI95_lo_FAM
    stats_FAM['CI 95_hi_FAM'] = CI95_hi_FAM
    stats_FAM['CV%_FAM'] = CV_run_FAM
    stats_FAM = stats_FAM.astype('int32')
    #stats_FAMM = stats_FAM.astype('int32')
    return(stats_FAM)

def stats_nFAM(df):
    
    
    stats_nFAM = df.groupby(['Run_ID', 'Result'])['norm_N_Cov'].agg(['count', 'mean','std', 'min', 'max']).astype(float)
    

    
    CI95_hi_nFAM = []
    CI95_lo_nFAM = []
    CV_run_nFAM = []


    for i in stats_nFAM.index:
        c,m,s,t,v =(stats_nFAM.loc[i])
        CI95_hi_nFAM.append(m + 1.96*s/math.sqrt(c))
        CI95_lo_nFAM.append(m - 1.96*s/math.sqrt(c))
        CV_run_nFAM.append(100 - (s/m*100))

    stats_nFAM['CI 95% low nFAM'] = CI95_lo_nFAM
    stats_nFAM['CI 95_hi_nFAM'] = CI95_hi_nFAM
    stats_nFAM['CV%_nFAM'] = CV_run_nFAM
    stats_nFAM = round(stats_nFAM, 1)
    #stats_nFAM['%Percent_detected'] = result['N1N2_detected'] / TOT*100
    return(stats_nFAM.fillna('-'))


def stats_CFO(df):
    
   
    stats_CFO = df.groupby(['Run_ID', 'Result'])['VIC_RFU'].agg(['count', 'mean','min', 'std', 'max']).fillna(0)
    

    
    CI95_hi_CFO = []
    CI95_lo_CFO = []
    CV_run_CFO = []


    for i in stats_CFO.index:
        c,m,s,t,v =(stats_CFO.loc[i])
        CI95_hi_CFO.append(m + 1.96*s/math.sqrt(c))
        CI95_lo_CFO.append(m - 1.96*s/math.sqrt(c))
        CV_run_CFO.append(s/m*100)

    stats_CFO['CI 95% low CFO'] = CI95_lo_CFO
    stats_CFO['CI 95_hi_CFO'] = CI95_hi_CFO
    stats_CFO['CV%_CFO'] = CV_run_CFO
    stats_CFO = stats_CFO.astype('int32')
    return(stats_CFO)
    
def stats_nCFO(df):
    
   
    stats_nCFO = df.groupby(['Run_ID', 'Result'])['norm_RNaseP'].agg(['count', 'mean','std', 'min', 'max']).astype(float)
    

    
    CI95_hi_nCFO = []
    CI95_lo_nCFO = []
    CV_run_nCFO = []


    for i in stats_nCFO.index:
        c,m,s,t,v =(stats_nCFO.loc[i])
        CI95_hi_nCFO.append(m + 1.96*s/math.sqrt(c))
        CI95_lo_nCFO.append(m - 1.96*s/math.sqrt(c))
        CV_run_nCFO.append(s/m*100)
        
    stats_nCFO['CI 95% low nCFO'] = CI95_lo_nCFO
    stats_nCFO['CI 95_hi_nCFO'] = CI95_hi_nCFO
    stats_nCFO['CV%_nCFO'] = CV_run_nCFO
    stats_nCFO = round(stats_nCFO, 1)
    #stats_nFAM['%Percent_detected'] = result['N1N2_detected'] / TOT*100
    return(stats_nCFO.fillna('-'))


  
@st.cache
def convert_df(df):
 # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv().encode('utf-8')

csv = convert_df(comp)




FAM_data = stats_FAM(comp)
CFO_data = stats_CFO(comp)
nFAM_data = stats_nFAM(comp)
nCFO_data = stats_nCFO(comp)
stats_ROX = ROXCV(comp)

st.subheader('ROX run stats')
st.dataframe(stats_ROX)
st.subheader('FAM run stats')
st.dataframe(FAM_data.astype(str))
st.subheader('CFO run stats')
st.dataframe(CFO_data.astype(str))
st.subheader('nFAM run stats')
st.dataframe(nFAM_data.astype(str))
st.subheader('nCFO run stats')
st.dataframe(nCFO_data.astype(str))



all_data = convert_df(comp)

CFO = convert_df(CFO_data)
     
nFAM = convert_df(nFAM_data)

nCFO = convert_df(nCFO_data)

ROX = convert_df(stats_ROX)

FAM = convert_df(FAM_data)

det = convert_df(detect)


st.sidebar.download_button(
    label="Download all data as CSV",
     data=all_data,
     file_name='araya_all_data.csv',
     mime='text/csv')
     


st.sidebar.download_button(
    label="Download FAM CSV",
     data=FAM,
     file_name='araya_FAM_out.csv',
     mime='text/csv')



st.sidebar.download_button(
    label="Download CFO CSV",
     data=CFO,
     file_name='araya_ox_CFO_out.csv',
     mime='text/csv')



st.sidebar.download_button(
    label="Download nFAM CSV",
     data=nFAM,
     file_name='araya_nFAM_out.csv',
     mime='text/csv')



st.sidebar.download_button(
    label="Download nCFO CSV",
     data=nCFO,
     file_name='araya_nCFO_out.csv',
     mime='text/csv')



st.sidebar.download_button(
    label="Download ROX CSV",
     data=ROX,
     file_name='araya_ROX_out.csv',
     mime='text/csv')

st.sidebar.download_button(
    label="Download detected CSV",
     data=det,
     file_name='araya_detected_out.csv',
     mime='text/csv')
