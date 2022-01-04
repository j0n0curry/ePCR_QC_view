

import streamlit as st
import os
from io import BytesIO
from io import StringIO
from io import TextIOWrapper
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

version = 'v1.1 - alpha'

def main():
    
    
    #Set up main page of application / Header info / data collection / file selection / remove files / Reset
    
    #Main landing page greating / info
    st.set_page_config(layout="wide")

    st.title('ePCR analysis tool ' +str(version))

    st.subheader("Upload Araya csv file directly for processing - either drag and drop or click to select files on your local machine")
    
    
    
    
    #Button to clear cache files and reselect - clearing all files
    if st.button('Clear Uploaded File(s)', help = 'press to clear all files and start a fresh session') and 'key' in st.session_state.keys():
        st.session_state.pop('key')
        st.experimental_rerun()
    
    #Select main function iinstatiat ArayaManager before for loop otherwise reads and closes files

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
    
    # start onwards with processing only if dataframe us greater than 0 
    if len(comp) > 0:
        #conversion of string to float in ArayaManager needed - test adding in a coersion function - as per below. 
        #coerce mixed float / int nummbers from somewhere. Add to Araya Manager a method coerce_numeric.
        comp['FAM_RFU'] = comp['FAM_RFU'].astype('float').astype('Int32')
        comp['VIC_RFU'] = comp['VIC_RFU'].astype('float').astype('Int32')
        comp['ROX_RFU'] = comp['ROX_RFU'].astype('float').astype('Int32')
        #will remove this function - to the bottom - allow files to be processed and user defined and downloaed as whole set and analysed sets. 
        
        #process file attributes in to parameters for QC. Essential information. 
        comp['Well'] = comp['Row_ID']+comp['Col_ID']
        comp['norm_RNaseP'] =  comp['VIC_RFU'].abs() / comp['ROX_RFU']
        comp['norm_N_Cov'] =  comp["FAM_RFU"]  / comp['ROX_RFU']
        comp.index.names=['order']
        comp.reset_index(inplace = True)
        comp['date_time'] = pd.to_datetime(comp['date_time'], format='%Y%m%d%H%M%S')
        #comp[['date', 'time']] = comp['date_time'].astype(str).str.split(' ', 1, expand=True)
        print(comp.norm_RNaseP)
        comp.to_csv('test_out.csv')
        
        controls = {'P19': 'A1500', 'O19' : 'A1500', 'O20': 'A1500',
                        'P21': 'NEG', 'O21': 'NEG', 'O22': 'NEG',
                        'O23':'S06', 'P23':'S06', 'O24': 'S06'}
           
        comp['control'] = comp['Well'].map(controls).fillna('paitent')
    def scoring(row):
    
        if row['norm_N_Cov'] < 3.0 and row['norm_RNaseP'] > 1.6:
            return('Negative Patient')
        elif row['norm_N_Cov'] > 3.0 and row['norm_N_Cov'] <= 10.0 and row['norm_RNaseP'] >1.1:
            return('PLOD')
        elif row['norm_N_Cov'] > 10.0 and row['norm_RNaseP'] >=1.0:
            return('N_Cov Paitent Positive')
        elif row['norm_N_Cov'] > 10.0 and row['norm_RNaseP']<= 1.0:
            return('Control_N_Cov')
        elif row['norm_N_Cov'] <= 3.0 and row['norm_RNaseP'] <=1.59:
            return('No_Call')
        elif row['norm_N_Cov'] > 3.0 and row['norm_N_Cov'] <= 10.0 and row['norm_RNaseP'] <1.0:
            return'CTRL_PLOD'
        else:
            return('missing')
    
    
    comp['Result'] = comp.apply(lambda row: scoring(row), axis = 1)   
    
    
        
    print(comp.head())
    st.table(comp.head())
    
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

        figN1.update_traces(marker_size=3)

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
        fig1bbnbb.update_traces(marker_size=3)
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
        fig3a.update_traces(marker_size=3)
        fig3a.update_yaxes(range=[0, 6000])
        
        st.plotly_chart(fig3a, use_container_width=True)
        
      
    
    #Display ROX vs FAM plot for over chemical performace vs dispense
    def roxfam(comp):
        figroxfam = px.scatter(comp, x= 'ROX_RFU', y = 'FAM_RFU' ,color = 'Result', title = 'N1 N2 detection Performance vs dispense')
        figroxfam.update_traces(marker_size=3)
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
        fig2b.update_traces(marker_size=3)
        st.plotly_chart(fig2b, use_container_width=True)
    
    
    
    
    
    
    #ctrl_qc_table = testdf.groupby(['date_time','Result'])['norm_N_Cov','FAM_RFU', 'ROX].agg('mean', 'std')
    #ctrl_qc_table = testdf.groupby(['date_time','Result'])['norm_N_Cov','FAM_RFU', 'ROX].agg('mean', 'std')
    
    #print(ctrl_qc_table)   
    
    def ctrl_view(testdf):
        
        figdt = px.scatter(testdf, x='date_time', y='norm_N_Cov', color = 'Result')
        figdt.update_yaxes(range=[0, 20])
        figdt.update_traces(marker_size=3)
        figdt.add_trace(go.Scatter(
            y=[12.5, 12.5],
            x=[testdf.date_time.min(), testdf.date_time.max()],
            mode="lines+markers+text",
            name="Val mean",
            text=["Mean"],
            textposition="top center",
            line=dict(color="red", dash = 'dash')))
        figdt.add_trace(go.Scatter(
            y=[10.2, 10.2],
            x=[testdf.date_time.min(), testdf.date_time.max()],
            mode="lines+markers+text",
            name="3SD low",
            text=["-3SD"],
            textposition="top center",
            line=dict(color="yellow", dash = 'dash')))
        figdt.add_trace(go.Scatter(
            y=[15.2, 15.2],
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
        
    
    st.subheader('Accuplex Control View - QC control data')
    
   


   
    testdf = comp[(comp.control == 'A1500')]
    
    ctrl_view(testdf)
    
    #st.dataframe(round(ctrl_qc_table),0)
    
    @st.cache
    def convert_df(df):
     # IMPORTANT: Cache the conversion to prevent computation on every rerun
        return df.to_csv().encode('utf-8')

    csv = convert_df(comp)

    st.download_button(
        label="Download process data as CSV",
         data=csv,
         file_name='araya_viewer.csv',
         mime='text/csv')
    
    

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
        
        


if __name__ == "__main__":
    main()
