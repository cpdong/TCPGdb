from flask import Flask,render_template,request,url_for
from lifelines import KaplanMeierFitter;
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt;
import matplotlib.gridspec as gridspec;
import seaborn as sns;
from scipy import stats
import seaborn as sns;
import pandas as pd;
import numpy as np;
import base64, os,io,urllib;
import pymongo, subprocess;

#os.chdir('/Users/cpdong/Desktop/flaskapp/TTTT')

app = Flask(__name__)


try:
    client = pymongo.MongoClient(host = ["localhost:27017"], serverSelectionTimeoutMS = 2000)
    client.server_info() # will throw an exception
    
    print ("MongoDB connection is on, turn off and re-start mongodb server!")
    ### stop the mongodb server
    #mongo_stop = subprocess.Popen(['mongod', '--dbpath', '/Users/cpdong/Documents/devel/mongoData/db', '--shutdown']).communicate()[0]; # stop service on Linux
    mongo_stop = subprocess.Popen('kill $(pgrep mongo)', shell=True).communicate()[0]; # stop service on Mac
    ### restart the mongodb server
    mongo_start = subprocess.Popen(['mongod', '--dbpath', '/Users/cpdong/Documents/devel/mongoData',
                                        '--logpath', '/Users/cpdong/Documents/devel/mongoLog/mongo.log', '--fork']).communicate()[0]
except:
    print ("MongoDB connection is off, start mongodb server!")
    ### restart the mongodb server
    mongo_start = subprocess.Popen(['mongod', '--dbpath', '/Users/cpdong/Documents/devel/mongoData',
                                    '--logpath', '/Users/cpdong/Documents/devel/mongoLog/mongo.log', '--fork']).communicate()[0]


client = pymongo.MongoClient()
dblist = client.list_database_names()



@app.route('/', methods = ['GET', 'POST'])
def index():
    mgdb = client['test']
    mgcollect = mgdb['test']
    for document in mgcollect.find({}):
        print(document)
    
    return render_template('index.html',title_name = 'welcome')    

@app.route('/datasets', methods = ['GET', 'POST'])
def datasets():
   return render_template('datasets.html')
   
@app.route('/crispr', methods = ['GET', 'POST'])
def crispr():
    if request.method == "GET":
        message = request.values.get('study')
        if message == None:
            message = "Belk_CRISPRi_CD8T"

    data = {"up_data" : [['IIIIIIIDDDDDDDDDDDC1', 0.99,'0.01'],
        ['IIIIIIIDDDDDDDDDDDC2', 0.92,'0.01'],
        ['IIIIIIIDDDDDDDDDDDC3', 0.88,'0.01'],
        ['IIIIIIIDDDDDDDDDDDC4', 0.83,'0.01'],
        ['IIIIIIIDDDDDDDDDDDC15IIIIIIIDDDDDDDDDDDC15IIIIIIIDDDDDDDDDDDC15',0.77,'0.01']],
    "down_data" :[['IIIIIIIDDDDDDDDDDD1IIIIIIIDDDDDDDDDDD1IIIIIIIDDDDDDDDDDD1', -0.79,'0.01'],
        ['IIIIIIIDDDDDDDDDDD2', -0.88,'0.01'],
        ['IIIIIIIDDDDDDDDDDD3', -0.9,'0.01'],
        ['IIIIIIIDDDDDDDDDDD4', -0.93,'0.01'],
        ['AAAAAAA', -0.93,'0.01']]
    }
    
    return render_template('crispr.html', study=message, data=data)


@app.route('/tscore', methods = ['GET', 'POST'])
def tscore():
    showdiv = "none"
    gene = None
    if request.method == "POST":
        gene = request.form.get('t_gene')
        if gene is not None:
            showdiv='bock';
        else:
            showdiv='none';
            
    return render_template('tscore.html', gene=gene, showDiv=showdiv)



@app.route('/search', methods = ['GET', 'POST'])
def search():    
    #title = 'Random Data - Two Groups'
    #ylabel = r'Measurement, $m$ / units'
    #groupCol = 'group'
    #geneCol = 'Gene'
    #violin_plot(survdata, groupCol, geneCol, title)

    search_gene = None
    if request.method == "GET":
        search_gene = request.form.get('search')

        if search_gene !='' and search_gene is not None:
            print(search_gene)
            print('---------+++++++++++++++++####################')
            #kmplt1_data = Kaplan(survdata, genename='Gene', dataset="dataset X1")
            return render_template('quicksearch.html', **locals())
        else:
            return render_template('search.html')
    else:
        return render_template('search.html')
    
# route a complex url https://stackoverflow.com/questions/37990395/complex-routing-for-get-request-from-html-form-in-flask
@app.route("/details", methods = ['GET', 'POST'])
def qsearch():
    
    def boxplot_tcell_exp(data_list):
        # https://rowannicholls.github.io/python/graphs/ax_based/boxplots_significance.html
        pltNum = len(data_list);
        if pltNum > 0:
            fig, axs = plt.subplots(nrows=1, ncols=pltNum, figsize=(8*pltNum, 5), sharey=False)
            for i, (key, value) in enumerate(data_list.items()):
                dataname = key;
                data = value;
                labels = list(data.columns)
                bxp_stats = data.apply(lambda x: {'med':x.med, 'q1':x.q1, 'q3':x.q3, 'whislo':x['min'], 'whishi':x['max']}, axis=0).tolist()
                for index, item in enumerate(bxp_stats):
                    item.update({'label':labels[index]})

                axs[i] = plt.subplot(1,pltNum, i+1)
                bxp =axs[i].bxp(bxp_stats, widths=0.6, showfliers=False, patch_artist=True);
                axs[i].spines['top'].set_visible(False)
                axs[i].spines['right'].set_visible(False)
                plt.xticks(rotation=30)
                plt.setp(bxp['medians'], color='k')
                for patch, color in zip(bxp['boxes'], ['dodgerblue','orange','green','orangered','purple','brown','pink','gray','olive','cyan']):
                    patch.set_facecolor(color)
                plt.title(dataname)
        
            #Kaplan(survdata, genename='Gene', dataset="dataset X")
            img = io.BytesIO();
            plt.savefig(img, format='png');
            img.seek(0);
            img_png = base64.b64encode(img.getvalue());
            img_asc2code= img_png.decode('ascii');
            #print(img_asc2code)
            return img_asc2code
        else:
            return None
    
    # KM survival in Python
    def Kaplan(survdata_sets, genename): # required OS.time/OS/geneExprs columns
        genename = genename;
        newdata_list={}
        for key, value in survdata_sets.items(): # precheck the dataset
            sdata = value;
            sdata['group'] = np.where(sdata[genename] >= sdata[genename].quantile(0.5), 1, 0)
            if len(sdata['group'].unique()) == 2:
                newdata_list[key] = pd.DataFrame(sdata)
                
        pltNum = len(newdata_list);
        if pltNum > 0:
            fig, axs = plt.subplots(nrows=1, ncols=pltNum, figsize=(4*pltNum, 4), sharey=False)
            for i, (key, value) in enumerate(newdata_list.items()):
                dataname = key;
                survdata = value;
                survdata['group'] = np.where(survdata[genename] >= survdata[genename].quantile(0.5), 1, 0)
                survdata['OS.time'] = survdata['OS.time'] / 365;
                lrt = logrank_test(durations_A=survdata[survdata['group']==0]['OS.time'],
                                   durations_B=survdata[survdata['group']==1]['OS.time'],
                                   event_observed_A=survdata[survdata['group']==0]['OS'],
                                   event_observed_B=survdata[survdata['group']==1]['OS'])
                if lrt.p_value < 0.001:
                    lrt_pval = ' < 0.001';
                else:
                    lrt_pval = ' = ' + str(round(lrt.p_value, 3));
            
                #fig, (ax, ax2) = plt.subplots(1, 2)
                axs[i] = plt.subplot(1,pltNum, i+1)
                kmf = KaplanMeierFitter()
                kmf.fit(durations=survdata[survdata['group']==0]['OS.time'], event_observed=survdata[survdata['group']==0]['OS'],label="<median")
                kmf.plot(ax=axs[i],ci_show=False)
                kmf.fit(durations=survdata[survdata['group']==1]['OS.time'], event_observed=survdata[survdata['group']==1]['OS'],label=">=median")
                kmf.plot(ax=axs[i], ci_show=False)
                plt.ylim(0, 1);
                plt.text(0.05, 0.2, 'log-rank P' + lrt_pval, transform = axs[i].transAxes)
                if i==0: #only label first y-axis label
                    plt.ylabel('Survival probability')
                plt.xlabel('Times (year)')
                plt.title(dataname);
           #Kaplan(survdata, genename='Gene', dataset="dataset X")
            img = io.BytesIO();
            plt.savefig(img, format='png');
            img.seek(0);
            img_png = base64.b64encode(img.getvalue());
            img_asc2code= img_png.decode('ascii');
            #print(img_asc2code)
            return img_asc2code
        else:
            return None;

    def boxplot_cart(data_list):
        # https://rowannicholls.github.io/python/graphs/ax_based/boxplots_significance.html
        pltNum = len(data_list);
        if pltNum > 0:
            fig, axs = plt.subplots(nrows=1, ncols=pltNum, figsize=(4*pltNum, 4), sharey=False)
            for i, (key, value) in enumerate(data_list.items()):
                dataname = key.split('@')[0]; # key: GSE94555@0.01
                pval = key.split('@')[1];
                data = value;
                labels = list(data.columns)
                bxp_stats = data.apply(lambda x: {'med':x.med, 'q1':x.q1, 'q3':x.q3, 'whislo':x['min'], 'whishi':x['max']}, axis=0).tolist()
                for index, item in enumerate(bxp_stats):
                    item.update({'label':labels[index]})
                valrange = data.max().max() - data.min().min(); # get range for label pvalue
                valmax = data.max().max();
                bar_height = (valrange * 0.05 * 2) + valmax
                bar_tips = (valrange * 0.05 * 2) + valmax - (valrange * 0.02)
                text_height = (valrange * 0.05 * 2) + valmax + (valrange * 0.01)

                axs[i] = plt.subplot(1,pltNum, i+1)
                bxp =axs[i].bxp(bxp_stats, widths=0.6, showfliers=False, patch_artist=True);
                axs[i].spines['top'].set_visible(False)
                axs[i].spines['right'].set_visible(False)
                plt.xticks(rotation=30)
                plt.setp(bxp['medians'], color='k')
                for patch, color in zip(bxp['boxes'], ['pink', 'lightblue']):
                    patch.set_facecolor(color)
                plt.plot([1, 1, 2, 2], [bar_tips, bar_height, bar_height, bar_tips], lw=1, c='k')
                plt.text(1.5, text_height, 'p=' + str(pval), ha='center', c='k')
                plt.title(dataname)
        
            #Kaplan(survdata, genename='Gene', dataset="dataset X")
            img = io.BytesIO();
            plt.savefig(img, format='png');
            img.seek(0);
            img_png = base64.b64encode(img.getvalue());
            img_asc2code= img_png.decode('ascii');
            #print(img_asc2code)
            return img_asc2code
        else:
            return None
        
        
    def boxplot_autoimmune(bigdata_list):
        bigdata_list = { k: v for k, v in bigdata_list.items() if len(v) >0 }; # remove 0 length element
        gridNum = len(bigdata_list);
        if gridNum > 0:
            pltNum = sum(len(v) for k, v in bigdata_list.items())
            width_ratio = [len(v) for k, v in bigdata_list.items()]
            
            fig = plt.figure(figsize=(pltNum*2, 4))
            grid = plt.GridSpec(1, gridNum, width_ratios=width_ratio)
            
            for i, (key, value) in enumerate(bigdata_list.items()):
                disName = key;
                data_list = value;
                fake = fig.add_subplot(grid[i])
                fake.set_title(f'{disName}\n', fontweight='semibold', size=14)
                fake.set_axis_off()
                
                if len(data_list) > 0:
                    gs = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=len(data_list), subplot_spec=grid[i])
                    # https://stackoverflow.com/questions/27426668/row-titles-for-matplotlib-subplot
                    for j, (k, v) in enumerate(data_list.items()):
                        # k: GSE94555@0.01
                        dataname = k.split('@')[0];
                        pval = k.split('@')[1];
                        data = v;
                        labels = list(data.columns)
                        bxp_stats = data.apply(lambda x: {'med':x.med, 'q1':x.q1, 'q3':x.q3, 'whislo':x['min'], 'whishi':x['max']}, axis=0).tolist()
                        for index, item in enumerate(bxp_stats):
                            item.update({'label':labels[index]})
                        valrange = data.max().max() - data.min().min(); # get range for label pvalue
                        valmax = data.max().max();
                        bar_height = (valrange * 0.05 * 2) + valmax
                        bar_tips = (valrange * 0.05 * 2) + valmax - (valrange * 0.02)
                        text_height = (valrange * 0.05 * 2) + valmax + (valrange * 0.01)

                        axs = fig.add_subplot(gs[j])
                        bxp =axs.bxp(bxp_stats, widths=0.6, showfliers=False, patch_artist=True);
                        axs.spines['top'].set_visible(False)
                        axs.spines['right'].set_visible(False)
                        axs.set_title(dataname)
                        plt.xticks(rotation=30)
                        plt.setp(bxp['medians'], color='k')
                        for patch, color in zip(bxp['boxes'], ['pink', 'lightblue']):
                            patch.set_facecolor(color)
                        plt.plot([1, 1, 2, 2], [bar_tips, bar_height, bar_height, bar_tips], lw=1, c='k')
                        plt.text(1.5, text_height, 'p=' + str(pval), ha='center', c='k')
                        plt.title(dataname)
    
            img = io.BytesIO();
            plt.savefig(img, format='png');
            img.seek(0);
            img_png = base64.b64encode(img.getvalue());
            img_asc2code= img_png.decode('ascii');
            return img_asc2code
        else:
            return None
    
    survdata=pd.read_csv('/Users/cpdong/Library/CloudStorage/Dropbox/project/TCPGdb/data/Tlymphoma/testData.csv', header=0)
    survdata['group'] = np.where(survdata['Gene'] >= survdata['Gene'].quantile(0.5), 1, 0)
    #print(survdata)
    #kmplt1_data = Kaplan(survdata, genename='Gene', dataset="dataset X1")
    #kmplot1_data='PD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0idXRmLTgiPz4NCjwhLS0gR2VuZXJhdG9yOiBBZG9iZSBJbGx1c3RyYXRvciAxNi4wLjAsIFNWRyBFeHBvcnQgUGx1Zy1JbiAuIFNWRyBWZXJzaW9uOiA2LjAwIEJ1aWxkIDApICAtLT4NCjwhRE9DVFlQRSBzdmcgUFVCTElDICItLy9XM0MvL0RURCBTVkcgMS4xLy9FTiIgImh0dHA6Ly93d3cudzMub3JnL0dyYXBoaWNzL1NWRy8xLjEvRFREL3N2ZzExLmR0ZCI+DQo8c3ZnIHZlcnNpb249IjEuMSIgaWQ9IkxheWVyXzEiIHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8yMDAwL3N2ZyIgeG1sbnM6eGxpbms9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkveGxpbmsiIHg9IjBweCIgeT0iMHB4Ig0KCSB3aWR0aD0iMTI2cHgiIGhlaWdodD0iMTI2cHgiIHZpZXdCb3g9IjAgMCAxMjYgMTI2IiBlbmFibGUtYmFja2dyb3VuZD0ibmV3IDAgMCAxMjYgMTI2IiB4bWw6c3BhY2U9InByZXNlcnZlIj4NCjxnPg0KCTxyZWN0IHg9IjEuMDk1IiB5PSI5OC4yMjQiIHdpZHRoPSIxMjMuODEiIGhlaWdodD0iMTkuMjc1Ii8+DQoJPHJlY3QgeD0iMS4wOTUiIHk9Ijg1Ljc0IiB3aWR0aD0iMTIzLjgxIiBoZWlnaHQ9IjUuMjA1Ii8+DQoJPHBhdGggZD0iTTE4LjQwNCw5NS43MjFjMC43NjcsMCwxLjM4OS0wLjYyMywxLjM4OS0xLjM5cy0wLjYyMi0xLjM4OC0xLjM4OS0xLjM4OEgzLjQ4MWMtMC43NjcsMC0xLjM4OCwwLjYyMS0xLjM4OCwxLjM4OA0KCQlzMC42MjIsMS4zOSwxLjM4OCwxLjM5SDE4LjQwNHoiLz4NCgk8cGF0aCBkPSJNNDQuNDMzLDk1LjcyMWMwLjc2NywwLDEuMzg4LTAuNjIzLDEuMzg4LTEuMzlzLTAuNjIyLTEuMzg4LTEuMzg4LTEuMzg4SDI5LjUxYy0wLjc2NywwLTEuMzg5LDAuNjIxLTEuMzg5LDEuMzg4DQoJCXMwLjYyMiwxLjM5LDEuMzg5LDEuMzlINDQuNDMzeiIvPg0KCTxwYXRoIGQ9Ik03MC40NjEsOTUuNzIxYzAuNzY3LDAsMS4zODgtMC42MjMsMS4zODgtMS4zOXMtMC42MjItMS4zODgtMS4zODgtMS4zODhINTUuNTM5Yy0wLjc2NywwLTEuMzg4LDAuNjIxLTEuMzg4LDEuMzg4DQoJCXMwLjYyMiwxLjM5LDEuMzg4LDEuMzlINzAuNDYxeiIvPg0KCTxwYXRoIGQ9Ik05Ni40OSw5NS43MjFjMC43NjcsMCwxLjM4OS0wLjYyMywxLjM4OS0xLjM5cy0wLjYyMi0xLjM4OC0xLjM4OS0xLjM4OEg4MS41NjdjLTAuNzY3LDAtMS4zODgsMC42MjEtMS4zODgsMS4zODgNCgkJczAuNjIyLDEuMzksMS4zODgsMS4zOUg5Ni40OXoiLz4NCgk8cGF0aCBkPSJNMTIyLjUxOSw5NS43MjFjMC43NjcsMCwxLjM4OS0wLjYyMywxLjM4OS0xLjM5cy0wLjYyMi0xLjM4OC0xLjM4OS0xLjM4OGgtMTQuOTIzYy0wLjc2NywwLTEuMzg4LDAuNjIxLTEuMzg4LDEuMzg4DQoJCXMwLjYyMiwxLjM5LDEuMzg4LDEuMzlIMTIyLjUxOXoiLz4NCgk8cGF0aCBkPSJNNy40MSw4MC45aDUzLjQ0MmMwLjg2MywwLDEuNTYyLTAuNjk5LDEuNTYyLTEuNTYyVjM5LjU0M2MwLTAuODYyLTAuNjk5LTEuNTYzLTEuNTYyLTEuNTYzSDQ1LjMxNHYtNi41MzkNCgkJYzAtMC44NjEtMC42OTgtMS41NjItMS41NjEtMS41NjJIMjMuNDI4Yy0wLjg2MywwLTEuNTYyLDAuNy0xLjU2MiwxLjU2MnY2LjU0SDcuNDFjLTAuODYyLDAtMS41NjIsMC43LTEuNTYyLDEuNTYzdjM5Ljc5NQ0KCQlDNS44NDgsODAuMjAxLDYuNTQ3LDgwLjksNy40MSw4MC45eiBNMzQuNDkyLDU3Ljg3NGgtMS43OTZ2LTYuNzY4aDEuNzk2VjU3Ljg3NHogTTI2LjU2MywzNC41NzRoMTQuMDU1djMuNDA2SDI2LjU2M1YzNC41NzR6DQoJCSBNMTAuNTQ0LDQyLjY3OGg0Ny4xNzN2MTEuOThIMzYuOTQydi00LjAwNmMwLTAuODYzLTAuNjk5LTEuNTYzLTEuNTYyLTEuNTYzaC0zLjU4MmMtMC44NjMsMC0xLjU2MiwwLjY5OS0xLjU2MiwxLjU2M3Y0LjAwNg0KCQlIMTAuNTQ0VjQyLjY3OHoiLz4NCgk8cGF0aCBkPSJNNjguNzM0LDgwLjloNDkuOTU4YzAuODA3LDAsMS40Ni0wLjY1MywxLjQ2LTEuNDZWMTcuNTM0YzAtMC44MDYtMC42NTMtMS40NTktMS40Ni0xLjQ1OWgtMTQuNTI0VjkuOTYxDQoJCWMwLTAuODA3LTAuNjUzLTEuNDYtMS40Ni0xLjQ2aC0xOWMtMC44MDcsMC0xLjQ2LDAuNjUzLTEuNDYsMS40NnY2LjExNUg2OC43MzRjLTAuODA3LDAtMS40NiwwLjY1My0xLjQ2LDEuNDU5Vjc5LjQ0DQoJCUM2Ny4yNzQsODAuMjQ3LDY3LjkyNyw4MC45LDY4LjczNCw4MC45eiBNODYuNjM4LDEyLjg5aDEzLjEzOXYzLjE4Nkg4Ni42MzhWMTIuODl6Ii8+DQo8L2c+DQo8L3N2Zz4NCg=='
    #kmplot2 = Kaplan(survdata, genename='Gene', dataset="dataset X2")
    #kmplot3 = Kaplan(survdata, genename='Gene', dataset="dataset X3")
    df = pd.DataFrame({'Control':[10,3,5,9,7.0],'Treat':[11,4,6,7,6.5]})
    df.index = ['max','min','q1','q3','med']
    data_list1={'data-1@0.05': pd.DataFrame(df),
               'data-2@0.33': pd.DataFrame(df),
               'data-3@0.43': pd.DataFrame(df),
               'data-4@1': pd.DataFrame(df)}
    data_list2={'data-5@0.23': pd.DataFrame(df),
               'data-6@0.1': pd.DataFrame(df),
               'data-7@0.001': pd.DataFrame(df)}
    data_list3={'data-8@0.5': pd.DataFrame(df),
                'data-9@0.5': pd.DataFrame(df),
                'data-10@0.5': pd.DataFrame(df)}
    
    bigdata_list= {'MS':data_list1,
                   'IBD': data_list2,
                   'SLE': data_list3}
    
    
    df2 = pd.DataFrame({'celltype1':[10,3,5,9,7.0],'celltype2':[11,4,6,7,6.5], 'celltype3':[11,4,6,7,6.5], 'celltype4':[11,4,6,7,6.5],
                   'celltype5':[11,4,6,7,6.5], 'celltype6':[11,4,6,7,6.5], 'celltype7':[11,4,6,7,6.5]})
    df2.index = ['max','min','q1','q3','med']
    data_list2={'DICE': pd.DataFrame(df2),
           'AIBS': pd.DataFrame(df2)}
    
    
    if request.method == "GET":
        search_gene = request.form.get('search')
        print('you are here aaaaaaaaaaaaa')
        print(search_gene)
        survdata_sets={'datasetX1': pd.DataFrame(survdata),
               'datasetX2': pd.DataFrame(survdata),
                'datasetX3': pd.DataFrame(survdata)}
        tcell_exp_plt = boxplot_tcell_exp(data_list2)
        kmplt1_data = Kaplan(survdata_sets, genename='Gene')
        cart_plt =  boxplot_cart(data_list1);
        autoimmune_plt = boxplot_autoimmune(bigdata_list)
        #kmplt1_data = None
        #cart_plt =  None
        #autoimmune_plt = None
        print(search_gene)
        print('---------+++++++++++++++++####################')
        #return f'Hello {name} !'
        return render_template('quicksearch.html', **locals())


@app.route('/document', methods = ['GET', 'POST'])
def document():
   return render_template('document.html')

@app.route('/contact', methods = ['GET', 'POST'])
def contact():
   return render_template('contact.html')


 
if __name__ == '__main__':
	app.run(host='0.0.0.0', port=5000)