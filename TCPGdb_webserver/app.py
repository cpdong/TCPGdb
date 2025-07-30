#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 24 10:39:52 2024

@author: cpdong
"""

from flask import Flask,render_template,request,url_for,render_template_string;
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
import json;

app = Flask(__name__)

mongodb_client = pymongo.MongoClient("yourmongodbtokens", serverSelectionTimeoutMS = 10000)
TCPGdb = mongodb_client['TCPGdb']

@app.route('/', methods = ['GET', 'POST'])
def index():
    try:
        mongodb_client.admin.command('ismaster')
        print('Atlas connection success!')
        return render_template('index.html', title_name = 'welcome')    
    except:
        print('Atlas connection failure!')
        return render_template_string('MongDb connection lost, please try again later!!')    

@app.route('/datasets', methods = ['GET', 'POST'])
def datasets():
   return render_template('datasets.html')
   
@app.route('/crispr', methods = ['GET', 'POST'])
def crispr():
    with open('./static/data/extdata/crispr.json') as f:
        crispr_summary = json.load(f)
    if request.method == "GET":
        message = request.values.get('study')
        if message == None:
            message = "Belk_CRISPRi_CD8T"
        else:
            message = message;

        mageck_sigNum_up = crispr_summary[message]['mageck_sigNum'][0]
        mageck_sigNum_down = crispr_summary[message]['mageck_sigNum'][1]
        
        gobp_sigNum_up = crispr_summary[message]['gobp_sigNum'][0]
        gobp_sigNum_down = crispr_summary[message]['gobp_sigNum'][1]
        go_bp_up_list = crispr_summary[message]['gobp_pos_list']
        go_bp_down_list = crispr_summary[message]['gobp_neg_list']
        gobp_totalNum = crispr_summary[message]['gobp_totalNum']
         
    return render_template('crispr.html', 
                           study=message, 
                           mageck_sigNum_up = mageck_sigNum_up,
                           mageck_sigNum_down = mageck_sigNum_down,
                           gobp_sigNum_up = gobp_sigNum_up,
                           gobp_sigNum_down = gobp_sigNum_down,
                           go_bp_up_list = go_bp_up_list,
                           go_bp_down_list = go_bp_down_list,
                           bpChartHeight =  max([gobp_totalNum * 45  + 55, 140]) )


@app.route('/score', methods = ['GET', 'POST'])
def score():
    showdiv = "none"
    gene = None;

    zlfc_cd8t_list = [];
    zlfc_cd4t_list = [];
    zlfc_treg_list = [];
    zlfc_cart_list = [];
    tscore_dash_cd8t = None,
    tscore_dash_cd4t = None,
    tscore_dash_treg = None,
    tscore_dash_cart = None,
    tscore_cd8t = [];
    tscore_cd4t =  [];
    tscore_treg =  [];
    realname = None;
    pert_method = 'activate';
    
    tscore_collection = TCPGdb['Tscore']
    with open('./static/data/extdata/boxplot_zlfc_frame.json') as f:
        boxplot_frame = json.load(f)
    
    with open('./static/data/extdata/tscore_top50.json') as f2:
        tscore_top50 = json.load(f2)

    searchable_tsdf = pd.read_csv('./static/data/extdata/TPS_genelist.tsv',header=0,sep='\t')
    searchable_tsgenelist = searchable_tsdf.to_numpy().flatten().tolist();


    tscore_top50_activate = tscore_top50['activate']
    cd8t_top50_activate = tscore_top50_activate['CD8T']
    cd4t_top50_activate = tscore_top50_activate['CD4T']
    treg_top50_activate = tscore_top50_activate['Treg']
    cart_top50_activate = tscore_top50_activate['CAR-T']
    
    tscore_top50_knockout = tscore_top50['knockout']
    cd8t_top50_knockout = tscore_top50_knockout['CD8T']
    cd4t_top50_knockout = tscore_top50_knockout['CD4T']
    treg_top50_knockout = tscore_top50_knockout['Treg']
    cart_top50_knockout = tscore_top50_knockout['CAR-T']
    
        
    if request.method == "POST":
        pert_method = request.form.get('pert_method')
        boxplot_frame_sele = boxplot_frame[pert_method]
    
        
        gene = request.form.get('t_gene')
        if gene in searchable_tsgenelist:
            if 'ENSG0' not in gene:

                gene = searchable_tsdf.loc[((searchable_tsdf['hgnc_symbol'] == gene) | (searchable_tsdf['mgi_symbol'] == gene)),
                                           'ensembl_gene_id'].tolist()[0];

                realname = searchable_tsdf.loc[(searchable_tsdf['ensembl_gene_id'] == gene), 'hgnc_symbol'].tolist()[0] ;

            try:
                tscore_collection.find({gene: {'$exists': 1}})[0];
                geneData = tscore_collection.find({gene: {'$exists': 1}})[0][ gene ]
                
                tscore_dash_cd8t = geneData['tscore']['score_' + pert_method][0];
                tscore_dash_cd4t = geneData['tscore']['score_' + pert_method][1];
                tscore_dash_treg = geneData['tscore']['score_' + pert_method][2];
                tscore_dash_cart = geneData['tscore']['score_' + pert_method][3];
                
                tscore_cd8t = geneData['tscore_sub']['score_' + pert_method]['CD8T']
                tscore_cd4t = geneData['tscore_sub']['score_' + pert_method]['CD4T']
                tscore_treg = geneData['tscore_sub']['score_' + pert_method]['Treg']
                
                boxplot_frame_sele = boxplot_frame[pert_method]
                # CD8T data get
                zlfc_cd8t_data = pd.DataFrame(boxplot_frame_sele['CD8T']).T
                zlfc_cd8t_data['zlfc'] = geneData['zlfc']['zlfc_' + pert_method]['CD8T'];
                zlfc_cd8t_data['gene'] = realname
                zlfc_cd8t_list = zlfc_cd8t_data.reset_index().values.tolist();
                # CD4T data get
                zlfc_cd4t_data = pd.DataFrame(boxplot_frame_sele['CD4T']).T
                zlfc_cd4t_data['zlfc'] = geneData['zlfc']['zlfc_' + pert_method]['CD4T'];
                zlfc_cd4t_data['gene'] = realname
                zlfc_cd4t_list = zlfc_cd4t_data.reset_index().values.tolist();
                # Treg data get
                zlfc_treg_data = pd.DataFrame(boxplot_frame_sele['Treg']).T
                zlfc_treg_data['zlfc'] = geneData['zlfc']['zlfc_' + pert_method]['Treg'];
                zlfc_treg_data['gene'] = realname
                zlfc_treg_list = zlfc_treg_data.reset_index().values.tolist();
                # CAR-T data get
                zlfc_cart_data = pd.DataFrame(boxplot_frame_sele['CART']).T
                zlfc_cart_data['zlfc'] = geneData['zlfc']['zlfc_' + pert_method]['CART'];
                zlfc_cart_data['gene'] = realname
                zlfc_cart_list = zlfc_cart_data.reset_index().values.tolist();
                
                showdiv='block';
            except:
                print("query gene not inside")
                showdiv='none';
        else:
            showdiv='none';
            
    return render_template('score.html', 
                           gene=gene,
                           showDiv=showdiv,
                           cd8t_top50_activate = cd8t_top50_activate,
                           cd4t_top50_activate = cd4t_top50_activate,
                           treg_top50_activate = treg_top50_activate,
                           cart_top50_activate = cart_top50_activate,
                           cd8t_top50_knockout = cd8t_top50_knockout,
                           cd4t_top50_knockout = cd4t_top50_knockout,
                           treg_top50_knockout = treg_top50_knockout,
                           cart_top50_knockout = cart_top50_knockout,
                           zlfc_cd8t_list = zlfc_cd8t_list,
                           zlfc_cd4t_list = zlfc_cd4t_list,
                           zlfc_treg_list = zlfc_treg_list,
                           zlfc_cart_list = zlfc_cart_list,
                           tscore_dash_cd8t = tscore_dash_cd8t,
                           tscore_dash_cd4t = tscore_dash_cd4t,
                           tscore_dash_treg = tscore_dash_treg,
                           tscore_dash_cart = tscore_dash_cart,
                           tscore_cd8t = tscore_cd8t,
                           tscore_cd4t = tscore_cd4t,
                           tscore_treg = tscore_treg,
                           realname = realname,
                           pert_method = pert_method )



@app.route('/search', methods = ['GET', 'POST'])
def search():    
    searchable_df = pd.read_csv('./static/data/extdata/tcellExpDb_genelist.tsv',header=0,sep='\t')
    searchable_genelist = searchable_df.to_numpy().flatten().tolist();
    search_gene = None
    if request.method == "GET":
        search_gene = request.values.get('search')
        if search_gene:
            print(search_gene)
    if search_gene in searchable_genelist:
        print(search_gene)
        if 'ENSG0' in search_gene:
            return render_template('quicksearch.html', **locals())
        else:
            search_gene = searchable_df.loc[searchable_df['hgnc_symbol'] == search_gene, 'ensembl_gene_id'].tolist()[0];
            #print('search_gene')
            return render_template('quicksearch.html', **locals())
    else:
        return render_template('search.html')
    
@app.route("/details", methods = ['GET', 'POST'])
def qsearch():
    
    def boxplot_tcell_exp(data_list):
        pltNum = len(data_list);
        if pltNum > 0:
            fig, axs = plt.subplots(nrows=1, ncols=pltNum, figsize=(8*pltNum, 5), sharey=False, dpi=72)
            for i, (key, value) in enumerate(data_list.items()):
                dataname = key;
                data = value;
                labels = list(data.columns)
                bxp_stats = data.apply(lambda x: {'med':x.med, 'q1':x.q1, 'q3':x.q3, 'whislo':x['min'], 'whishi':x['max']}, axis=0).tolist()
                for index, item in enumerate(bxp_stats):
                    item.update({'label':labels[index]})
                
                if pltNum>1:
                    axs[i] = plt.subplot(1,pltNum, i+1)
                    bxp =axs[i].bxp(bxp_stats, widths=0.6, showfliers=False, patch_artist=True);
                    axs[i].spines['top'].set_visible(False)
                    axs[i].spines['right'].set_visible(False)
                else:
                    axs = plt.subplot(1,pltNum, i+1)
                    bxp =axs.bxp(bxp_stats, widths=0.6, showfliers=False, patch_artist=True);
                    axs.spines['top'].set_visible(False)
                    axs.spines['right'].set_visible(False)
                    
                plt.xticks(rotation=30, ha='right')
                plt.setp(bxp['medians'], color='k')
                for patch, color in zip(bxp['boxes'], ['dodgerblue','orange','green','orangered','purple','brown','pink','gray','olive','cyan']):
                    patch.set_facecolor(color)
                plt.title(dataname)
        
            img = io.BytesIO();
            plt.subplots_adjust(bottom=0.3)
            plt.savefig(img, format='png');
            img.seek(0);
            img_png = base64.b64encode(img.getvalue());
            img_asc2code= img_png.decode('ascii');
            #print(img_asc2code)
            return img_asc2code
        else:
            return None
    
    def Kaplan(survdata_sets, genename): # required OS.time/OS/geneExprs columns
        genename = genename;
        newdata_list={}
        for key, value in survdata_sets.items(): # precheck the dataset
            sdata = value;
            if len(sdata['group'].unique()) == 2:
                newdata_list[key] = pd.DataFrame(sdata)
                
        pltNum = len(newdata_list);
        if pltNum > 0:
            fig, axs = plt.subplots(nrows=1, ncols=pltNum, figsize=(5*pltNum, 4), sharey=False, dpi=72)
            for i, (key, value) in enumerate(newdata_list.items()):
                dataname = key;
                survdata = value;

                lrt = logrank_test(durations_A=survdata[survdata['group']==0]['OS.time'],
                                   durations_B=survdata[survdata['group']==1]['OS.time'],
                                   event_observed_A=survdata[survdata['group']==0]['OS'],
                                   event_observed_B=survdata[survdata['group']==1]['OS'])
                if lrt.p_value < 0.001:
                    lrt_pval = ' < 0.001';
                else:
                    lrt_pval = ' = ' + str(round(lrt.p_value, 3));
            
                if pltNum>1:
                    axs[i] = plt.subplot(1,pltNum, i+1)
                kmf = KaplanMeierFitter()
                kmf.fit(durations=survdata[survdata['group']==0]['OS.time'], event_observed=survdata[survdata['group']==0]['OS'],label="<median")
                if pltNum >1:
                    kmf.plot(ax=axs[i],ci_show=False)
                else:
                    kmf.plot(ax=axs,ci_show=False)
                kmf.fit(durations=survdata[survdata['group']==1]['OS.time'], event_observed=survdata[survdata['group']==1]['OS'],label=">=median")
                if pltNum >1:
                    kmf.plot(ax=axs[i], ci_show=False)
                else:
                    kmf.plot(ax=axs, ci_show=False)
                plt.ylim(0, 1);
                if pltNum >1:
                    plt.text(0.05, 0.2, 'log-rank P' + lrt_pval, transform = axs[i].transAxes)
                else:
                    plt.text(0.05, 0.2, 'log-rank P' + lrt_pval, transform = axs.transAxes)
                if i==0: #only label first y-axis label
                    plt.ylabel('Survival probability')
                plt.xlabel('Times (year)')
                plt.title(dataname);

            plt.subplots_adjust(wspace=0.5)
            img = io.BytesIO();
            plt.savefig(img, format='png');
            img.seek(0);
            img_png = base64.b64encode(img.getvalue());
            img_asc2code= img_png.decode('ascii');
            return img_asc2code
        else:
            return None;

    def boxplot_cart(data_list):
        pltNum = len(data_list);
        if pltNum > 0:
            fig, axs = plt.subplots(nrows=1, ncols=pltNum, figsize=(3*pltNum, 6), sharey=False, dpi=72)
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

                if pltNum>1:
                    axs[i] = plt.subplot(1,pltNum, i+1)
                    bxp =axs[i].bxp(bxp_stats, widths=0.6, showfliers=False, patch_artist=True);
                    axs[i].spines['top'].set_visible(False)
                    axs[i].spines['right'].set_visible(False)
                else:
                    axs = plt.subplot(1,pltNum, i+1)
                    bxp =axs.bxp(bxp_stats, widths=0.6, showfliers=False, patch_artist=True);
                    axs.spines['top'].set_visible(False)
                    axs.spines['right'].set_visible(False)
                    
                plt.xticks(rotation=30)
                plt.setp(bxp['medians'], color='k')
                for patch, color in zip(bxp['boxes'], ['pink', 'lightblue']):
                    patch.set_facecolor(color)
                plt.plot([1, 1, 2, 2], [bar_tips, bar_height, bar_height, bar_tips], lw=1, c='k')
                plt.text(1.5, text_height, 'p=' + str(pval), ha='center', c='k')
                plt.title(dataname)
        
            #Kaplan(survdata, genename='Gene', dataset="dataset X")
            plt.subplots_adjust(left=0.1, right=0.9, bottom=0.15, wspace=0.3)
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
            
            fig = plt.figure(figsize=(pltNum*1.9, 4), dpi=72)
            grid = plt.GridSpec(1, gridNum, width_ratios=width_ratio)
            
            for i, (key, value) in enumerate(bigdata_list.items()):
                disName = key;
                data_list = value;
                fake = fig.add_subplot(grid[i])
                fake.set_title(f'{disName}\n', fontweight='semibold', size=14)
                fake.set_axis_off()
                
                if len(data_list) > 0:
                    gs = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=len(data_list), subplot_spec=grid[i])
                    for j, (k, v) in enumerate(data_list.items()):
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
                        plt.xticks(rotation=20)
                        plt.setp(bxp['medians'], color='k')
                        for patch, color in zip(bxp['boxes'], ['pink', 'lightblue']):
                            patch.set_facecolor(color)
                        plt.plot([1, 1, 2, 2], [bar_tips, bar_height, bar_height, bar_tips], lw=1, c='k')
                        plt.text(1.5, text_height, 'p=' + str(pval), ha='center', c='k')
                        plt.title(dataname)
            
            plt.subplots_adjust(left=0.05, right=0.95, wspace=0.2)
            img = io.BytesIO();
            plt.savefig(img, format='png');
            img.seek(0);
            img_png = base64.b64encode(img.getvalue());
            img_asc2code= img_png.decode('ascii');
            return img_asc2code
        else:
            return None
    
    # pre-load the survival information
    tcancer_collection = TCPGdb['Tcancer'];
    TcellExp_collection =  TCPGdb['TcellExp'];
    Tcoexp50gene_collection =  TCPGdb['Tcoexp50gene'];
    Tcoexp30pathway_collection =  TCPGdb['Tcoexp30pathway'];
    cartExp_collection =  TCPGdb['cartExp'];
    tAutoimmune_collection =  TCPGdb['tAutoimmune'];
    geneInfo_collection =  TCPGdb['geneInfo'];
    
    try:
        tcancer_collection.find({'survdata': {'$exists': 1}})[0];
        survdataSet = tcancer_collection.find({'survdata': {'$exists': 1}})[0][ 'survdata' ]
        TcellExp_cell_label = TcellExp_collection.find({'celltype': {'$exists': 1}})[0][ 'celltype' ]
        
    except:
        print("Check your databse connections!")
        return render_template('search.html')
    
    with open('./static/data/extdata/GO_BP_terms.json') as f:
        go_bp_terms = json.load(f)
    go_bp_terms = pd.DataFrame(go_bp_terms['Terms'], columns=['ID', 'GO'])
        
    searchable_df = pd.read_csv('./static/data/extdata/tcellExpDb_genelist.tsv',header=0,sep='\t')
    searchable_genelist = searchable_df.to_numpy().flatten().tolist();
    
    search_gene = request.values.get('search')
    if search_gene in searchable_genelist:
        if 'ENSG0' not in search_gene:
            realname = search_gene
            search_gene = searchable_df.loc[searchable_df['hgnc_symbol'] == search_gene, 'ensembl_gene_id'].tolist()[0];
        else:
            realname = searchable_df.loc[searchable_df['ensembl_gene_id'] == search_gene, 'hgnc_symbol'].tolist()[0];
    
        try:
            gene_info_data = geneInfo_collection.find({search_gene : {'$exists': 1}})[0][ search_gene ]
            entrezid = gene_info_data[2];
            if entrezid is not None:
                entrezid = str(int(entrezid))
            gene_type=gene_info_data[1];
            if gene_type == 0:
                gene_type = 'protein_coding'
            elif gene_type == 1:
                gene_type = 'lnRNA'
            elif gene_type == 2:
                gene_type = 'processed_pseudogene'
            elif gene_type == 3:
                gene_type = 'miRNA'
            gene_BP = gene_info_data[3];
            gene_CC = gene_info_data[4];
            gene_MF = gene_info_data[5];
        except:
            entrezid = gene_type = gene_BP = gene_CC =  gene_MF =None;
    
        try:
            geneBinary = tcancer_collection.find({search_gene : {'$exists': 1}})[0][ search_gene ]
            survData = {}
            for key, value in geneBinary.items():
                survData[key] =pd.DataFrame({'group': geneBinary[key], 'OS': survdataSet[key]['event'], 'OS.time': survdataSet[key]['time'] })
        except:
            survData = {};
        
        try:
            geneExpData = TcellExp_collection.find({search_gene : {'$exists': 1}})[0][ search_gene ]
            TexpData = {}
            for key, value in geneExpData.items():
                TexpData[key] = pd.DataFrame(geneExpData[key][0:], columns=['med', 'min', 'q1', 'q3', 'max'], index=TcellExp_cell_label[key]).T
        except:
            TexpData = {};
            
        try:
            coexpData = Tcoexp50gene_collection.find({search_gene : {'$exists': 1}})[0][ search_gene ]
            coexp_pos = coexpData['posGene']; coexp_neg = coexpData['negGene']; 
            coexpPathwayData = Tcoexp30pathway_collection.find({search_gene : {'$exists': 1}})[0][ search_gene ]
            coexpPathway_pos = pd.DataFrame(coexpPathwayData['posGset'][0:], columns=['GO','NES','pvalue']);
            coexpPathway_pos = pd.merge(go_bp_terms, coexpPathway_pos, on='GO', how='inner')
            coexpPathway_pos = coexpPathway_pos.drop('GO', axis=1);
            coexpPathway_pos = coexpPathway_pos.sort_values(by=['NES'], ascending=False);
            coexpPathway_pos_data = coexpPathway_pos.values.tolist();
            coexpPathway_neg = pd.DataFrame(coexpPathwayData['negGset'][0:], columns=['GO','NES','pvalue']);
            coexpPathway_neg = pd.merge(go_bp_terms, coexpPathway_neg, on='GO', how='inner')
            coexpPathway_neg = coexpPathway_neg.drop('GO', axis=1);
            coexpPathway_neg = coexpPathway_neg.sort_values(by=['NES'], ascending=False);
            coexpPathway_neg_data = coexpPathway_neg.values.tolist();
        except:
            coexp_pos = coexp_neg = coexpPathway_pos_data = coexpPathway_neg_data = []

        try:
            cartQuery = cartExp_collection.find({search_gene : {'$exists': 1}})[0][ search_gene ]
            cartExpData = {}
            studylist = []
            for key, value in cartQuery.items():
                for key2, value2 in value.items():
                    studylist.append(key2)
                    
            if 'GSE147991' in studylist:
                cartExpData['GSE147991_RNAseq@' + str(cartQuery['bulkRNAseq']['GSE147991'][0]) ] = pd.DataFrame({'Ctrl': cartQuery['bulkRNAseq']['GSE147991'][1:6],
                                                                                                          'GM-CS-KO': cartQuery['bulkRNAseq']['GSE147991'][6:11] },
                                                                                                         index=['med', 'min', 'q1', 'q3', 'max'])
            if 'GSE160154' in studylist:
                cartExpData['GSE160154_RNAseq@' + str(cartQuery['bulkRNAseq']['GSE160154'][0]) ] = pd.DataFrame({'Exhaust': cartQuery['bulkRNAseq']['GSE160154'][1:6],
                                                                                                          'Ctrl': cartQuery['bulkRNAseq']['GSE160154'][6:11] },
                                                                                                         index=['med', 'min', 'q1', 'q3', 'max'])
            if 'GSE173918' in studylist:
                cartExpData['GSE173918_RNAseq@' + str(cartQuery['bulkRNAseq']['GSE173918'][0]) ] = pd.DataFrame({'Ctrl': cartQuery['bulkRNAseq']['GSE173918'][1:6],
                                                                                                          'PRDM1-KO': cartQuery['bulkRNAseq']['GSE173918'][6:11] },
                                                                                                                index=['med', 'min', 'q1', 'q3', 'max'])

            if 'GSE175981' in studylist:
                cartExpData['GSE175981RNAseq@' + str(cartQuery['bulkRNAseq']['GSE175981'][0]) ] = pd.DataFrame({'Inactive': cartQuery['bulkRNAseq']['GSE175981'][1:6],
                                                                                                          'Active': cartQuery['bulkRNAseq']['GSE175981'][6:11] },
                                                                                                         index=['med', 'min', 'q1', 'q3', 'max'])
            if 'GSE223655' in studylist:
                cartExpData['GSE223655_RNAseq@' + str(cartQuery['bulkRNAseq']['GSE223655'][0]) ] = pd.DataFrame({'PD': cartQuery['bulkRNAseq']['GSE223655'][1:6],
                                                                                                          'cR': cartQuery['bulkRNAseq']['GSE223655'][6:11] },
                                                                                                                index=['med', 'min', 'q1', 'q3', 'max'])

            if 'GSE241456' in studylist:
                cartExpData['GSE241456_RNAseq@' + str(cartQuery['bulkRNAseq']['GSE241456'][0]) ] = pd.DataFrame({'Ctrl': cartQuery['bulkRNAseq']['GSE241456'][1:6],
                                                                                                          'NR4a-KO': cartQuery['bulkRNAseq']['GSE241456'][6:11] },
                                                                                                         index=['med', 'min', 'q1', 'q3', 'max'])
            if 'GSE245184' in studylist:
                cartExpData['GSE245184_RNAseq@' + str(cartQuery['bulkRNAseq']['GSE245184'][0]) ] = pd.DataFrame({'Ctrl': cartQuery['bulkRNAseq']['GSE245184'][1:6],
                                                                                                          'SUV39H1': cartQuery['bulkRNAseq']['GSE245184'][6:11] },
                                                                                                                index=['med', 'min', 'q1', 'q3', 'max'])

            if 'GSE248382' in studylist:
                cartExpData['GSE248382_RNAseq@' + str(cartQuery['bulkRNAseq']['GSE248382'][0]) ] = pd.DataFrame({'Control': cartQuery['bulkRNAseq']['GSE248382'][1:6],
                                                                                                          'D133p53': cartQuery['bulkRNAseq']['GSE248382'][6:11] },
                                                                                                         index=['med', 'min', 'q1', 'q3', 'max'])


            if 'GSE151511' in studylist:
                cartExpData['GSE151511_scRNAseq@' + str(cartQuery['scRNAseq']['GSE151511'][0]) ] = pd.DataFrame({'PD': cartQuery['scRNAseq']['GSE151511'][1:6],
                                                                                                          'CR_PR': cartQuery['scRNAseq']['GSE151511'][6:11] },
                                                                                                                        index=['med', 'min', 'q1', 'q3', 'max'])
            if 'GSE160160' in studylist:
                cartExpData['GSE160160_scRNAseq@' + str(cartQuery['scRNAseq']['GSE160160'][0]) ] = pd.DataFrame({'Exhaust': cartQuery['scRNAseq']['GSE160160'][1:6],
                                                                                                          'Ctrl': cartQuery['scRNAseq']['GSE160160'][6:11] },
                                                                                                                        index=['med', 'min', 'q1', 'q3', 'max'])
            if 'GSE243325' in studylist:
                cartExpData['GSE243325_scRNAseq@' + str(cartQuery['scRNAseq']['GSE243325'][0]) ] = pd.DataFrame({'PD': cartQuery['scRNAseq']['GSE243325'][1:6],
                                                                                                          'CR': cartQuery['scRNAseq']['GSE243325'][6:11] },
                                                                                                                        index=['med', 'min', 'q1', 'q3', 'max'])
            if 'GSE246960' in studylist:
                cartExpData['GSE246960_scRNAseq@' + str(cartQuery['scRNAseq']['GSE246960'][0]) ] = pd.DataFrame({'Ctrl': cartQuery['scRNAseq']['GSE246960'][1:6],
                                                                                                          'SUVKO_enhance': cartQuery['scRNAseq']['GSE246960'][6:11] },
                                                                                                                        index=['med', 'min', 'q1', 'q3', 'max'])
                
        except:
            cartExpData = {};
            
            
        try:
            autoimmuneQuery = tAutoimmune_collection.find({search_gene : {'$exists': 1}})[0][ search_gene ]
            autoimmuneData = {}
            for key, value in autoimmuneQuery.items():
                layer2 = {}
                for key2, value2 in value.items():
                    layer2[ key2 + '@' + str(value2[0]) ] = pd.DataFrame({'HC': value2[1:6], key : value2[6:11] },
                                                                   index=['med', 'min', 'q1', 'q3', 'max'])
                autoimmuneData[key] = layer2
        except:
            autoimmuneData = {};


        tcell_exp_plt = boxplot_tcell_exp(TexpData)
        kmplt1_data = Kaplan(survData, genename='Gene')
        cart_plt =  boxplot_cart(cartExpData);
        autoimmune_plt = boxplot_autoimmune(autoimmuneData)
        print(search_gene)
        return render_template('quicksearch.html', **locals())
    else:
        #return "input not valid!! \n Please check again!!!"
        return render_template('quicksearch_error.html', search_gene=search_gene)



@app.route('/document', methods = ['GET', 'POST'])
def document():
   return render_template('document.html')


@app.route('/contact', methods = ['GET', 'POST'])
def contact():
   return render_template('contact.html')


 
if __name__ == '__main__':
	app.run(host='127.0.0.1', port=5000)
