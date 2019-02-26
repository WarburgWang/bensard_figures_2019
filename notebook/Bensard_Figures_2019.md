
<b><u>BENSARD ET AL, 2019 NOTEBOOK</u></b>   
Paper: Loss of the Mitochondrial Pyruvate Carrier Promotes Tumor Initiation  
Paper Author: Claire L. Bensard, et. al.   
Code by Jordan A. Berg and Alex J. Bott   
Description: The following notebook contains the code required to replicate specified plots and supplements

<b><u>IMPORT DEPENDENCIES</u></b> 


```python
import os
import pandas as pd 
import xpresstools as xp

#Some data was also processed using Alt Analyze

#Set path to this directory for accessing and saving files
__path__  = os.getcwd() + '/'
```

<b><u>IMPORT HUMAN GSE8671 DATASET</u></b>   
The GEO-accessible data is not normalized. We therefore used Alt Analyze (http://www.altanalyze.org/) to RMA normalize probe signal.


```python
#Get data
df_GSE8671 = xp.get_df(__path__ + 'data/GSE8671_rma_normalized.csv') #RMA normalized with Alt Analyze
info_GSE8671 = xp.get_info(__path__ + "data/sample_info_gse8671.csv")
df_GSE8671_c = xp.keep_labels(df_GSE8671, info_GSE8671, label_list=['Normal','Adenoma'])
df_GSE8671_clean = xp.clean_df(df_GSE8671_c)

#Collapse multi-mapping probes
df_GSE8671_collapsed = xp.probe_collapse(df_GSE8671_clean, __path__ + "data/GPL570.txt")
df_GSE8671_collapsed.to_csv(__path__ + "data/collapsed_GSE8671.csv",sep=',')

#Scale dataset
df_GSE8671_scaled, df_GSE8671_labeled = xp.prep_data(df_GSE8671_collapsed, info_GSE8671)

#sort by sample labels
info_GSE8671_sorted = info_GSE8671.copy()
info_GSE8671_sorted = info_GSE8671_sorted.loc[info_GSE8671_sorted[1].isin(['Adenoma', 'Normal'])]
info_GSE8671_sorted = info_GSE8671_sorted.sort_values([1], ascending=False)
info_GSE8671_sorted_list = info_GSE8671_sorted[0].tolist()
df_GSE8671_scaled_sorted = df_GSE8671_scaled[info_GSE8671_sorted_list]

#set palette
gse8671_colors = {'Adenoma': (0.8705882352941177, 0.5607843137254902, 0.0196078431372549),
        'Normal': (0.00784313725490196, 0.6196078431372549, 0.45098039215686275)}
```


```python
xp.check_samples(df_GSE8671_clean)
```


![png](output_5_0.png)


<b><u>IMPORT HUMAN GSE20916 DATASET</u></b> 


```python
#Get data
df_GSE20916, info_GSE20916 = xp.get_geo('GSE20916')
info_GSE20916[1] = info_GSE20916[1].str.capitalize() #Make sample types look nice
info_GSE20916 = info_GSE20916.replace('Normal_colon', 'Normal')
df_GSE20916_c = xp.keep_labels(df_GSE20916, info_GSE20916, label_list=['Normal','Adenoma','Adenocarcinoma'])
df_GSE20916_clean = xp.clean_df(df_GSE20916_c)

#Collapse multi-mapping probes
df_GSE20916_collapsed = xp.probe_collapse(df_GSE20916_clean, __path__ + "data/GPL570.txt")
df_GSE20916_collapsed.to_csv(__path__ + "data/collapsed_GSE20916.txt",sep='\t')

#Scale sorted dataset
df_GSE20916_scaled, df_GSE20916_labeled = xp.prep_data(df_GSE20916_collapsed, info_GSE20916)

#sort by sample labels
info_GSE20916_sorted = info_GSE20916.copy()
info_GSE20916_sorted = info_GSE20916_sorted.loc[info_GSE20916_sorted[1].isin(['Adenoma', 'Adenocarcinoma','Normal'])]
info_GSE20916_sorted = info_GSE20916_sorted.sort_values([1], ascending=False)
info_GSE20916_sorted_list = info_GSE20916_sorted[0].tolist()
df_GSE20916_scaled_sorted = df_GSE20916_scaled[info_GSE20916_sorted_list]

#Drop Adenocarcinomas
df_GSE20916_collapsed_noac = xp.drop_label(df_GSE20916_collapsed, info_GSE20916, 'Adenocarcinoma')

gse20916_colors = {'Adenocarcinoma': (0.5725490196078431, 0.5843137254901961, 0.5686274509803921),
        'Adenoma': (0.8705882352941177, 0.5607843137254902, 0.0196078431372549),
        'Normal': (0.00784313725490196, 0.6196078431372549, 0.45098039215686275)}
```

    26-Feb-2019 15:54:45 DEBUG utils - Directory ./ already exists. Skipping.
    26-Feb-2019 15:54:45 INFO GEOparse - File already exist: using local version.
    26-Feb-2019 15:54:45 INFO GEOparse - Parsing ./GSE20916_family.soft.gz: 
    26-Feb-2019 15:54:45 DEBUG GEOparse - DATABASE: GeoMiame
    26-Feb-2019 15:54:45 DEBUG GEOparse - SERIES: GSE20916
    26-Feb-2019 15:54:45 DEBUG GEOparse - PLATFORM: GPL570
    /anaconda3/lib/python3.6/site-packages/GEOparse/GEOparse.py:84: DtypeWarning:
    
    Columns (2) have mixed types. Specify dtype option on import or set low_memory=False.
    
    26-Feb-2019 15:54:47 DEBUG GEOparse - SAMPLE: GSM523242
    26-Feb-2019 15:54:47 DEBUG GEOparse - SAMPLE: GSM523243
    26-Feb-2019 15:54:47 DEBUG GEOparse - SAMPLE: GSM523244
    26-Feb-2019 15:54:47 DEBUG GEOparse - SAMPLE: GSM523245
    26-Feb-2019 15:54:47 DEBUG GEOparse - SAMPLE: GSM523246
    26-Feb-2019 15:54:47 DEBUG GEOparse - SAMPLE: GSM523247
    26-Feb-2019 15:54:47 DEBUG GEOparse - SAMPLE: GSM523248
    26-Feb-2019 15:54:47 DEBUG GEOparse - SAMPLE: GSM523249
    26-Feb-2019 15:54:47 DEBUG GEOparse - SAMPLE: GSM523250
    26-Feb-2019 15:54:47 DEBUG GEOparse - SAMPLE: GSM523251
    26-Feb-2019 15:54:47 DEBUG GEOparse - SAMPLE: GSM523252
    26-Feb-2019 15:54:47 DEBUG GEOparse - SAMPLE: GSM523253
    26-Feb-2019 15:54:47 DEBUG GEOparse - SAMPLE: GSM523254
    26-Feb-2019 15:54:47 DEBUG GEOparse - SAMPLE: GSM523255
    26-Feb-2019 15:54:48 DEBUG GEOparse - SAMPLE: GSM523256
    26-Feb-2019 15:54:48 DEBUG GEOparse - SAMPLE: GSM523257
    26-Feb-2019 15:54:48 DEBUG GEOparse - SAMPLE: GSM523258
    26-Feb-2019 15:54:48 DEBUG GEOparse - SAMPLE: GSM523259
    26-Feb-2019 15:54:48 DEBUG GEOparse - SAMPLE: GSM523260
    26-Feb-2019 15:54:48 DEBUG GEOparse - SAMPLE: GSM523261
    26-Feb-2019 15:54:48 DEBUG GEOparse - SAMPLE: GSM523262
    26-Feb-2019 15:54:48 DEBUG GEOparse - SAMPLE: GSM523263
    26-Feb-2019 15:54:48 DEBUG GEOparse - SAMPLE: GSM523264
    26-Feb-2019 15:54:48 DEBUG GEOparse - SAMPLE: GSM523265
    26-Feb-2019 15:54:48 DEBUG GEOparse - SAMPLE: GSM523266
    26-Feb-2019 15:54:48 DEBUG GEOparse - SAMPLE: GSM523267
    26-Feb-2019 15:54:48 DEBUG GEOparse - SAMPLE: GSM523268
    26-Feb-2019 15:54:48 DEBUG GEOparse - SAMPLE: GSM523269
    26-Feb-2019 15:54:48 DEBUG GEOparse - SAMPLE: GSM523270
    26-Feb-2019 15:54:48 DEBUG GEOparse - SAMPLE: GSM523271
    26-Feb-2019 15:54:48 DEBUG GEOparse - SAMPLE: GSM523272
    26-Feb-2019 15:54:48 DEBUG GEOparse - SAMPLE: GSM523273
    26-Feb-2019 15:54:48 DEBUG GEOparse - SAMPLE: GSM523274
    26-Feb-2019 15:54:48 DEBUG GEOparse - SAMPLE: GSM523275
    26-Feb-2019 15:54:48 DEBUG GEOparse - SAMPLE: GSM523276
    26-Feb-2019 15:54:48 DEBUG GEOparse - SAMPLE: GSM523277
    26-Feb-2019 15:54:49 DEBUG GEOparse - SAMPLE: GSM523278
    26-Feb-2019 15:54:49 DEBUG GEOparse - SAMPLE: GSM523279
    26-Feb-2019 15:54:49 DEBUG GEOparse - SAMPLE: GSM523280
    26-Feb-2019 15:54:49 DEBUG GEOparse - SAMPLE: GSM523281
    26-Feb-2019 15:54:49 DEBUG GEOparse - SAMPLE: GSM523282
    26-Feb-2019 15:54:49 DEBUG GEOparse - SAMPLE: GSM523283
    26-Feb-2019 15:54:49 DEBUG GEOparse - SAMPLE: GSM523284
    26-Feb-2019 15:54:49 DEBUG GEOparse - SAMPLE: GSM523285
    26-Feb-2019 15:54:49 DEBUG GEOparse - SAMPLE: GSM523286
    26-Feb-2019 15:54:49 DEBUG GEOparse - SAMPLE: GSM523287
    26-Feb-2019 15:54:49 DEBUG GEOparse - SAMPLE: GSM523288
    26-Feb-2019 15:54:49 DEBUG GEOparse - SAMPLE: GSM523289
    26-Feb-2019 15:54:49 DEBUG GEOparse - SAMPLE: GSM523290
    26-Feb-2019 15:54:49 DEBUG GEOparse - SAMPLE: GSM523291
    26-Feb-2019 15:54:49 DEBUG GEOparse - SAMPLE: GSM523292
    26-Feb-2019 15:54:49 DEBUG GEOparse - SAMPLE: GSM523293
    26-Feb-2019 15:54:49 DEBUG GEOparse - SAMPLE: GSM523294
    26-Feb-2019 15:54:49 DEBUG GEOparse - SAMPLE: GSM523295
    26-Feb-2019 15:54:49 DEBUG GEOparse - SAMPLE: GSM523296
    26-Feb-2019 15:54:50 DEBUG GEOparse - SAMPLE: GSM523297
    26-Feb-2019 15:54:50 DEBUG GEOparse - SAMPLE: GSM523298
    26-Feb-2019 15:54:50 DEBUG GEOparse - SAMPLE: GSM523299
    26-Feb-2019 15:54:50 DEBUG GEOparse - SAMPLE: GSM523300
    26-Feb-2019 15:54:50 DEBUG GEOparse - SAMPLE: GSM523301
    26-Feb-2019 15:54:50 DEBUG GEOparse - SAMPLE: GSM523302
    26-Feb-2019 15:54:50 DEBUG GEOparse - SAMPLE: GSM523303
    26-Feb-2019 15:54:50 DEBUG GEOparse - SAMPLE: GSM523304
    26-Feb-2019 15:54:50 DEBUG GEOparse - SAMPLE: GSM523305
    26-Feb-2019 15:54:50 DEBUG GEOparse - SAMPLE: GSM523306
    26-Feb-2019 15:54:50 DEBUG GEOparse - SAMPLE: GSM523307
    26-Feb-2019 15:54:50 DEBUG GEOparse - SAMPLE: GSM523308
    26-Feb-2019 15:54:50 DEBUG GEOparse - SAMPLE: GSM523309
    26-Feb-2019 15:54:50 DEBUG GEOparse - SAMPLE: GSM523310
    26-Feb-2019 15:54:50 DEBUG GEOparse - SAMPLE: GSM523311
    26-Feb-2019 15:54:50 DEBUG GEOparse - SAMPLE: GSM523312
    26-Feb-2019 15:54:50 DEBUG GEOparse - SAMPLE: GSM523313
    26-Feb-2019 15:54:50 DEBUG GEOparse - SAMPLE: GSM523314
    26-Feb-2019 15:54:50 DEBUG GEOparse - SAMPLE: GSM523315
    26-Feb-2019 15:54:50 DEBUG GEOparse - SAMPLE: GSM523316
    26-Feb-2019 15:54:51 DEBUG GEOparse - SAMPLE: GSM523317
    26-Feb-2019 15:54:51 DEBUG GEOparse - SAMPLE: GSM523318
    26-Feb-2019 15:54:51 DEBUG GEOparse - SAMPLE: GSM523319
    26-Feb-2019 15:54:51 DEBUG GEOparse - SAMPLE: GSM523320
    26-Feb-2019 15:54:51 DEBUG GEOparse - SAMPLE: GSM523321
    26-Feb-2019 15:54:51 DEBUG GEOparse - SAMPLE: GSM523322
    26-Feb-2019 15:54:51 DEBUG GEOparse - SAMPLE: GSM523323
    26-Feb-2019 15:54:51 DEBUG GEOparse - SAMPLE: GSM523324
    26-Feb-2019 15:54:51 DEBUG GEOparse - SAMPLE: GSM523325
    26-Feb-2019 15:54:51 DEBUG GEOparse - SAMPLE: GSM523326
    26-Feb-2019 15:54:51 DEBUG GEOparse - SAMPLE: GSM523327
    26-Feb-2019 15:54:51 DEBUG GEOparse - SAMPLE: GSM523328
    26-Feb-2019 15:54:51 DEBUG GEOparse - SAMPLE: GSM523329
    26-Feb-2019 15:54:51 DEBUG GEOparse - SAMPLE: GSM523330
    26-Feb-2019 15:54:51 DEBUG GEOparse - SAMPLE: GSM523331
    26-Feb-2019 15:54:51 DEBUG GEOparse - SAMPLE: GSM523332
    26-Feb-2019 15:54:51 DEBUG GEOparse - SAMPLE: GSM523333
    26-Feb-2019 15:54:51 DEBUG GEOparse - SAMPLE: GSM523334
    26-Feb-2019 15:54:51 DEBUG GEOparse - SAMPLE: GSM523335
    26-Feb-2019 15:54:51 DEBUG GEOparse - SAMPLE: GSM523336
    26-Feb-2019 15:54:52 DEBUG GEOparse - SAMPLE: GSM523337
    26-Feb-2019 15:54:52 DEBUG GEOparse - SAMPLE: GSM523338
    26-Feb-2019 15:54:52 DEBUG GEOparse - SAMPLE: GSM523339
    26-Feb-2019 15:54:52 DEBUG GEOparse - SAMPLE: GSM523340
    26-Feb-2019 15:54:52 DEBUG GEOparse - SAMPLE: GSM523341
    26-Feb-2019 15:54:52 DEBUG GEOparse - SAMPLE: GSM523342
    26-Feb-2019 15:54:52 DEBUG GEOparse - SAMPLE: GSM523343
    26-Feb-2019 15:54:52 DEBUG GEOparse - SAMPLE: GSM523344
    26-Feb-2019 15:54:52 DEBUG GEOparse - SAMPLE: GSM523345
    26-Feb-2019 15:54:52 DEBUG GEOparse - SAMPLE: GSM523346
    26-Feb-2019 15:54:52 DEBUG GEOparse - SAMPLE: GSM523347
    26-Feb-2019 15:54:52 DEBUG GEOparse - SAMPLE: GSM523348
    26-Feb-2019 15:54:52 DEBUG GEOparse - SAMPLE: GSM523349
    26-Feb-2019 15:54:52 DEBUG GEOparse - SAMPLE: GSM523350
    26-Feb-2019 15:54:52 DEBUG GEOparse - SAMPLE: GSM523351
    26-Feb-2019 15:54:52 DEBUG GEOparse - SAMPLE: GSM523352
    26-Feb-2019 15:54:52 DEBUG GEOparse - SAMPLE: GSM523353
    26-Feb-2019 15:54:52 DEBUG GEOparse - SAMPLE: GSM523354
    26-Feb-2019 15:54:52 DEBUG GEOparse - SAMPLE: GSM523355
    26-Feb-2019 15:54:52 DEBUG GEOparse - SAMPLE: GSM523356
    26-Feb-2019 15:54:53 DEBUG GEOparse - SAMPLE: GSM523357
    26-Feb-2019 15:54:53 DEBUG GEOparse - SAMPLE: GSM523358
    26-Feb-2019 15:54:53 DEBUG GEOparse - SAMPLE: GSM523359
    26-Feb-2019 15:54:53 DEBUG GEOparse - SAMPLE: GSM523360
    26-Feb-2019 15:54:53 DEBUG GEOparse - SAMPLE: GSM523361
    26-Feb-2019 15:54:53 DEBUG GEOparse - SAMPLE: GSM523362
    26-Feb-2019 15:54:53 DEBUG GEOparse - SAMPLE: GSM523363
    26-Feb-2019 15:54:53 DEBUG GEOparse - SAMPLE: GSM523364
    26-Feb-2019 15:54:53 DEBUG GEOparse - SAMPLE: GSM523365
    26-Feb-2019 15:54:53 DEBUG GEOparse - SAMPLE: GSM523366
    26-Feb-2019 15:54:53 DEBUG GEOparse - SAMPLE: GSM523367
    26-Feb-2019 15:54:53 DEBUG GEOparse - SAMPLE: GSM523368
    26-Feb-2019 15:54:53 DEBUG GEOparse - SAMPLE: GSM523369
    26-Feb-2019 15:54:53 DEBUG GEOparse - SAMPLE: GSM523370
    26-Feb-2019 15:54:53 DEBUG GEOparse - SAMPLE: GSM523371
    26-Feb-2019 15:54:53 DEBUG GEOparse - SAMPLE: GSM523372
    26-Feb-2019 15:54:53 DEBUG GEOparse - SAMPLE: GSM523373
    26-Feb-2019 15:54:53 DEBUG GEOparse - SAMPLE: GSM523374
    26-Feb-2019 15:54:53 DEBUG GEOparse - SAMPLE: GSM523375
    26-Feb-2019 15:54:53 DEBUG GEOparse - SAMPLE: GSM523376
    26-Feb-2019 15:54:54 DEBUG GEOparse - SAMPLE: GSM523377
    26-Feb-2019 15:54:54 DEBUG GEOparse - SAMPLE: GSM523378
    26-Feb-2019 15:54:54 DEBUG GEOparse - SAMPLE: GSM523379
    26-Feb-2019 15:54:54 DEBUG GEOparse - SAMPLE: GSM523380
    26-Feb-2019 15:54:54 DEBUG GEOparse - SAMPLE: GSM523381
    26-Feb-2019 15:54:54 DEBUG GEOparse - SAMPLE: GSM523382
    26-Feb-2019 15:54:54 DEBUG GEOparse - SAMPLE: GSM523383
    26-Feb-2019 15:54:54 DEBUG GEOparse - SAMPLE: GSM523384
    26-Feb-2019 15:54:54 DEBUG GEOparse - SAMPLE: GSM523385
    26-Feb-2019 15:54:54 DEBUG GEOparse - SAMPLE: GSM523386


    Data processing summary:
    {'Cell intensity files were generated using GCOS (Affymetrix). The probe set data was generated using R/BioConductor (version 2.8.1) packages affy (version 1.20.2), gcrma (version 2.14.1), and FLUSH.LVS.bundle (version 1.2.1, proportion=0.6). For data filtration, we selected the probe sets with signal intensity above the threshold limit in at least 5% of samples. The threshold was established at the 98th percentile of the expression levels from Y-chromosome–linked probe set signals detected detectable in female samples', 'Cell intensity files were generated using GCOS (Affymetrix). The probe set data was generated using R/BioConductor (version 2.8.1) packages affy (version 1.20.2), gcrma (version 2.14.1), and FLUSH.LVS.bundle (version 1.2.1, proportion=0.6). For data filtration, we selected the probe sets with signal intensity above the threshold limit in at least 5% of samples. The threshold was established at the 98th percentile of the expression levels from Y-chromosome–linked probe set signals detected detectable in female samples. In addition, the probe sets with signal FC higher than 1.5 (in relation to median) in less than 6 samples were removed.'}



```python
xp.check_samples(df_GSE20916_clean)
```


![png](output_8_0.png)


<b><u>IMPORT MOUSE APC DATASET</u></b> 


```python
apc_df = pd.read_csv(__path__ + "data/custom_apc_normalized.csv", sep=",", index_col=0, low_memory=False)
apc_info = xp.get_info(__path__ + "data/apc_info.csv")
apc_clean = xp.clean_df(apc_df)
apc_scaled, apc_labeled = xp.prep_data(apc_clean, apc_info)

apc_colors = {'APC-KO': (0.00392156862745098, 0.45098039215686275, 0.6980392156862745),
         'APC-WT': (0.8705882352941177, 0.5607843137254902, 0.0196078431372549),
         'Normal-KO': (0.8352941176470589, 0.3686274509803922, 0.0),
         'Normal-WT': (0.00784313725490196, 0.6196078431372549, 0.45098039215686275)}
```

<b><u>IMPORT MOUSE AOMDSS DATASET</u></b> 


```python
aomdss_df = pd.read_csv(__path__ + "data/custom_aomdss_normalized.csv", sep=",", index_col=0, low_memory=False)
aomdss_info = xp.get_info(__path__ + "data/aomdss_info.csv")
aomdss_clean = xp.clean_df(aomdss_df)
aomdss_scaled, aomdss_labeled = xp.prep_data(aomdss_clean, aomdss_info)

aomdss_colors = {'AOMDSS-KO': (0.00392156862745098, 0.45098039215686275, 0.6980392156862745),
         'AOMDSS-WT': (0.8705882352941177, 0.5607843137254902, 0.0196078431372549),
         'Normal-KO': (0.8352941176470589, 0.3686274509803922, 0.0),
         'Normal-WT': (0.00784313725490196, 0.6196078431372549, 0.45098039215686275)}
```

<b><u>GET GENE SETS</u></b> 


```python
with open(__path__ + "data/custom22.csv", 'r') as f:
    for line in f:
        custom22_mouse = line.split(",")
custom22_human = [x.upper() for x in custom22_mouse]

stem_genes = ['Ctnnb1','Notch1','Ascl2','Myc','Hopx','Sox9','Lgr5','Lef1','Mmp7','Axin2','Cd44','Ccnd1','Bmi1','Tert']
diff_genes = ['Fabp2','Atoh1','Muc2','Krt20','Chga','Vil1','MPC1','MPC2']

pyru = pd.read_csv('data/custom_pyruvate_list.csv',header=None)
pyru_list = pyru[0].tolist()
```

<b><u>FIGURE 1</u></b> 


```python
"""
LDH and MPC Relationship in Human Dataset GSE8671
"""
xp.gene_overview(df_GSE8671_labeled, info_GSE8671, gene_name='MPC1',palette=gse8671_colors, 
                  order=['Normal','Adenoma'], grid=True, whitegrid=True,
                 save_fig=__path__ + 'plots/MPC1_boxswarm_GSE8671.pdf')
```


    <Figure size 432x288 with 0 Axes>



![png](output_16_1.png)



```python
xp.gene_overview(df_GSE8671_labeled, info_GSE8671, gene_name='MPC2',palette=gse8671_colors, 
                  order=['Normal','Adenoma'],
                 save_fig=__path__ + 'plots/MPC2_boxswarm_GSE8671.pdf', grid=True, whitegrid=True)
```


    <Figure size 432x288 with 0 Axes>



![png](output_17_1.png)



```python
xp.gene_overview(df_GSE8671_labeled, info_GSE8671, gene_name='LDHA',palette=gse8671_colors, 
                  order=['Normal','Adenoma'],
                 save_fig=__path__ + 'plots/LDHA_boxswarm_GSE8671.pdf', grid=True, whitegrid=True)
```


    <Figure size 432x288 with 0 Axes>



![png](output_18_1.png)



```python
xp.gene_overview(df_GSE8671_labeled, info_GSE8671, gene_name='LDHB',palette=gse8671_colors, 
                  order=['Normal','Adenoma'],
                 save_fig=__path__ + 'plots/LDHB_boxswarm_GSE8671.pdf', grid=True, whitegrid=True)
```


    <Figure size 432x288 with 0 Axes>



![png](output_19_1.png)



```python
"""
LDH and MPC Relationship in Human Dataset GSE20916
"""
xp.gene_overview(df_GSE20916_labeled, info_GSE20916, gene_name='MPC1',palette=gse20916_colors, 
                  order=['Normal','Adenoma','Adenocarcinoma'],
                 save_fig=__path__ + 'plots/MPC1_boxswarm_GSE20916.pdf', grid=True, whitegrid=True)
```


    <Figure size 432x288 with 0 Axes>



![png](output_20_1.png)



```python
xp.gene_overview(df_GSE20916_labeled, info_GSE20916, gene_name='MPC2',palette=gse20916_colors, 
                  order=['Normal','Adenoma','Adenocarcinoma'],
                 save_fig=__path__ + 'plots/MPC2_boxswarm_GSE20916.pdf', grid=True, whitegrid=True)
```


    <Figure size 432x288 with 0 Axes>



![png](output_21_1.png)



```python
xp.gene_overview(df_GSE20916_labeled, info_GSE20916, gene_name='LDHA',palette=gse20916_colors, 
                  order=['Normal','Adenoma','Adenocarcinoma'],
                 save_fig=__path__ + 'plots/LDHA_boxswarm_GSE20916.pdf', grid=True, whitegrid=True)
```


    <Figure size 432x288 with 0 Axes>



![png](output_22_1.png)



```python
xp.gene_overview(df_GSE20916_labeled, info_GSE20916, gene_name='LDHB',palette=gse20916_colors, 
                  order=['Normal','Adenoma','Adenocarcinoma'],
                 save_fig=__path__ + 'plots/LDHB_boxswarm_GSE20916.pdf', grid=True, whitegrid=True)
```


    <Figure size 432x288 with 0 Axes>



![png](output_23_1.png)


<b><u>FIGURE 1 SUPPLEMENT</u></b> 


```python
"""
Pyruvate-related Gene Expression Changes in Human Dataset GSE8671
"""
xp.heatmap(df_GSE8671_scaled_sorted, info_GSE8671, sample_palette=gse8671_colors, gene_list=pyru_list, 
            figsize=(14,16), save_fig=(__path__ + 'plots/GSE8671_customPyruvate_heatmap_colclustered.pdf'),
            row_cluster=True, col_cluster=True,
            cbar_kws={'label': 'z-score'})
```


    <Figure size 432x288 with 0 Axes>



![png](output_25_1.png)



```python
xp.heatmap(df_GSE8671_scaled_sorted, info_GSE8671, sample_palette=gse8671_colors, gene_list=pyru_list, 
            figsize=(14,16), save_fig=(__path__ + 'plots/GSE8671_customPyruvate_heatmap.pdf'),
            row_cluster=True, col_cluster=False,
            cbar_kws={'label': 'z-score'})
```


    <Figure size 432x288 with 0 Axes>



![png](output_26_1.png)



```python
"""
Pyruvate-related Gene Expression Changes in Human Dataset GSE20916
"""
xp.heatmap(df_GSE20916_scaled_sorted, info_GSE20916, sample_palette=gse20916_colors, gene_list=pyru_list, 
            figsize=(30,16), save_fig=(__path__ + 'plots/GSE20916_customPyruvate_heatmap_colclustered.pdf'),
            row_cluster=True, col_cluster=True,
            cbar_kws={'label': 'z-score'})
```


    <Figure size 432x288 with 0 Axes>



![png](output_27_1.png)



```python
xp.heatmap(df_GSE20916_scaled_sorted, info_GSE20916, sample_palette=gse20916_colors, gene_list=pyru_list, 
            figsize=(30,16), save_fig=(__path__ + 'plots/GSE20916_customPyruvate_heatmap.pdf'),
            row_cluster=True, col_cluster=False,
            cbar_kws={'label': 'z-score'})
```


    <Figure size 432x288 with 0 Axes>



![png](output_28_1.png)



```python
"""
Pyruvate-related Gene Expression Changes in Human Dataset GSE8671 -- PCA
"""
xp.pca(df_GSE8671_scaled, info_GSE8671, palette=gse8671_colors,
        gene_list=pyru_list,
       save_scree=__path__ + 'plots/GSE8671_pyru_list_PCA_scree.pdf',
       title='GSE8671 Pyruvate-related Genes PCA', 
        save_fig=__path__ + 'plots/GSE8671_pyru_list_2CI_PCA.pdf')
```


![png](output_29_0.png)



```python
"""
Overall Gene Expression Changes in Human Dataset GSE8671 -- PCA
"""
xp.pca(df_GSE8671_scaled, info_GSE8671, palette=gse8671_colors,
       save_scree=__path__ + 'plots/GSE8671_all_genes_PCA_scree.pdf',
       title='GSE8671 All Genes PCA', 
        save_fig=__path__ + 'plots/GSE8671_all_genes_2CI_PCA.pdf')
```


![png](output_30_0.png)



```python
"""
Pyruvate-related Gene Expression Changes in Human Dataset GSE20916 -- PCA
"""
xp.pca(df_GSE20916_scaled, info_GSE20916, palette=gse20916_colors, order_legend=[1,3,2],
        gene_list=pyru_list,
       save_scree=__path__ + 'plots/GSE20916_pyru_list_PCA_scree.pdf',
       title='GSE20916 Pyruvate-related Genes PCA', 
        save_fig=__path__ + 'plots/GSE20916_pyru_list_2CI_PCA.pdf')
```


![png](output_31_0.png)



```python
"""
Overall Gene Expression Changes in Human Dataset GSE20916 -- PCA
"""
xp.pca(df_GSE20916_scaled, info_GSE20916, palette=gse20916_colors, order_legend=[1,3,2],
       save_scree=__path__ + 'plots/GSE20916_all_genes_PCA_scree.pdf',
       title='GSE20916 All Genes PCA', 
        save_fig=__path__ + 'plots/GSE20916_all_genes_2CI_PCA.pdf')
```


![png](output_32_0.png)


<b><u>FIGURE 4</u></b> 


```python
"""
Mouse AOMDSS Heatmap Custom22 Probes
"""
xp.heatmap(aomdss_scaled, aomdss_info, sample_palette=aomdss_colors, gene_list=custom22_mouse, 
            figsize=(16,6.5), save_fig=(__path__ + 'plots/AOMDSS_custom22_heatmap.pdf'), 
            cbar_kws={'label': 'z-score'})
```


    <Figure size 432x288 with 0 Axes>



![png](output_34_1.png)



```python
"""
Mouse APC Heatmap Custom22 Probes
"""
xp.heatmap(apc_scaled, apc_info, sample_palette=apc_colors, gene_list=custom22_mouse, 
            figsize=(16,6.5), save_fig=(__path__ + 'plots/APC_custom22_heatmap.pdf'), 
            cbar_kws={'label': 'z-score'})
```


    <Figure size 432x288 with 0 Axes>



![png](output_35_1.png)



```python
"""
Human GSE8671 Heatmap Custom22 Probes
"""
xp.heatmap(df_GSE8671_scaled, info_GSE8671, sample_palette=gse8671_colors, gene_list=custom22_human, 
            figsize=(22,8), save_fig=(__path__ + 'plots/GSE8671_custom22_heatmap.pdf'), 
            cbar_kws={'label': 'z-score'})
```


    <Figure size 432x288 with 0 Axes>



![png](output_36_1.png)



```python
"""
Human GSE20916 Heatmap Custom22 Probes
"""
xp.heatmap(df_GSE20916_scaled, info_GSE20916, sample_palette=gse20916_colors, gene_list=custom22_human, 
            figsize=(40,8), save_fig=(__path__ + 'plots/GSE20916_custom22_heatmap.pdf'), 
            cbar_kws={'label': 'z-score'})
```


    <Figure size 432x288 with 0 Axes>



![png](output_37_1.png)


<b><u>FIGURE 4 SUPPLEMENT</u></b> 


```python
"""
Mouse AOMDSS Heatmap Custom22 Probes -- PCA
"""
xp.pca(aomdss_scaled, aomdss_info, palette=aomdss_colors, order_legend=[3,4,2,1], n_components=2,
        gene_list=custom22_human, whitegrid=True, grid=False,
       save_scree=__path__ + 'plots/AOMDSS_custom22_PCA_scree.pdf',
       title='AOMDSS Custom22 Probes PCA', 
        save_fig=__path__ + 'plots/AOMDSS_custom22_2CI_PCA.pdf')
```


![png](output_39_0.png)



```python
"""
Mouse APC Heatmap Custom22 Probes -- PCA
"""
xp.pca(apc_scaled, apc_info, palette=apc_colors, order_legend=[2,3,1,4], n_components=2,
        gene_list=custom22_human, whitegrid=True, grid=False,
       save_scree=__path__ + 'plots/APC_custom22_PCA_scree.pdf',
       title='APC Custom22 Probes PCA', 
        save_fig=__path__ + 'plots/APC_custom22_2CI_PCA.pdf')
```


![png](output_40_0.png)



```python
"""
Human GSE8671 Custom22 Probes -- PCA
"""
xp.pca(df_GSE8671_scaled, info_GSE8671, palette=gse8671_colors, n_components=2,
        gene_list=custom22_human, whitegrid=True, grid=False,
       save_scree=__path__ + 'plots/GSE8671_custom22_PCA_scree.pdf',
       title='GSE8671 Custom22 Probes PCA', 
        save_fig=__path__ + 'plots/GSE8671_custom22_2CI_PCA.pdf')
```


![png](output_41_0.png)



```python
"""
Human GSE20916 Custom22 Probes -- PCA
"""
xp.pca(df_GSE20916_scaled, info_GSE20916, palette=gse20916_colors, order_legend=[1,3,2], n_components=2,
        gene_list=custom22_human, whitegrid=True, grid=False,
       save_scree=__path__ + 'plots/GSE20916_custom22_PCA_scree.pdf',
       title='GSE20916 Custom22 Probes PCA', 
        save_fig=__path__ + 'plots/GSE20916_custom22_2CI_PCA.pdf')
```


![png](output_42_0.png)



```python
"""
Human GSE8671 Jointplots -- MPC1 vs AXIN2
"""
xp.jointplot(df_GSE8671_collapsed, info_GSE8671, 'MPC1', 'AXIN2', palette=gse8671_colors, 
              order=['Normal','Adenoma'],
              save_fig=__path__ + 'plots/GSE8671_MPC1_AXIN2_jointplot.pdf', 
              dpi=600, bbox_inches='tight', whitegrid=True, grid=True)
```

    /anaconda3/lib/python3.6/site-packages/scipy/stats/stats.py:1713: FutureWarning:
    
    Using a non-tuple sequence for multidimensional indexing is deprecated; use `arr[tuple(seq)]` instead of `arr[seq]`. In the future this will be interpreted as an array index, `arr[np.array(seq)]`, which will result either in an error or a different result.
    



    <Figure size 432x288 with 0 Axes>



![png](output_43_2.png)



```python
"""
Human GSE8671 Jointplots -- MPC1 vs SOX9
"""
xp.jointplot(df_GSE8671_collapsed, info_GSE8671, 'MPC1', 'SOX9', palette=gse8671_colors, 
              order=['Normal','Adenoma'],
              save_fig=__path__ + 'plots/GSE8671_MPC1_SOX9_jointplot.pdf', 
              dpi=600, bbox_inches='tight', whitegrid=True, grid=True)
```


    <Figure size 432x288 with 0 Axes>



![png](output_44_1.png)



```python
"""
Human GSE20916 Jointplots -- MPC1 vs AXIN2, no adenocarcinomas
"""
xp.jointplot(df_GSE20916_collapsed_noac, info_GSE20916, 'MPC1', 'AXIN2', palette=gse20916_colors, 
              order=['Normal','Adenoma'],
              save_fig=__path__ + 'plots/GSE20916_MPC1_AXIN2_jointplot.pdf', 
              dpi=600, bbox_inches='tight', whitegrid=True, grid=True)
```


    <Figure size 432x288 with 0 Axes>



![png](output_45_1.png)



```python
"""
Human GSE20916 Jointplots -- MPC1 vs SOX9, no adenocarcinomas
"""
xp.jointplot(df_GSE20916_collapsed_noac, info_GSE20916, 'MPC1', 'SOX9', palette=gse20916_colors, 
              order=['Normal','Adenoma'],
              save_fig=__path__ + 'plots/GSE20916_MPC1_SOX9_jointplot.pdf', 
              dpi=600, bbox_inches='tight', whitegrid=True, grid=True)
```


    <Figure size 432x288 with 0 Axes>



![png](output_46_1.png)



```python
"""
Human GSE20916 Jointplots -- MPC1 vs AXIN2
"""
xp.jointplot(df_GSE20916_collapsed, info_GSE20916, 'MPC1', 'AXIN2', palette=gse20916_colors, 
              order=['Normal','Adenoma'],
              save_fig=__path__ + 'plots/GSE20916_MPC1_AXIN2_jointplot.pdf', 
              dpi=600, bbox_inches='tight', whitegrid=True, grid=True)
```


    <Figure size 432x288 with 0 Axes>



![png](output_47_1.png)



```python
"""
Human GSE20916 Jointplots -- MPC1 vs SOX9
"""
xp.jointplot(df_GSE20916_collapsed, info_GSE20916, 'MPC1', 'SOX9', palette=gse20916_colors, 
              order=['Normal','Adenoma'],
              save_fig=__path__ + 'plots/GSE20916_MPC1_SOX9_jointplot.pdf', 
              dpi=600, bbox_inches='tight', whitegrid=True, grid=True)
```


    <Figure size 432x288 with 0 Axes>



![png](output_48_1.png)


<b><u>FIGURE 4 CORRELATIONS</u></b> 


```python
#GSE20916/MPC1/ALL
xp.linreg(df_GSE20916_collapsed, 'MPC1', __path__ + 'data/GSE20916_MPC1_ALL_correlations.csv')

#GSE20916/MPC2/ALL
xp.linreg(df_GSE20916_collapsed, 'MPC2', __path__ + 'data/GSE20916_MPC2_ALL_correlations.csv')

#GSE20916/MPC1/NoAdenocarcinoma
df_GSE20916_collapsed_noAC = xp.drop_label(df_GSE20916_collapsed, info_GSE20916, 'Adenocarcinoma')
xp.linreg(df_GSE20916_collapsed_noAC, 'MPC1', __path__ + 'data/GSE20916_MPC1_noAC_correlations.csv')

#GSE20916/MPC2/NoAdenocarcinoma
df_GSE20916_collapsed_noAC = xp.drop_label(df_GSE20916_collapsed, info_GSE20916, 'Adenocarcinoma')
xp.linreg(df_GSE20916_collapsed_noAC, 'MPC2', __path__ + 'data/GSE20916_MPC2_noAC_correlations.csv')

#GSE8671/MPC1/ALL
xp.linreg(df_GSE8671_collapsed, 'MPC1', __path__ + 'data/GSE8671_MPC1_ALL_correlations.csv')

#GSE8671/MPC2/ALL
xp.linreg(df_GSE8671_collapsed, 'MPC2', __path__ + 'data/GSE8671_MPC2_ALL_correlations.csv')
```


```python

```
