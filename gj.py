
import os
import sys
#from time import time
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.colors as mcolors
# from mpl_toolkits.axes_grid.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
from matplotlib.patches import Polygon
"""
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (12, 8),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
plt.rcParams.update(params)
"""
import time
import datetime
import math
import re
import collections
from collections import defaultdict,Counter
import itertools
import pandas as pd
from nested_dict import nested_dict
#import goatools_genelist_enrich
pd.set_option('max_rows', 10)
pd.set_option('expand_frame_repr', False)
from pandas.tools.plotting import table
from pyfasta import Fasta
from string import digits
from matplotlib_venn import venn3
from matplotlib_venn import venn2
import subprocess
import traceback
from multiprocessing import Pool
from inspect import getargvalues, stack
from scipy import stats
import scipy
import seaborn as sns
import json,pickle,cPickle
#import gffutils
from sklearn.cluster import KMeans
#from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn import metrics
from sklearn.decomposition import PCA
#sns.set(style="white", color_codes=True)
sns.set(style="darkgrid")
sns.set_context("poster")
#sns.set_context(font_scale=2.5)

#plt.style.use('ggplot')
#mpl.style.use('ggplot')




##################################################
### manuplate file path/name related
def get_project_scripts_pwd():
    scripts_pwd=os.getcwd()
    #project_pwd=scripts_pwd.replace("/scripts","")
    project_pwd=rreplace(scripts_pwd,"/scripts","")
    return project_pwd,scripts_pwd

def add_file_pwd(input_f,*args):
    project_pwd,scripts_pwd=get_project_scripts_pwd()
    if len(args) > 0:
        project_pwd=args[0]
    input_f = input_f if input_f.startswith("/Share") else project_pwd+'/'+input_f
    return input_f

def get_dir_files(input_dir):
    fn_ls = os.listdir(input_dir)
    fn_path_ls=[input_dir+'/'+fn for fn in fn_ls]
    return input_dir,fn_ls,fn_path_ls

##################################################
### time stamp realted
def timestamp():
    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    return st


##################################################
### matplotlib plot related 

def ax_legend_set_text_line_color(ax,color_ls=None):
    """ set legend text color as same as lines
    """
    leg = ax.legend(loc='best')
    # leg = plt.legend()
    for i,(line,text) in enumerate(zip(leg.get_lines(),leg.get_texts())):
        line_color = line.get_color()
        if color_ls is not None:
            line_color = color_ls[i]
        #text.set_color(line_color)
        plt.setp(text,color=line_color)


# plot dict as bar graph
def dict_bar_plot(a={'A>T': 1187, 'C>G': 6090, 'T>G': 2726},x_label="",y_label="",title_str="",save_fn="",add_h=0,percent=0,by=10,cal_ylim=0,text_rotation="vertical",fig_size_x=8,fig_size_y=6,x_rotation=45):
    x_pos=np.arange(len(a)) # x tick location
    tickLocation=x_pos
    width=0.8
    rectLocation=tickLocation-(width/2.0)
    #x_lab=sorted(a.keys())
    x_lab=sort_str_num_ls(a.keys())
    #print x_lab
    y_pos=[a[i] for i in x_lab]
    if percent > 0:
        y_pos=[i/float(by) for i in y_pos]
        y_pos=["%.1f"%(i) for i in y_pos]
    sub_all=sum(a.values())
    if save_fn == "":
        file_str=title_str+'.png'
    elif not save_fn.endswith(".png"):
        print "save fig should with .png end.(%s)"%(save_fn)
        sys.exit()
    else:
        file_str=save_fn
    if len(title_str.split("/")) > 3:
         title_str='/'.join(title_str.split("/")[-3:])
    title_str=title_str+' '+str(sub_all)
    w_all=len(a)  # 20-(8,6)
    w,h=8,6
    w = w*(w_all/20)
    h = h*(w_all/20)
    fig=plt.figure(1,figsize=(fig_size_x,fig_size_y), dpi=80)
    plt.bar(rectLocation,y_pos,width,alpha=0.4)
    ax=plt.gca()
    ax.set_xticklabels(x_lab, minor=False)
    ax.set_xticks(ticks=tickLocation, minor=False)
    #ax.set_yticks(ticks=[1,10,100,1000,10000,100000,1000000],minor=False)
    if cal_ylim:
        low = min(a.values())
        high=max(a.values())
        plt.ylim([math.ceil(low-0.5*(high-low)), math.ceil(high+0.5*(high-low))])
    plt.xticks(rotation=x_rotation)
    plt.xlabel(x_label) 
    plt.ylabel(y_label) 
    plt.title(title_str)
    for x,y in zip(tickLocation,y_pos):
        ax.text(x-(width/6.0),y,str(y),color='blue',rotation=text_rotation,size="medium")
    if add_h > 0:
        unit_n=4
        h_l_s=tickLocation[0:len(y_pos):unit_n]
        print "h_l_s: ",h_l_s
        h_l_e=tickLocation[unit_n-1::unit_n]
        print "h_l_e: ",h_l_e
        unit=len(y_pos)/unit_n
        print "num of unit: ",unit
        print y_pos
        unit_local_sum_ls=[sum(y_pos[i*unit_n:i*unit_n+unit_n]) for i in xrange(unit)]
        print "unit_sum_local: ",unit_local_sum_ls
        unit_local_max_ls=[max(y_pos[i*unit_n:i*unit_n+unit_n]) for i in xrange(unit)]
        print "unit_local_max: ",unit_local_max_ls
        for h_s,h_e,unit_local_max,unit_sum in zip(h_l_s,h_l_e,unit_local_max_ls,unit_local_sum_ls):
            print "add line: %s,%s,%s"%(unit_local_max+2,h_s,h_e)
            plt.hlines(y=max(y_pos)*1.1,xmin=h_s-width/2.0,xmax=h_e+width/2.0,color="red",linestyles='dashed')
            plt.axvline(x=h_e+width/2.0+width/8.0,ls='--')
            text_p_x=(h_s+h_e)/2.0
            text_p_y=max(y_pos)*1.1+10
            ax.text(text_p_x-width/2.0,text_p_y,str(unit_sum),color='red',size="medium")
    fig.tight_layout()
    plt.savefig(file_str)
    plt.close(fig)

# bar plot, truncated
def bar_plot_truncated(df,yLowMax=200,yHighMin=800,ax=None,savefn=None,figsize_x=16,figsize_y=8):

    status,yLowMax,yHighMin=get_df_trucated_val(df)

    if status == 0:
        if ax == None: f, ax = plt.subplots(figsize=(figsize_x,figsize_y))
        df.plot(kind='barh',ax=ax)
        if savefn != None: plt.savefig(savefn)
        return 

    if ax == None:
        f, ax = plt.subplots(figsize=(figsize_x,figsize_y))
    
    ax.axis('off')
    
    axis=[0,0]
    axis[0] = ax_add_sub_ax(fig=None,ax=ax,width="100%",height="30%",loc=1,axis_bgcol='white')
    axis[1] = ax_add_sub_ax(fig=None,ax=ax,width="100%",height="55%",loc=4,axis_bgcol='white')
    
    df.plot(kind='bar', ax=axis[0])
    df.plot(kind='bar', ax=axis[1])

    axis[0].set_ylim(yHighMin, )
    axis[1].set_ylim(0, yLowMax)
    axis[0].legend().set_visible(False)
    axis[1].legend().set_visible(False)
    
    plt.legend(bbox_to_anchor=(1.3, 1.1), bbox_transform=axis[0].transAxes)  # put legend on the right of ax[0]

    axis[0].spines['bottom'].set_visible(False)
    axis[1].spines['top'].set_visible(False)
    axis[0].xaxis.tick_top()
    axis[0].tick_params(labeltop='off')
    axis[1].xaxis.tick_bottom()
    d = .015
    kwargs = dict(transform=axis[0].transAxes, color='k', clip_on=False)
    axis[0].plot((-d,+d),(-d,+d), **kwargs)
    axis[0].plot((1-d,1+d),(-d,+d), **kwargs)
    kwargs.update(transform=axis[1].transAxes)
    axis[1].plot((-d,+d),(1-d,1+d), **kwargs)
    axis[1].plot((1-d,1+d),(1-d,1+d), **kwargs)
    
    if savefn != None:
        plt.savefig(savefn)
    plt.close()

def ax_bar_add_value(ax,orient='v',value_ls=None,map_color=0):
    for n,p in enumerate(ax.patches):
        x = p.get_bbox().get_points()[:,0]
        y = p.get_bbox().get_points()[1,1]
        z = p.get_bbox().get_points()[:,1]
        w = p.get_bbox().get_points()[1,0]

        if map_color :
            p_color = p.get_facecolor()
        else:
            p_color = "black"

        if orient == 'v':
            xloc = x.mean()
            yloc = y
            ha, va = 'center', 'bottom'
            if value_ls is None:
                loc_val = yloc
            else:
                loc_val = value_ls[n]
        elif orient == 'h':
            xloc = w
            yloc = z.mean()
            ha, va = 'left', 'center'
            if value_ls is None:
                loc_val = yloc
            else:
                loc_val = value_ls[n]
    
        ax.annotate('{:.0f}'.format(loc_val), (xloc, yloc), ha=ha, va=va, size='large', color=p_color)

        """
        if value_ls is None:
            ax.annotate('{:.0f}'.format(yloc), (xloc, yloc), ha=ha, va=ha, size='large')
            # text alighment center,bottom
        else:
            ax.annotate('{:.0f}'.format(value_ls[n]), (xloc, yloc), ha=ha, va=va, size='large')
        """

# pie plot   
def plot_ls_pie(labels="",val="",dic="",title_str="",file_str=None):
    fig=plt.figure(1,figsize=(8,8), dpi=80)
    file_str=title_str+'_pie.png'
    if len(title_str.split("/")) > 3:
         title_str='/'.join(title_str.split("/")[-3:])
    title_str=title_str+' '+str(sum(val))
    plt.pie(val,labels=labels,autopct='%1.1f%%',shadow=False)
    plt.title(title_str)
    fig.tight_layout()
    plt.savefig(file_str)
    plt.close(fig)

# plot dict as bar graph by converting to dataframe first
def dict_bar_plot_by_dataframe(d={'A>T': 1187, 'C>G': 6090, 'T>G': 2726},x_label="mutation type",y_label="",title_str="",save_fn="test22.png",add_h=0,percent=0,by=10,cal_ylim=0,text_rotation="vertical",fig_size_x=8,fig_size_y=6,x_rotation=45):
    all_num=sum(d.values())
    df=pd.DataFrame.from_dict(d,orient='index')
    df.columns=['count']
    df=df.sort_index(axis=0)
    print "[dict_bar_plot_by_dataframe]"
    print df
    fig=plt.figure(2,figsize=(fig_size_x,fig_size_y),dpi=80)
    df.plot(kind='barh',legend=False,figsize=(fig_size_x,fig_size_y),color='g')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title_str+' '+str(all_num)) 
    for i,label in enumerate(list(df.index)):
        s=df.ix[label]['count']
        y_pos=s
        plt.gca().annotate(str(s),(y_pos,i),color='r')
    plt.tight_layout()
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.savefig(save_fn)
    plt.close()
    return df

# read dict to df, then plot bar
def dict_bar_plot_by_dataframe_mutlti(d={'66':[1,2,3],'67':[6,2,3],'68':[4,3,4]},columns=["1-2","1&2","2-1"],x_label="sample",y_label="value",title_str="",save_fn="test22.png",orient='v',stacked=1,sum_all="",fig_size_x=8,fig_size_y=6):
    df=pd.DataFrame.from_dict(d,orient='index')
    if sum_all != "":
        sum_all=sum(sum(df.values))
    df.columns=columns
    df=df.sort_index(axis=0)
    print "[dict_bar_plot_by_dataframe_mutlti]"
    print df
    fig=plt.figure(2,figsize=(fig_size_x,fig_size_y))
    legend= 1 if df.shape[1] > 1 else 0
    if orient == 'h':
        df.plot(kind='barh',legend=legend,stacked=stacked,figsize=(fig_size_x,fig_size_y))
        plt.xlabel(y_label)
        plt.ylabel(x_label)
        plt.gca().yaxis.grid(False)
    if orient == 'v':
        df.plot(kind='bar',legend=legend,stacked=stacked,figsize=(fig_size_x,fig_size_y))
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.gca().xaxis.grid(False)
    plt.title(title_str+' '+str(sum_all))
    ax=plt.gca()
    for t in ax.get_xticklabels():
                #t.set(rotation=45,fontsize=15)
                t.set(rotation=45)

    #plt.gcf().subplots_adjust(bottom=0.15,right=0.95,left=0.1)
    plt.tight_layout()
    plt.savefig(save_fn)
    plt.close()
    return df

# venn plot
def venn3plot(mode,subsets_ls,labels_ls,title_str,save_fn,axis=None):
    fig=plt.figure(figsize=(10,10))
    if axis is not None:
        ax = axis
    else:
        fig,ax=plt.subplots()
    if mode == "string": # set1 = set(['A', 'B', 'C', 'D'])
        if len(subsets_ls) == 2: # 2 sets of string
            venn2(subsets = subsets_ls, set_labels = labels_ls,ax=ax)
        elif len(subsets_ls) == 3: # 3 sets of string
            venn3(subsets = subsets_ls, set_labels = labels_ls,ax=ax,set_colors=('b', 'g', 'r'))
        else:
            print "[Venn plot error] cannot handle string list len: %s"%(len(subsets_ls)) 
            return
    if mode == "number": # (1, 1, 1, 2, 1, 2, 2)
	    if len(subsets_ls) == 3: # 2 sets
		venn2(subsets = subsets_ls, set_labels = labels_ls)
	    elif len(subsets_ls) == 7: # 3 sets
		venn3(subsets = subsets_ls, set_labels = labels_ls)
	    else:
		print "[Venn plot error] cannot handle list len: %s"%(len(subsets_ls))
		return
    plt.title(title_str)
    if axis is not None: return axis
    plt.tight_layout()
    plt.savefig(save_fn)
    plt.close()
    return save_fn




# plot weblogo from fasta file
def plotFaWeblogo(input_fa="scripts/test_dataset/weblogo.fa",weblogo_f_title="", mode='dna'):
    #project_pwd,scripts_pwd=get_project_scripts_pwd()
    #input_fa = add_file_pwd(input_fa)
    if not input_fa.endswith('fa'):
        print "Input file should be fasta: %s"%(input_fa)
        sys.exit()
    weblogo="/Share/home/zhangqf/gongjing/software/weblogo/weblogo"
    #weblogo_f_str=input_fa.replace("fa","weblogo.png")
    weblogo_f_str=rreplace(input_fa,"fa","weblogo.png")
    if weblogo_f_title == "":
            seq_num = 0
	    with open(input_fa,'r') as IN:
		for line in IN:
		    if line.startswith(">"): seq_num += 1
	    weblogo_f_title=input_fa.split("/")[-1]+':'+str(seq_num)
    print "[plotFaWeblogo]"
    print "  - input_fa: %s"%(input_fa)
    print "  - weblogo_f_title: %s"%(weblogo_f_title)
    print "  - save to: %s"%(weblogo_f_str)
    subprocess.call(["%s -f %s -D fasta -o %s -F png -A %s -s large -U probability -c classic --title %s --reverse-stacks YES --stack-width 40.8 -n 200 "%("weblogo",input_fa,weblogo_f_str,mode,weblogo_f_title)],shell=True)    
    print 
    return weblogo_f_str

# display .png/.jpg file in specific axes
def ax_show_img(fig=None,ax=None,img_fn=None,savefn=None,aspect_auto=None):
    if ax == None:
        fig,ax=plt.subplots()
    if not img_fn:
        print "[error] img_fn not provided"
        sys.exit()
    img_fn=add_file_pwd(img_fn)
    printFuncRun('ax_show_img')
    printFuncArgs()
    image = mpimg.imread(img_fn)
    ax.axis("off")
    if aspect_auto:
        ax.imshow(image,aspect='auto')
    ax.imshow(image)
    if savefn and fig != None:
        fig.savefig(savefn)
        plt.close()
    printFuncRun('ax_show_img')
    print

def ax_add_sub_ax(fig=None,ax=None,width="50%",height="50%",loc=1,axis_bgcol='white',xticks=None,yticks=None):
    printFuncRun('ax_add_sub_ax')
    printFuncArgs()
    inset_ax = inset_axes(ax,
                    width=width, # width = 30% of parent_bbox
                    height=height, # height : 1 inch
                    loc=loc,
                    axes_kwargs={'axis_bgcolor':axis_bgcol,'frame_on':1,'visible':1})
    if xticks is None: inset_ax.set_xtick([])
    if yticks is None: inset_ax.set_ytick([])
    printFuncRun('ax_add_sub_ax')
    return inset_ax

# accumulate plot
def cumulate_dist_plot(ls_ls,ls_ls_label,bins=40,title=None,ax=None,savefn=None,xlabel=None,ylabel=None,add_vline=None,add_hline=None,log2transform=0,xlim=None,ylim=None):
    printFuncRun('cumulate_dist_plot')
    if ax is None:
        with sns.axes_style("ticks"):
            fig,ax = plt.subplots(figsize=(8,8))
    color_ls = sns_color_ls()
    ls_ls_label = [j+' ('+str(len(i))+')' for i,j in zip(ls_ls,ls_ls_label)]
    if log2transform:
        ls_ls = [np.log2(i) for i in ls_ls]
    for n,ls in enumerate(ls_ls):
        values,base = np.histogram(ls,bins=bins)
        cumulative = np.cumsum(values)
        cumulative_norm = [i/float(len(ls)) for i in cumulative]
        ax.plot(base[:-1],cumulative_norm,color=color_ls[n],label=ls_ls_label[n])
        print "plot line num: %s"%(n)
    ax.legend(loc='best')
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    else:
        ax.set_ylabel("Accumulate percent over total")
    if title is not None:
        ax.set_title(title)
    ax_legend_set_text_line_color(ax)
    if add_vline is not None:
        for vline in add_vline:
            ax.axvline(vline,ls="--", color='lightgrey')
    if add_hline is not None:
        for hline in add_hline:
            ax.axhline(hline,ls="--", color='lightgrey')
    if xlim is not None:
        ax.set_xlim(xlim[0], xlim[1])
    if ylim is not None:
        ax.set_ylim(ylim[0], ylim[1])
    ax.legend(loc="best")
    if savefn is not None:
        plt.savefig(savefn, dpi=200)
        plt.close()
    printFuncRun('cumulate_dist_plot')


##################################################
### seaborn plot

def sns_color_ls():
    return sns.color_palette("Set1", n_colors=8, desat=.5)*2
    #return sns.color_palette()

# violin plot
def df_sns_violinplot(df,col_str,savefn,orient='v',class_col="",order=None):
    assert savefn.endswith("png"), "[savefn error] should be .png: %s"%(savefn)
    savefn=savefn.replace('png','violin.png')

    print "[df_sns_violinplot]"
    col_str = None if col_str == "" else col_str
    class_col= None if class_col == "" else class_col
    print "  - col_str: %s"%(col_str)
    print "  - class_col: %s"%(class_col)
    print df.head(2)
    print df.dtypes

    #fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(8,8)) # plot multiple violin in one fig
    if class_col:
        s=sns.violinplot(y=col_str,x=class_col,data=df,order=order) # ax=ax[0,1,...]
    else:
        print "plot col_str: ",col_str
        #s=sns.violinplot(x='p_val_poisson',data=df) # ax=ax[0,1,...]
        s=sns.violinplot(df[col_str],order=order) # ax=ax[0,1,...]

    title_str=savefn.split("/")[-1]+':'+str(df.shape[0])
    fig=s.get_figure()
    fig.suptitle(title_str)
    fig.savefig(savefn)
    print "  - savefig: %s"%(savefn)
    plt.close()
    return savefn


# strip plot                
def df_sns_stripplot(df,col_str,savefn,orient='v',class_col=""):             
    assert savefn.endswith("png"), "[savefn error] should be .png: %s"%(savefn) 
    savefn=savefn.replace('png','strip.png')                                 
                             
    print "[df_sns_stripplot]"                                               
    col_str = None if col_str == "" else col_str                              
    class_col= None if class_col == "" else class_col                         
    print "  - col_str: %s"%(col_str)                                         
    print "  - class_col: %s"%(class_col)                                     
    print df.head(2)         
    print df.dtypes          
                             
    #fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(8,8)) # plot multiple violin in one fig                
    if class_col:            
        s=sns.stripplot(y=col_str,x=class_col,data=df,jitter=True) # ax=ax[0,1,...]      
    else:                    
        print "plot col_str: ",col_str                                        
        #s=sns.violinplot(x='p_val_poisson',data=df) # ax=ax[0,1,...]         
        s=sns.stripplot(df[col_str],jitter=True) # ax=ax[0,1,...]                        
                             
    title_str=savefn.split("/")[-1]+':'+str(df.shape[0])                      
    fig=s.get_figure()       
    fig.suptitle(title_str)  
    fig.savefig(savefn)      
    print "  - savefig: %s"%(savefn)                                          
    plt.close()              
    return savefn            

# swarm plot                
def df_sns_swarmplot(df,col_str,savefn,orient='v',class_col=""):
    assert savefn.endswith("png"), "[savefn error] should be .png: %s"%(savefn)
    savefn=savefn.replace('png','swarm.png')
            
    print "[df_sns_swarmplot]"
    col_str = None if col_str == "" else col_str
    class_col= None if class_col == "" else class_col
    print "  - col_str: %s"%(col_str)
    print "  - class_col: %s"%(class_col)
    print df.head(2)
    print df.dtypes
            
    #fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(8,8)) # plot multiple violin in one fig                
    if class_col:
        s=sns.swarmplot(y=col_str,x=class_col,data=df) # ax=ax[0,1,...]      
    else:   
        print "plot col_str: ",col_str
        #s=sns.violinplot(x='p_val_poisson',data=df) # ax=ax[0,1,...]         
        s=sns.swarmplot(df[col_str]) # ax=ax[0,1,...]                        
            
    title_str=savefn.split("/")[-1]+':'+str(df.shape[0])
    fig=s.get_figure()
    fig.suptitle(title_str)
    fig.savefig(savefn)
    print "  - savefig: %s"%(savefn)
    plt.close()
    return savefn


# jointplot
def df_sns_jointplot(col_str_x,col_str_y,savefn,df=None,list1='list1',list2='list2',xlim=None,ylim=None,x_y_lim_same=0,title_str=None,title_suptitle='right',use_scale_x_y_lim=1,color=None,xlabel=None,ylabel=None):
    print "[df_sns_jointplot]"
    assert savefn.endswith(("png","pdf")), "[savefn error] should be .png or .pdf: %s"%(savefn)
    #savefn=savefn.replace('png',col_str_x+'_'+col_str_y+'_joint.png')
    #savefn=rreplace(savefn,'png',str(col_str_x)+'_'+str(col_str_y)+'_joint.png')
    if df is None:
       print "df is None, generate df from supplied lists"
       print 'list1: %s, list2: %s'%(len(col_str_x),len(col_str_y))
       df=pd.DataFrame({list1:col_str_x,list2:col_str_y})
       col_str_x,col_str_y = list1,list2
       print df.info
       print df.head()
    s=sns.jointplot(x=col_str_x,y=col_str_y,data=df,size=8,kind='reg',xlim=xlim,ylim=ylim,color=color)
    #s=sns.jointplot(x=col_str_x,y=col_str_y,data=df,size=8,edgecolor="b")

    #ax_right = s.fig.gca() # get right axis
    #ax_right.set_title('')

    #ax_scatter = s.ax_joint.get_axes() # get axis of scatter

    a = s.fig.get_axes() # ax0,scatter; ax1,uppper; ax2:right
    ax_scatter = a[0]
    ax_upper = a[1]
    ax_right = a[2]


    xlim = (min(df[col_str_x]),max(df[col_str_x]))
    ylim = (min(df[col_str_y]),max(df[col_str_y]))
    if use_scale_x_y_lim:
       xlim = ax_scatter.get_xlim()
       ylim = ax_scatter.get_ylim()
    ax_scatter_min = min(xlim[0],ylim[0])
    ax_scatter_max = max(xlim[1],ylim[1])
   
    if x_y_lim_same>0:
        ax_scatter.set_xlim(ax_scatter_min,ax_scatter_max)
        ax_scatter.set_ylim(ax_scatter_min,ax_scatter_max)
        x = np.linspace(*ax_scatter.get_xlim())
        ax_scatter.plot(x, x,linestyle='--')  # add x=y line

    #s=(sns.jointplot(x=col_str_x,y=col_str_y,data=df,size=14).plot_joint(sns.regplot,x_jitter=.2,y_jitter=.2))
    if title_str is None:
        title_str=savefn.split("/")[-1]+':'+str(df.shape[0])
    #fig=s.get_figure() # 'JointGrid' object has no attribute 'get_figure'
    if title_suptitle == 'suptitle':
        ax_upper.set_title(title_str)
    if title_suptitle == 'title':
        ax_scatter.set_title(title_str)
    if title_suptitle == 'right':
        ax_right.set_title(title_str)
    if xlabel is not None: ax_scatter.set_xlabel(xlabel)
    if ylabel is not None: ax_scatter.set_ylabel(ylabel)
    s.savefig(savefn)
    plt.close()
    return savefn

# barplot
def df_sns_barplot(df,col_str_x,col_str_y,savefn,title_str,orient):
    s=sns.barplot(x=col_str_x,y=col_str_y,data=df,orient=orient)
    if title_str == "":
        title_str=savefn.split("/")[-1]+':'+str(df.shape[0])
    print "[df_sns_barplot]"
    fig=s.get_figure()
    fig.suptitle(title_str)
    fig.savefig(savefn,dpi=199)
    plt.close()
    return savefn

# pointplot
def df_sns_pointplot(col_str_x,col_str_y,savefn,title_str=None,df=None,linestyles=''): 
    s=sns.pointplot(x=col_str_x,y=col_str_y,data=df,linestyles=linestyles)
    point_num = df.shape[0] if df else len(col_str_x)
    if title_str:
        title_str=savefn.split("/")[-1]+':'+str(point_num)
    else:
        title_str = title_str + '('+str(point_num)+')'
    print "[df_sns_pointplot]"
    fig=s.get_figure()
    fig.suptitle(title_str)
    fig.savefig(savefn,dpi=199)
    plt.close()
    return savefn

# facet_grid
def df_sns_facetgrid(df="",col="",col_wrap="",row="",hue="",plot_col="",plot_type="plt.hist",sharex=False,sharey=True,size=4,aspect=1,palette="",margin_titles=1,hue_order="",row_order="",col_order="",savefn="",titles_ls="",suptitle_str=""):
    #for i in [col_wrap,hue,plot_col,palette,hue_order,row_order,col_order,col,row]:
        #i = str_empty_assign(i)
        #i = None if i == "" else i
    col = str_empty_assign(col)
    col_wrap = str_empty_assign(col_wrap)
    row = str_empty_assign(row)
    hue = str_empty_assign(hue)
    plot_col = str_empty_assign(plot_col)
    palette = str_empty_assign(palette)
    hue_order = str_empty_assign(hue_order)
    row_order = str_empty_assign(row_order)
    col_order = str_empty_assign(col_order)
    plot_type = str_empty_assign(plot_type)
    printFuncArgs()
    
    #g=sns.FacetGrid(data=df,col=col,col_wrap=eval(col_wrap),sharex=sharex,sharey=sharey,size=size,aspect=aspect,margin_titles=margin_titles,hue_order=hue_order,row_order=row_order,col_order=col_order,row=row)
    g=sns.FacetGrid(data=df,col=col,col_wrap=col_wrap,row=row,sharex=sharex,sharey=sharey,size=size,aspect=aspect,margin_titles=margin_titles,hue=hue,hue_order=hue_order,row_order=row_order,col_order=col_order)

    # if no col_order, works well;;; [BIG bug boss]
    # set col_order: fig, axes = plt.subplots(nrow, ncol, **kwargs); IndexError: index out of range

    if 'hist' in plot_type:
        g=g.map(eval(plot_type),plot_col,bins=np.arange(0,1.02,0.02),color='r') 
    if 'count' in plot_type:
        g=(g.map(eval(plot_type),plot_col).add_legend())
    if 'bar' in plot_type:
        g=g.map(eval(plot_type),plot_col)
    if titles_ls != "": 
        for ax, title in zip(g.axes.flat, titles_ls):
            ax.set_title(title)
    for ax in g.axes.flatten():
        for t in ax.get_xticklabels():
                #t.set(rotation=45,fontsize=15)
                t.set(rotation=45)
    plt.subplots_adjust(top=0.95)
    g.fig.suptitle(suptitle_str,fontsize=29)   
    g.savefig(savefn)
    plt.close()
    return savefn

# factorplot
def df_sns_factorplot(df,count_str,col_str,row_str,col_wrap,row_wrap,suptitle,savefn):
    s=sns.factorplot(count_str,col=col_str,col_wrap=col_wrap,data=df,kind='count',size=8,sharex=False)
    for ax in s.axes.flatten():
        for t in ax.get_xticklabels():
                t.set(rotation=45)
    plt.subplots_adjust(top=0.95)
    s.fig.suptitle(suptitle,fontsize=29)
    s.savefig(savefn)
    plt.close()

# matrix correlation plot called
def corrfunc(x, y, coor='pearson', **kws):
    if coor == 'pearson':
        #r, p = stats.pearsonr(list(np.exp2(x)), list(np.exp2(y)))
        r, p = stats.pearsonr(x, y)
    if coor == 'spearman':
        #r, p = stats.spearmanr(list(np.exp2(x)), list(np.exp2(y)))
        r, p = stats.spearmanr(x, y)
    ax = plt.gca()
    ax.annotate("{} r = {:.2f}\np = {:.2e}".format(coor,r,p),
                xy=(.5, .5),ha='center',va='center',fontsize=20, xycoords=ax.transAxes)

def sample_label(x,y,**kws):
    ax = plt.gca()
    ax.annotate('xy',xy=(.5,.5),ha='center',va='center',fontsize=20, xycoords=ax.transAxes)    

# matrix correlation plot
def df_corr_matrix_plot(df,savefn=None,size=4,rot=30,share_x_y=1,hue=None,diag='kde'):
    print df.info
    print df.describe() 
    #g = sns.PairGrid(df.fillna(0,inplace=1), palette=["red"],size=size)
    g = sns.PairGrid(df, palette=["red"],size=size,hue=hue)
    g.map_lower(plt.scatter, s=10)
    g.map_upper(corrfunc,coor='pearsonr')
    #g.map_diag(sns.distplot, hist=False, rug=True)
    if diag == 'kde':
        g.map_diag(sns.kdeplot)
    if diag == 'violin':
        g.map_diag(sns.violinplot)
    if diag == 'hist':
        g.map_diag(plt.hist)
    
    x_min_ls,x_max_ls,y_min_ls,y_max_ls = [],[],[],[]
    for ax in g.axes.flat:  
        plt.setp(ax.get_xticklabels(), rotation=rot)
        x_min_ls.append(ax.get_xlim()[0])
        x_max_ls.append(ax.get_xlim()[1])
        y_min_ls.append(ax.get_ylim()[0])
        y_max_ls.append(ax.get_ylim()[1])
    if share_x_y > 0:
       
        axis_min=min(min(x_min_ls),min(y_min_ls))
        axis_max=max(max(x_max_ls),max(y_max_ls))
        df_min=df.min().min()
        df_max=df.max().max() # df.values.max()
        df_min=df.describe().loc['min',:].min()
        df_max=df.describe().loc['max',:].max()
        for ax in g.axes.flat:
            ax.set_xlim(df_min,df_max)
            ax.set_ylim(df_min,df_max)
            #ax.set_xlim(axis_min,axis_max)
            #ax.set_ylim(axis_min,axis_max)
            #ax.set_xlim(min(x_min_ls),max(x_max_ls))
            #ax.set_ylim(min(y_min_ls),max(y_max_ls))

    """
    for i in xrange(df.shape[0]):
        label=g.axes[-1,i].get_xlabel()
        g.axes[-1,i].set_xlabel("") # set xlable invisiable
        g.axes[i,0].set_ylabel("")
        ax=g.axes[i,i]
        ax.set_title(label)
        ax.annotate("%s"%(label),
                    xy=(.5, .9),ha='center',va='center',fontsize=20, xycoords=ax.transAxes)
    """

    if savefn:
        g.savefig(savefn)
    plt.close()

# read df or from df_fn
def read_df_or_fn(df=None,df_fn=None,cols_keep=None,cols_drop=None):
    printFuncRun('read_df_or_fn')
    printFuncArgs()
    if df_fn is None and df is None:
        print "[error] both df & df_fn are None, please provide either of them!"
        return
    if df_fn is None:
        print "df_fn is None, use df directly!"
    if df is None:
        print "df is None, will read df from: %s"%(df_fn)
        df = pd.read_csv(df_fn,sep='\t',header=0)                 
    if cols_keep is not None:                            
        cols = df.columns  
        for col in cols_keep:   
            if col not in cols: 
                print "cols: %s not in columns"%(col),cols
                return                                   
        df = df[cols_keep]
    if cols_drop is not None:
        for col in cols_drop:
            del df[col]
    printFuncRun('read_df_or_fn')
    return df
   


# seaborn heatmap
def sns_heatmap(df=None,df_fn=None,cols_keep=None,cols_drop=None,savefn=None,figsize_x=80,figsize_y=100,cbar=True,square=True,annot=None,fmt='d',x_rotate=45,title=None, mask=None):
    printFuncRun('sns_heatmap')
    printFuncArgs()
    df = read_df_or_fn(df=df,df_fn=df_fn,cols_keep=cols_keep,cols_drop=cols_drop)
    rows,cols = df.shape
    #figsize_x = max(cols,8)
    #figsize_x = min(figsize_x,60)
    #figsize_y = max(rows/4,8)
    #if df.shape[0] > 1000:
    #    figsize_y = rows/200
    #figsize_y = min(figsize_y,180)
    print "figsize: (%s,%s)"%(figsize_x,figsize_y)
    figure = plt.figure(figsize=(figsize_x,figsize_y))
    annot = True if annot is None else annot
    ax = sns.heatmap(df.applymap(lambda x:int(x)),
                                linewidths=0.1, linecolor='white',annot=annot,fmt=fmt,cbar=cbar,square=square,mask=mask,
                                cbar_kws={"orientation": "vertical"})  # ,cmap="YlGnBu")
    ax.set_yticklabels(ax.yaxis.get_majorticklabels(), rotation=0) 
    ax.set_xticklabels(ax.xaxis.get_majorticklabels(), rotation=x_rotate) 
    if title is not None:
        ax.set_title(title)
    ht_f = ax.get_figure()
    plt.tight_layout()
    ht_f.savefig(savefn)
    plt.close()
    printFuncRun('sns_heatmap')

def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)
    #c = mcolors.ColorConverter().to_rgb
    #rvb = make_colormap([c('green'),c('black'), 0.5, c('black'), c('red'), 0.8, c('red')])


def sns_heatmap_annotate(df,cols_annotate=None,savefn=None):
    x, y = df.shape
    fig_size_x,fig_size_y = 8,8
    fig_size_y = x/7
    fig_size_x = y
        
    print fig_size_x,fig_size_y
    plt.subplots(figsize=(fig_size_x,fig_size_y))
    
    gs = gridspec.GridSpec(1, 2,width_ratios=[1,7])
    gs.update(wspace=0.35)
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    cols = df.columns
    print df.shape
    if cols_annotate is not None:
        cols_plot = [i for i in cols if i not in cols_annotate]
    df_plot = df[cols_plot]
    df_annotated = df[cols_annotate]
    
    annotated_cols_n = df_annotated.shape[1]
    print "haetmap cols:",cols_plot
    print "annotated cols(%s):"%(annotated_cols_n),cols_annotate
    print df_annotated
    
    #fig,ax = plt.subplots(1,2,figsize=(8,8))
    print 
    g1 = sns.heatmap(df_plot,ax=ax2,linewidths=0.1,yticklabels=False)
    
    df_annotated_scaled = df_annotated_scale(df_annotated)
    current_palette = sns.color_palette()
    g2 = sns.heatmap(df_annotated_scaled,ax=ax1,cbar=False,linewidths=0.01,cmap='Dark2',yticklabels=False)
    
   
    ax1.set_xticklabels(ax1.xaxis.get_majorticklabels(), rotation=45)
    ax1.get_yaxis().set_visible(False)

    plt.savefig(savefn)
    plt.close()

def df_annotated_scale(df):
    print df
    cols = df.columns
    dt = pd.DataFrame()
    for col in cols:
        col_value = list(df[col].value_counts().index)
        print col_value
        col_dict = dict(zip(col_value,range(1,len(col_value)+1)))
        print col_dict
        dt[col] = df[col].apply(lambda x:col_dict[x])
    print dt
    dt = dt.applymap(lambda x:int(x))
    return dt

# seaborn heatmap plot, along with bar plot for each axis
def sns_plot_heatmap_bar(df,input_f="",f_size_x=None,f_size_y=None,ratio=3,title='',top_bar_x=1,right_bar_y=1,savefn=None,pct_by_row=0,pct_by_col=0,ratio_w=3,ratio_h=3):
        df=df.fillna(0)
        print df
        if input_f != "":
                df=pd.read_csv(input_f)
        column_labels=list(df.columns)
        row_labels=list(df.index)
        df_max=df.values.max()
        df_min=df.values.min()
        row_n,col_n=df.shape

        """
        if not ratio and (col_n <= 5 or row_n <= 5): 
            ratio = 2
        else:
            ratio = 4
        """

        if not f_size_x: f_size_x = df.shape[1]
        if not f_size_y: f_size_y = df.shape[0]
        if f_size_x <= 6:
            f_size_x = f_size_x*2*ratio_w
        if f_size_y <= 6:
            f_size_y = f_size_y*2*ratio_h
        f_size_x = int(f_size_x * 1.5)
        f_size_y = int(f_size_y * 1.5)
        fig=plt.figure(figsize=(f_size_x, f_size_y))
        ax1 = plt.subplot2grid((ratio_h+2,ratio_w+1), (0,0), colspan=ratio_w,rowspan=1)
        ax2 = plt.subplot2grid((ratio_h+2,ratio_w+1), (1,0), colspan=ratio_w,rowspan=ratio_h)
        ax3 = plt.subplot2grid((ratio_h+2,ratio_w+1), (1,ratio_w), colspan=1,rowspan=ratio_h)
        ax4 = plt.subplot2grid((ratio_h+2,ratio_w+1), (ratio_h+1,0),colspan=ratio_w,rowspan=1)

        inset_ax = inset_axes(ax4,
                    width="100%", # width = 30% of parent_bbox
                    height="30%", # height : 1 inch
                    loc=1,
                    axes_kwargs={'axis_bgcolor':'white','frame_on':1,'visible':1,'xticks':[],'yticks':[]})
        ax4.axis('off')

        y_value=list(df.sum()) # column sum
        x_label=row_labels
        x_pos=range(len(x_label))
        #ax1.bar(x_pos,y_value,width=1,color='grey')
        #ax1.set_xticks(np.arange(df.shape[0])+0.5, minor=False)
        df.sum(axis=0).plot(kind='bar',stacked=True,legend=False,ax=ax1)
        x0,x1 = ax1.get_xlim()
        #ax1.set_xlim(x0+0.25,x1-0.25)
        if top_bar_x:
        #        ax1.set_xticklabels(x_label, minor=False)
                 ax1.get_xaxis().set_visible(False)

        #df.T.plot(kind='bar',stacked=True,ax=ax1,legend=False,color='grey')

        r=float(ratio)
        l=float(r/(r+1)/10)
        b=float(1/(3*r+6))
        w=float(3*r/(r+1)/5)
        h=float(1/(3*r+6))
        #print l,b,w,h
        #cbar_ax = fig.add_axes([l, b, w, h])

        if pct_by_row:
            df_pct = df.apply(lambda c: c / c.sum() * 100, axis=1)
        elif pct_by_col:
            df_pct = df.apply(lambda c: c / c.sum() * 100, axis=0)
        else:
            df_pct = df

        heatmap = sns.heatmap(df_pct.applymap(lambda x:int(x)),
                                linewidths=0.1, linecolor='white',annot=True,fmt="d",ax=ax2,
                                cbar_ax=inset_ax,cbar_kws={"orientation": "horizontal"},cmap="YlGnBu")

        x_value=list(df.sum(axis=1))[::-1] # row sum
        y=range(len(df.index))
        y_lable=column_labels[::-1]
        #ax3.barh(y,x_value,height=1,color="grey")
        #ax3.set_yticks(np.arange(df.shape[1])+0.5, minor=False)
        label_rev=df.index[::-1]
        df.sum(axis=1)[label_rev].plot(kind='barh',stacked=True,legend=False,ax=ax3)
        y0,y1=ax3.get_ylim()
        print y0,y1
        #ax3.set_ylim(y0+0.25,y1-0.25)
        if right_bar_y:
        #        ax3.set_yticklabels(y_lable, minor=False)
                 ax3.get_yaxis().set_visible(False)
        #ax3.set_xticklabels(rotation=45)
        for tick in ax3.get_xticklabels():
                tick.set_rotation(45)
        if title != '':
            plt.suptitle(title,fontsize=25)

        #fig.tight_layout(pad=0)
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        #fig.set_tight_layout(True)

        if savefn:
                fig.savefig(savefn)
        plt.close()

def ax_ticklabel_rot(ax,x_rot=90,y_rot=0):
    for tick in ax.get_xticklabels():
        tick.set_rotation(x_rot)
    for tick in ax.get_yticklabels():
        tick.set_rotation(y_rot)

def sns_heatmap_plus_annotate(df,ax2_x,ax2_y,ax2_count_val,ax1_x,ax1_y,ax3_x,ax3_y,savefn,ratio_h=2,ratio_w=2):
        fig=plt.figure(figsize=(20, 20))
        ax1 = plt.subplot2grid((ratio_h+2,ratio_w+1), (0,0), colspan=ratio_w,rowspan=1)
        ax2 = plt.subplot2grid((ratio_h+2,ratio_w+1), (1,0), colspan=ratio_w,rowspan=ratio_h)
        ax3 = plt.subplot2grid((ratio_h+2,ratio_w+1), (1,ratio_w), colspan=1,rowspan=ratio_h)
        ax4 = plt.subplot2grid((ratio_h+2,ratio_w+1), (ratio_h+1,0),colspan=ratio_w,rowspan=1)

        inset_ax = inset_axes(ax4,
                    width="100%", # width = 30% of parent_bbox
                    height="30%", # height : 1 inch
                    loc=1,
                    axes_kwargs={'axis_bgcolor':'white','frame_on':1,'visible':1,'xticks':[],'yticks':[]})
        ax4.axis('off')

        df_pivot_table = df.pivot_table(index=ax2_y,columns=ax2_x,values=ax2_count_val,aggfunc='count')

        heatmap = sns.heatmap(df_pivot_table.applymap(lambda x:int(x)),
                                linewidths=0.1, linecolor='white',annot=True,fmt="d",ax=ax2,
                                cbar_ax=inset_ax,cbar_kws={"orientation": "horizontal"},cmap="YlGnBu")

        # ax1
        sns.violinplot(x=ax1_x,y=ax1_y,data=df,ax=ax1)
 
        # ax2
        sns.violinplot(x=ax3_x,y=ax3_y,data=df,ax=ax3,orient='h')

        plt.savefig(savefn)
        plt.close()

### sns cluster map
def sns_clustermap(df=None,df_fn=None,col_as_color=None,cols_keep=None,cols_drop=None,colors=None,savefn=None,figsize_x=4,figsize_y=4):
    printFuncRun('sns_clustermap')
    printFuncArgs()
    if df_fn is None and df is None:
        print "[error] both df & df_fn are None, please provide either of them!"
        return
    if df_fn is None:
        print "df_fn is None, use df directly!"
        if savefn is None:
            print "savefn is None, please provide one for later saving"
            return
    if df is None:
        print "df is None, will read df from: %s"%(df_fn)
        df = pd.read_csv(df_fn,sep='\t',header=0)
        savefn = df_fn
    df_copy = df.copy()
    if not col_as_color is None:
        colors = df_copy[col_as_color]
        del df_copy[col_as_color]
    if not cols_keep is None:
        cols = df_copy.columns
        for col in cols_keep:
            if not col in cols:
                print "cols: %s not in columns"%(col),cols
                return
        df_copy = df_copy[cols_keep]
        print "only keep for clustering columns: ",cols_keep
    if not savefn.endswith(('.png','.pdf','.jpg')):
        savefn += '.clustermap.pdf'

    g = sns.clustermap(df_copy,figsize=(figsize_x,figsize_y),linewidths=0.1,row_colors=colors,col_cluster=False)
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0) # rotate horizental text label
    plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=45)
    g.savefig(savefn)
    plt.close()

    print "save fig to %s"%(savefn) 
    printFuncRun('sns_clustermap')


def scatter_corr(x, y, coor='pearson',**kws):
    if coor == 'pearson':
        r, p = stats.pearsonr(x, y)
    if coor == 'spearman':
        r, p = stats.spearmanr(x, y)
    ax = plt.gca()
    ax.scatter(x,y)
    ax.annotate("{} r = {:.2f}\np = {:.2e}".format(coor,r,p),
                xy=(.5, .5),ha='center',va='center',fontsize=20, xycoords=ax.transAxes)

def feature_count_gini_pct_df_or_fn_plot(pct_df=None,pct_fn=None,savefn=None,row_str=None,col_str=None,col_order=None,hue_str=None,map_x=None,map_y=None,type='violin',pval=0,show_xticks=1,show_xlabel=1,show_title=1,show_yticks=1,show_ylabel=1,size=5,aspect=1,col_wrap=None,sharex=True,sharey=True,p_val_loc=3,add_legend=False,x_rot=45,ytick_right=0,add_strip_df=None,add_strip_str=None,add_strip_top=95,margin_titles=False):
    printFuncRun('feature_count_gini_pct_df_or_fn_plot from zhangqf5')
    printFuncArgs()
    if pct_fn is None and pct_df is None:
        print "both pct_df & pct_fn are None, please provide either of them"
        return
    if pct_fn is None:
        print "pct_fn is None, will use pct_df"
        if savefn is None:
            print "savefn is None, please provide one"
            return
    if pct_df is None:
        print "pct_df is None, will read df from pct_fn: %s"%(pct_fn)
        pct_df = pd.read_csv(pct_fn,sep='\t',header=0)
        savefn=pct_fn

    pct_df_cols = pct_df.columns
    #use_cols = ['pct','feature','TE_type','gini']
    #for col in use_cols:
    #    if col not in pct_df_cols:
    #        print "%s not in pct_df_cols"%(col),pct_df_cols
    #        return
    print pct_df.info
    print pct_df.describe()

    if col_order is None:
        map_x_order = sorted(list(pct_df[map_x].value_counts().index))
    else:
        map_x_order = col_order
    #pct_df[col_str] = pct_df[col_str].astype(str)

    pct_df[map_x] = pd.Categorical(pct_df[map_x],categories=map_x_order,ordered=True)
    pct_df.sort_values(by=map_x,axis=0,inplace=True)


    #g=sns.FacetGrid(data=pct_df,row='pct',col='feature',size=5,aspect=1, margin_titles=False,hue='feature')
    print "start sns.FacetGrid ",timestamp()

    if col_wrap is not None:
        g=sns.FacetGrid(data=pct_df,row=row_str,col=col_str,size=size,aspect=aspect, margin_titles=margin_titles,hue=hue_str,col_wrap=col_wrap,sharex=sharex,sharey=sharey)
    else:
        g=sns.FacetGrid(data=pct_df,row=row_str,col=col_str,size=size,aspect=aspect, margin_titles=margin_titles,hue=hue_str,sharex=sharex,sharey=sharey)

    #g.map(sns.violinplot,'TE_type','gini',cut=0)
    #g.map(sns.violinplot,'TE_type','gini')
    if type == 'violin':
        g.map(sns.violinplot,map_x,map_y,order=map_x_order)
    if type == 'scatter':
        #g.map(sns.pointplot,map_x,map_y,linestyles='')
        g.map(plt.scatter,map_x,map_y,alpha=0.5)
        #g.map(sns.regplot,map_x,map_y,truncate=False)
        #g.map(sns.lmplot,map_x,map_y,hue_str)
    if type == "point":
        g.map(sns.pointplot,map_x,map_y,linestyles="",scale=0.1)
    if type == "bar":
        g.map(sns.barplot,map_x,map_y)
    if type == 'count':
        g.map(sns.countplot,map_x,order=map_x_order)
    if type == 'dist':
        #g.map(sns.distplot,map_x)
        g.map(sns.kdeplot,map_x)
    if type == 'strip':
        g.map(sns.stripplot,map_x,map_y,jitter=True,order=map_x_order)
    print "finish sns.FacetGrid ",timestamp()

    ax_arr = g.axes.flat

    if not show_xticks: g.set(xticks=[])
    if not show_yticks: g.set(yticks=[])
    if not show_xlabel: g.set(xlabel='')
    if not show_ylabel: g.set(ylabel='')
    if not show_title: g.set(title='')

    if add_legend:
        g = g.add_legend()

    if x_rot:
        for ax in ax_arr:
            plt.setp(ax.get_xticklabels(), rotation=x_rot)

    if ytick_right > 0:
        print "[detected]ytick_right(%s) >= 0"%(ytick_right)
        for ax in ax_arr:
            ax.yaxis.tick_right()

    #if add_strip_df is not None:
    #    sns.stripplot(x=map_x,y=map_y,data=add_strip_df,order=map_x_order,ax=ax_arr[0],jitter=True,color='blue',size=3)

    p_val_loc_dict = {1:(0.4,0.95),2:(0.1,0.95),3:(0.4,0.01),4:(0.1,0.01)}

    for n,ax in enumerate(ax_arr):
        ax_title=ax.get_title()
        if ax_title == "": continue
        pct = ax_title.split('|')[0].split('=')[1].strip()
        row_str_val = ax_title.split('|')[0].split('=')[1].strip()
        feature = ax_title.split('|')[1].split('=')[1].strip()
        col_str_val = ax_title.split('|')[1].split('=')[1].strip()
        df_dtypes = pct_df.dtypes
        pct_dtype = df_dtypes[row_str]
        feature_dtype = df_dtypes[col_str]
        try:
             pct = float(pct)
        except:
             pass
        print 'annotate axes: %s'%(n),ax,'ax_title: %s'%(ax_title),'feature:%s pct:%s'%(feature,pct)
        try:
            df_plot = pct_df[(pct_df[col_str]==feature) & (pct_df[row_str]==pct)]
        except:
            print traceback.print_exc()
            print "col_str: %s, feature: %s ; row_str: %s, pct: %s"%(col_str,feature,row_str,pct)
            df_plot = pct_df[(pct_df[col_str]==feature) & (pct_df[row_str]==float(pct))]
        #print df_plot
        if pval:
            uniq_val_ls = list(df_plot[map_x].value_counts().index)
            ha = 'left' if p_val_loc == 1 or p_val_loc == 3 else 'right'
            va = 'top' if p_val_loc == 1 or p_val_loc == 2 else 'bottom'
            if len(uniq_val_ls) == 0:
                annotate_str='n=%s'%(df_plot.shape[0])
                ax.annotate(annotate_str,xy=(.01, .01),ha='center',va='bottom',fontsize=10, xycoords=ax.transAxes)
            else:
                #x = df_plot[df_plot['TE_type']=='slow']['gini']
                #y = df_plot[df_plot['TE_type']=='fast']['gini']
                p_val_str_ls = []
                for map_x_val1,map_x_val2 in itertools.combinations(map_x_order,2):
                    map_x_val1_ls = df_plot[df_plot[map_x]==map_x_val1][map_y]
                    map_x_val2_ls = df_plot[df_plot[map_x]==map_x_val2][map_y]
                    p = scipy.stats.ks_2samp(map_x_val1_ls,map_x_val2_ls)[1]
                    p_val_str_ls.append("%s vs %s: p(ks)=%.2e"%(map_x_val1,map_x_val2,p))
                #p=scipy.stats.ttest_ind(x,y)[1]
                #p=scipy.stats.ks_2samp(x,y)[1]
                #num = str(len(x)) if len(x) == len(y) else str(len(x))+'|'+str(len(y))
                #print "slow vs fast: pval ",p
                #annotate_str='n=%s,p=%.2e'%(df_plot.shape[0],p)
                #ax.annotate('n=%s,p=%.2e'%(num,p),xy=(.5, .01),ha='center',va='bottom',fontsize=20, xycoords=ax.transAxes)
                ax.annotate('\n'.join(p_val_str_ls),xy=p_val_loc_dict[p_val_loc],ha=ha,va=va,fontsize=8, xycoords=ax.transAxes)
        if add_strip_str is not None:
              add_strip_str_top = np.percentile(df_plot[df_plot[map_x]==add_strip_str][map_y],add_strip_top)
              df_add_strip_str = df_plot[(df_plot[map_x]==add_strip_str) & (df_plot[map_y]>=add_strip_str_top)]
              df_add_strip_str_id_ls = list(df_add_strip_str['id'])
              df_add_strip = df_plot[df_plot['id'].isin(df_add_strip_str_id_ls)]
              sns.stripplot(x=map_x,y=map_y,data=df_add_strip,order=map_x_order,ax=ax,jitter=True,color='blue',size=3)
        if add_strip_df is not None:
              df_add_strip = add_strip_df[(add_strip_df[col_str]==feature) & (add_strip_df[row_str]==pct)]
              sns.stripplot(x=map_x,y=map_y,data=df_add_strip,order=map_x_order,ax=ax,jitter=True,color='blue',size=3)



    plt.close()
    if not savefn.endswith(('png','.pdf','.jpg')):
        savefn += '.'+type+'.pdf'
    g.savefig(savefn)

    print "save fig to: %s"%(savefn)
    printFuncRun('feature_count_gini_pct_df_or_fn_plot')
    return pct_df,savefn

def tableview(txt=None, parse_table=None, make_conf=None, conf_dir=None, circos_conf=None, save_dir=None, parsed_conf=None):
        if parse_table is None:
                parse_table = '/Share/home/zhangqf5/gongjing/software/circos-tools-0.22/tools/tableviewer/bin/parse-table'
        if make_conf is None:
                make_conf = '/Share/home/zhangqf5/gongjing/software/circos-tools-0.22/tools/tableviewer/bin/make-conf'
        if conf_dir is None:
                conf_dir = '/Share/home/zhangqf5/gongjing/software/circos-tools-0.22/tools/tableviewer/data'
        if circos_conf is None:
                circos_conf = '/Share/home/zhangqf5/gongjing/software/circos-tools-0.22/tools/tableviewer/etc/circos.conf'
        if save_dir is None:
                save_dir = os.path.dirname(txt)
        if parsed_conf is None:
                parsed_conf = '/Share/home/zhangqf5/gongjing/software/circos-tools-0.22/tools/tableviewer/samples/parse-table-02a.conf'
        print "[tableview start] file: %s"%(txt)
        printFuncArgs()
        subprocess.call(["cat {txt} | {parse_table} -conf {parsed_conf} | {make_conf} -dir {conf_dir}".format(txt=txt, parse_table=parse_table, parsed_conf=parsed_conf, make_conf=make_conf, conf_dir=conf_dir)],shell=True)
        file_png = txt.split('/')[-1].replace('txt', 'png')
        subprocess.call(["circos -conf {circos_conf} -outputdir {save_dir} -outputfile {file_png} | grep created".format(circos_conf=circos_conf, save_dir=save_dir, file_png=file_png)],shell=True)
        print "[tableview end] file: %s"%(txt)

def tableview_call(txt):
    tableview_dir = '/Share/home/zhangqf5/gongjing/software/circos-tools-0.22/tools/tableviewer2'
    tableview_scirpt = '/Share/home/zhangqf5/gongjing/software/circos-tools-0.22/tools/tableviewer2/makeimage.py'
    subprocess.call(["cd %s; python %s %s; cd -"%(tableview_dir, tableview_scirpt, txt)],shell=True)


def plot_color_ls(rgb_ls=None, label_ls=None, savefn=None, label_font_size=30, label_linewidth=28):
    if rgb_ls is None:
        rgb_ls = ['202,75,78', '83,169,102', '205,185,111', '98,180,208', '129,112,182', '255,0,0', '0,128,0', '74,113,178', '169,169,169']
    if label_ls is None:
        label_ls = ['proteinCoding', 'lncRNA', 'miRNA', 'rRNA', 'snoRNA', 'snRNA', 'tRNA', 'NoncanonicalRNA', 'others']
    rgb_ls = [map(lambda x:x/float(255), map(int, i.split(','))) for i in rgb_ls]
    
    with sns.axes_style('white'):
        fig, ax = plt.subplots(figsize=(10,10))
    for n,(i,j) in enumerate(zip(label_ls,rgb_ls)):
        x = [1,2,3]
        y = [n,0,0]
        ax.plot(x, y, color=j, label=" "*2+i)
    legend = ax.legend()
    for n,label in enumerate(legend.get_texts()):
        label.set_fontsize(label_font_size)
        label.set_color(rgb_ls[n])
    for label in legend.get_lines():
        label.set_linewidth(label_linewidth)
    if savefn is not None:
        plt.savefig(savefn)
    plt.close()

##################################################
### dict related
# 
def read_file_2_dict(input_f="/Share/home/zhangqf/gongjing/paris-2016-05/result/SRR2814700/Aligned.out.native.cache_sub/mutation.pos.txt",line_split_sep='\t',key_col=[0],key_join=' ',val_col=[0],val_join=''):
    print "[read_file_2_dict]",timestamp()
    print "  - input_f: %s"%(input_f)
    print "  - key_col:",key_col
    print "  - val_col:",val_col
    file_dict={} # key: key_col value: count
    with open(input_f,'r') as IN:
        for line in IN:
            line=line.rstrip()
            if line and not line.startswith("#"):
                arr=line.split(line_split_sep)
                key=key_join.join([arr[i].strip() for i in key_col])
                #key=arr[key_col].rstrip()
                val =  val_join.join([arr[i] for i in val_col])
            if not file_dict.has_key(key):  # keep the first value if a key appears multile times
                file_dict[key] = val
    key0 = file_dict.keys()[0]
    print "  - dict len: %s,   key:val -> %s:%s"%(len(file_dict),key0,file_dict[key0])
    print "[read_file_2_dict]",timestamp()
    print
    return file_dict

def read_txt_2_dict(input_f=None, line_split_sep='\t', header=0, key_col=[0]):
	print "[read_txt_2_dict]",timestamp()
	d = nested_dict()
	with open(input_f, 'r') as IN:
		for n,line in enumerate(IN):
			line = line.strip()
			if n == 0 and header == 0:
				header = line.split(line_split_sep)
				continue
			if not line or line.startswith('#'): continue
			arr = line.split(line_split_sep)
			key=' '.join([arr[i].strip() for i in key_col])
			for i,j in zip(header, arr):
				d[key][i] = j
	key0 = d.keys()[0]
	print "  - dict len: %s,   key:val -> %s:%s"%(len(d),key0,d[key0])
	print "[read_txt_2_dict]",timestamp()
	return d.to_dict()

# print a dict key,val
# used for print fucntion args by passing locals() to d
def printLocal(d):
    for key,val in d.items():
        print "  - %s: %s"%(key,val)

# print func name & run timestamp
def printFuncRun(name,loc=0):
    if loc == 0:
        print "[%s]"%(name),timestamp()
    elif loc == 1:
        print "[%s]"%(name),timestamp()+'\n'
    else:
        pass

# print function args
def printFuncArgs():
        """Returns tuple containing dictionary of calling function's
           named arguments and a list of calling function's unnamed
           positional arguments.
        """
        #from inspect import getargvalues, stack
        posname, kwname, args = getargvalues(stack()[1][0])[-3:]
        posargs = args.pop(posname, [])
        args.update(args.pop(kwname, []))
        if len(args)>0:
            printLocal(args)
        if len(posargs)>0:
            print posargs
        return args, posargs

def print_dict_item0(mydict,num=1):
    key0_ls = mydict.keys()[0:num]
    for key0 in key0_ls:
        print "key0: val0 -> ",key0,mydict[key0]
    return key0_ls[0]

def print_dict(d):
    print json.dumps(d, indent=4)

##################################################
### list related
def sort_str_num_ls(ls=[1,2,3]):
    if isinstance(ls[0],int):
        return sorted(ls)
    if isinstance(ls[0],str):
        try:
            return map(str,sorted([int(i) for i in ls]))
        except:
            return sorted(ls)

def find_all_value_index_in_list(lst=[1,2,3,4,5,1],f=1):
    return [i for i, x in enumerate(lst) if x == f]

def list_list_sum(lists=[[1,2],[3,4]],mode='count_sum'):
    if mode == 'count_sum':
        total=sum(sum(ls) for ls in lists)
    if mode == "len_sum":
        total=sum(len(ls) for ls in lists)
    return total

def ls_ls_flat(ls_ls):
    return list(itertools.chain.from_iterable(ls_ls))

# get a number of elements from a list by index
# i: 0-based(python way index)
def getElementsFromList(ls,i): 
    ls_sub=map(ls.__getitem__,i)
    return ls_sub

# list to percent list
def list_pct(ls):
    ls=map(float,ls)
    ls_sum=sum(ls)
    ls_pct=[i/ls_sum for i in ls]
    return ls_pct

# refine a list: 
# input: [1,2,2,6]
# output: [1,2,2.001,6]
def refine_list(ls):
    ls_refine=[]
    for n,i in enumerate(ls):
        if n == 0:
            ls_refine.append(i)
        else:
            if i == ls[n-1]:
                i=ls_refine[n-1]+1/float(1000)
            ls_refine.append(i)
    return ls_refine

# remove non-digital elements from list
def list_rm_non_digital(ls):
    new_ls=[]
    for i in ls:
        try:
            i=float(i)
            status=1
        except:
            status=0
        if status: new_ls.append(i)
    return new_ls

# sort str list by separator
def str_num_ls_sorted(ls,sep,str=0):
    if str:
        return sep.join(sorted(ls,key=lambda x:(float(x.split(sep)[0]),float(x.split(sep)[1]))))
    else: 
        return sorted(ls,key=lambda x:(float(x.split(sep)[0]),float(x.split(sep)[1])))

# plot list of sequence as 'logo'
# ['atc','agc']
def ls_weblogo_plot(ls_ls,pos_labels=None,ax=None,title=''):
    ls_ls=[list(i) for i in ls_ls]
    pos_len=len(ls_ls[0])
    pos_df=pd.DataFrame(ls_ls)
    pos_df_t=pos_df.apply(lambda x:x.value_counts()).T.fillna(0)
    if pos_labels != None: pos_df_t.index=pos_labels
    if ax == None: ax=plt.gca()
    pos_df_t.plot(kind='bar',stacked=True,ax=ax,title=title+' '+str(len(ls_ls)))  

# convert list to .fa
def list2fa(seq_ls,savefn):
    if type(seq_ls) != list:
        print "[error] provided seq_ls are not List type",seq_ls
    if not savefn.endswith('fa'):
        print "[error] provided savefn Not .fa: %s"%(savefn)
        sys.exit()
    savefn=add_file_pwd(savefn)
    printFuncRun('list2fa')
    #printFuncArgs()
    print "  - savefn: %s"%(savefn)
    FA=open(savefn,'w')
    for n,seq in enumerate(seq_ls):
        print >>FA,">"+str(n)
        print >>FA,seq.strip()
    FA.close()
    printFuncRun('list2fa')
    return savefn

# describe a list of list, return the describe df
def list_list_describe(ls_ls,ls_labels,round=0,sample_col_string='sample',val_col_string='value'):
    ls_ls_describe=[]
    ls_ls_df_ls=[]
    for ls,ls_label in zip(ls_ls,ls_labels):
        ss=pd.Series(ls)
        ss_describe=ss.describe()
        ls_ls_describe.append(ss_describe)
        
        ls_df = pd.DataFrame({val_col_string:ls,sample_col_string:ls_label})
        ls_ls_df_ls.append(ls_df)
        
    df_describe= pd.concat(ls_ls_describe,axis=1)
    df_describe.columns = ls_labels
    df_describe=df_describe.round(round)
    
    print df_describe
    
    ls_ls_df = pd.concat(ls_ls_df_ls,axis=0)
    
    return df_describe,ls_ls_df

# sns plot box+strip
def list_list_box_strip_plot(ls_ls,ls_labels,sample_col_string='sample',val_col_string='value',fig=None,ax=None,savefn=None,title=None,violin_box='box',strip=True):
    df_describe,ls_ls_df = list_list_describe(ls_ls,ls_labels,round=0,sample_col_string=sample_col_string,val_col_string=val_col_string)
    
    ax1=ax
    if ax1 == None:
        fig_y = len(ls_ls) if len(ls_ls) > 8 else 6
        fig,ax1=plt.subplots(figsize=(8,fig_y))
    else:
        fig=plt.gcf()

    #g = sns.boxplot(x=val_col_string, y=sample_col_string, hue=sample_col_string,data=ls_ls_df,whis=np.inf,ax=ax1) # each sample -> different color, cover max/min
    if violin_box == 'box':
        g = sns.boxplot(x=val_col_string, y=sample_col_string, order=ls_labels,data=ls_ls_df,whis=np.inf,ax=ax1) # each sample -> different color, cover max/min
    if violin_box == 'violin':
        g = sns.violinplot(x=val_col_string, y=sample_col_string, order=ls_labels,data=ls_ls_df,ax=ax1) # each sample -> different color, cover max/min
    if strip:
        g = sns.stripplot(x=val_col_string, y=sample_col_string, data=ls_ls_df,jitter=True, size=3, color=".3", linewidth=0,ax=ax1)
 
    if savefn:
        try:
            plt.tight_layout()
        except:
            pass
        fig.savefig(savefn)
        plt.close()
    else:
        return g


# plot box from list of list values
def list_list_box_plot(list_list,list_label,fig=None,ax=None,savefn=None,title=None,xlabel=None,ylabel=None):
        ax1=ax

        if ax1 == None:
            fig_x = len(list_list) if len(list_list) > 6 else 5 
            fig,ax1=plt.subplots(figsize=(fig_x,6))
        else:
            fig=plt.gcf()

        data=list_list
	data_tickname=list_label

	fig.canvas.set_window_title('A Boxplot Example')
	plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

	bp = plt.boxplot(data, notch=0, sym='+', vert=1, whis=1.5)
	plt.setp(bp['boxes'], color='black')
	plt.setp(bp['whiskers'], color='black')
	plt.setp(bp['fliers'], color='red', marker='+')

	# Add a horizontal grid to the plot, but make it very light in color
	# so we can use it for reading data values but not be distracting
	ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',alpha=0.5)

	# Hide these grid behind plot objects
	ax1.set_axisbelow(True)
        title = title if title else "bar plot"
        xlabel = xlabel if xlabel else 'distribution'
        ylabel = ylabel if ylabel else 'value'
	ax1.set_title(title)
	ax1.set_xlabel(xlabel)
	ax1.set_ylabel(ylabel)

	numDists=len(list_list)

	# Now fill the boxes with desired colors
	boxColors = ['darkkhaki', 'royalblue','darkred','darkgreen','darkgrey','lightpink','lightgreen','lightsteelblue','c','orange','lightcoral','chocolate','mediumslateblue']
	numBoxes = numDists
	medians = list(range(numBoxes))
	for i in range(numDists):
	    box = bp['boxes'][i]
	    boxX = []
	    boxY = []
	    for j in range(5):
		boxX.append(box.get_xdata()[j])
		boxY.append(box.get_ydata()[j])
	    boxCoords = list(zip(boxX, boxY))
	    # Alternate between Dark Khaki and Royal Blue
	    k = i % len(boxColors)
	    boxPolygon = Polygon(boxCoords, facecolor=boxColors[k])
	    ax1.add_patch(boxPolygon)
	    # Now draw the median lines back over what we just filled in
	    med = bp['medians'][i]
	    medianX = []
	    medianY = []
	    for j in range(2):
		medianX.append(med.get_xdata()[j])
		medianY.append(med.get_ydata()[j])
		plt.plot(medianX, medianY, 'k')
		medians[i] = medianY[0]
	    # Finally, overplot the sample averages, with horizontal alignment
	    # in the center of each box
	    plt.plot([np.average(med.get_xdata())], [np.average(data[i])],
		     color='w', marker='*', markeredgecolor='k')

	# Set the axes ranges and axes labels
	ax1.set_xlim(0.5, numBoxes + 0.5)
	top = max([j for i in data for j in i]) * 1.3
	bottom = 0 if min([j for i in data for j in i]) > 0 else min([j for i in data for j in i])*1.1
	ax1.set_ylim(bottom, top)
	xtickNames = plt.setp(ax1, xticklabels=data_tickname)
        #xtickcolors=[i%len(boxColors) for i,j in enumerate(data_tickname)]
        #for i,t in enumerate(ax1.get_xticklabels()): t.set_color(xtickcolors[i])
	plt.setp(xtickNames, rotation=45, fontsize=16)

	# Due to the Y-axis scale being different across samples, it can be
	# hard to compare differences in medians across the samples. Add upper
	# X-axis tick labels with the sample medians to aid in comparison
	# (just use two decimal places of precision)
	pos = np.arange(numBoxes) + 1
	upperLabels = [str(np.round(s, 2))+'\n'+'('+str(l)+')' for s,l in zip(medians,[len(i) for i in data])]
	weights = ['bold', 'semibold']
	for tick, label in zip(range(numBoxes), ax1.get_xticklabels()):
	    k = tick % len(boxColors)
	    #k=tick
	    ax1.text(pos[tick], top - (top*0.05), upperLabels[tick],
		     horizontalalignment='center', va='center',size='x-small', weight='bold',
		     color=boxColors[k])

	# Finally, add a basic legend
	#plt.figtext(0.80, 0.08, ' Random Numbers',backgroundcolor=boxColors[0], color='black', weight='roman',size='x-small')

        if savefn:
            try:
                 plt.tight_layout()
            except:
                 pass
	    plt.close()
	    fig.savefig(savefn)


def ls_ls_common(ls_ls,return_ls=1):
    if return_ls:
        return list(set(ls_ls[0]).intersection(*ls_ls))

def ls_ls_union(ls_ls,return_ls=1):
    if return_ls:
	    return list(set(ls_ls[0]).union(*ls_ls))

def ls_ls_uniq(ls_ls):
    """ [[1,2,3],[1,2],[1,2]] -> [[1,2],[1,2,3]] """
    return [list(x) for x in set(tuple(x) for x in ls_ls)]

def ls_ls_overlap(ls_ls,ls_ls_labels,savefn=None,title_str=None,figsize_x=8,figsize_y=8,keep_diagonal=1):
    ol_dict =defaultdict(dict)
    for n,(ls,ls_label) in enumerate(zip(ls_ls,ls_ls_labels)):
        #other_ls_ls = [(i,j) for i,j in zip(ls_ls,ls_ls_labels) if j != ls_label]
        other_ls_ls = [(i,j) for i,j in zip(ls_ls,ls_ls_labels) if j != '']
        for other_ls,other_ls_label in other_ls_ls:
            ol = set(ls) & set(other_ls)       
            print ls_label,other_ls_label,len(ol)
            ol_dict[ls_label][other_ls_label] = len(ol)
    ol_df = pd.DataFrame(ol_dict)
    ol_df.replace(np.NaN,0,inplace=True)
    ol_df = ol_df.loc[ls_ls_labels,ls_ls_labels]
    print ol_df
    if savefn is not None:
        mask = np.zeros_like(ol_df)
        mask[np.triu_indices_from(mask)] = True
        if keep_diagonal:
            np.fill_diagonal(mask,0)  # keep diagonal
        with sns.axes_style("white"): # seaborn set temp axes style firdt, then create subplots
            fig,ax = plt.subplots(figsize=(figsize_x,figsize_y))
            sns.heatmap(ol_df,annot=True,fmt='.0f',cbar_kws={"orientation": "vertical"},linewidths=.5,square=True,mask=mask,cbar=False)
        if title_str is not None: plt.title(title_str)
        plt.tight_layout()
        plt.savefig(savefn)
        plt.close()
        print "save heatmap to %s"%(savefn)
    return ol_df



def list_smooth_gaussian(ls,degree=5):
    ls = np.pad(ls,(degree-1,degree-1),mode='edge')
    window = degree*2-1
    weight = np.arange(-degree+1,degree)/window
    weight = np.exp(-16*weight**2)
    weight /= sum(weight)
    smoothed = np.convolve(ls,weight,mode='valid')
    return smoothed

def list_smooth_Triangle(ls,degree=5,dropVals=False):
	"""performs moving triangle smoothing with a variable degree."""
	"""note that if dropVals is False, output length will be identical
	to input length, but with copies of data at the flanking regions"""
	triangle=numpy.array(range(degree)+[degree]+range(degree)[::-1])+1
	smoothed=[]
	for i in range(degree,len(data)-degree*2):
		point=data[i:i+len(triangle)]*triangle
		smoothed.append(sum(point)/sum(triangle))
	if dropVals: return smoothed
	smoothed=[smoothed[0]]*(degree+degree/2)+smoothed
	while len(smoothed)<len(data):smoothed.append(smoothed[-1])
	return smoothed

##################################################
### string related
# remove from end of string by occurrence, eg: replace filename subfix
def rreplace(s="1232425", old='2', new=' ', occurrence=1):
    li=s.rsplit(old,occurrence)
    return new.join(li)

# assign empty string a new value, default None
def str_empty_assign(string,assign=None):
    string = assign if string == "" else string
    return string

##################################################
### json,pickle,cpickle related, for later fast loading, especially for frequently used objects(str,list,dict et.)
def dump_json(out_json,mydict):    
    output = open(out_json, 'wb')
    json.dump(mydict, output)
    output.close()    

def dump_pickle(out_pickle,mydict):    
    output = open(out_pickle, 'wb')
    pickle.dump(mydict, output,protocol=cPickle.HIGHEST_PROTOCOL)
    output.close()

def dump_cpickle(out_cpickle,mydict):    
    output = open(out_cpickle, 'wb')
    cPickle.dump(mydict, output,protocol=cPickle.HIGHEST_PROTOCOL)
    output.close()

def load_json(in_json):
    input = open(in_json, 'rb')
    mydict = json.load(input)
    input.close()
    return mydict

def load_pickle(in_pickle):
    input = open(in_pickle, 'rb')
    mydict = pickle.load(input)
    input.close()
    return mydict

def load_cpickle(in_cpickle):
    input = open(in_pickle, 'rb')
    mydict = pickle.load(input)
    input.close()
    return mydict


##################################################
### pandas dataframe related
def dataFrameGroupBy(df,group_ls):
    t=df.groupby(group_ls).size()
    td=pd.DataFrame(t)
    return td

# label dataframe column by percentile
def df_label_col(df,col_str,per_n):
    bin=100/float(per_n)
    try:
        df[col_str]=df[col_str].astype(float)
    except:
        print "%s: not all can convert to float"%(col_str)
    print "label %s by %s"%(col_str,per_n)
    per_ls=[bin*i for i in xrange(per_n+1)]
    step_ls=[np.percentile(df[col_str],i) for i in per_ls]
    step_ls=refine_list(step_ls)
    print "step_ls: ",step_ls
    col_name=col_str+'_label'
    df[col_name]=df.apply(lambda row:col_label(row,col_str,step_ls,per_ls)[1],axis=1)
    print "finish label %s by %s"%(col_str,per_n)
    
def col_label(row,col_str,step_ls,per_ls):
    val=float(row[col_str])
    if np.isnan(val):
        step="nan"
    if not type(val) is float:
        return "."
    label=""
    val_label=""
    demical = 1 if val > 1 else 3
    # "{1:.{0}f}".format(demical, 20000 + 1/float(3)) # specify demical precision num(nest)
    for i,j in enumerate(step_ls):
        if i == 0: continue
        if step_ls[i-1] <= val < j: 
            label=str(per_ls[i-1])+'-'+str(per_ls[i])
            #val_label="%.3f"%(step_ls[i-1])+'-'+'%.3f'%(step_ls[i])
            val_label="{1:.{0}f}".format(demical,step_ls[i-1])+'-'+"{1:.{0}f}".format(demical,step_ls[i])
            break
    if label == "" and val == step_ls[-1]: # as to biggest elements, add the label
        label=str(per_ls[-2])+'-'+str(per_ls[-1])
        #val_label="%.3f"%(step_ls[-2])+'-'+'%.3f'%(step_ls[-1])
        val_label="{1:.{0}f}".format(demical,step_ls[-2])+'-'+"{1:.{0}f}".format(demical,step_ls[-1])
    return label,val_label
    
# save df as a table picture
def df_to_table_png(df,fig=None,ax=None,colLabels=None,rowLabels=None,cellLoc='center',colwidth=0.17,savefn=None,scale_x=1.2,scale_y=1.2,fontsize=14):
    if fig == None and ax == None:
        ax=plt.subplot(111,frame_on=False)
    if not rowLabels:
        rowLabels=['']*df.shape[0] 
    if rowLabels == 1:
        rowLabels = None
    #ax = plt.subplot(111, frame_on=False) # no visible frame
    colWidths=[colwidth]*df.shape[1] if colwidth else None 
    printFuncRun('df_to_table_png')
    printFuncArgs()
    tb=table(ax,df,loc='center',cellLoc=cellLoc,rowLabels=rowLabels,colLabels=colLabels,colWidths=colWidths)
    ax.axis('tight')
    ax.xaxis.set_visible(False)  # hide the x axis
    ax.yaxis.set_visible(False)  # hide the y axis
    for key, cell in tb.get_celld().items():
        if key[0] == 0: # column name
            #cell.set_facecolor('red') # fill cell with red
            cell.set_linewidth(0)
            cell.set_text_props(color='red',weight='bold')
        if key[0] > 0:  # value rows
            cell.set_linewidth(0)
    tb.auto_set_font_size(False)
    tb.set_fontsize(fontsize)
    tb.scale(scale_x,scale_y)
    if savefn:
        plt.gcf().savefig(savefn)
    printFuncRun('df_to_table_png')

# for a df, set column >= 76 and index >=76 as one entery(combine columns/indexs)
def df_part_add_reformat(df,index_loc,column_loc):
    print df
    
    index_ix=list(df.index).index(index_loc)
    column_ix=list(df.columns).index(column_loc)
    
    df_1=df.iloc[:index_ix,column_ix:]
    df_2=df.iloc[index_ix:,:column_ix]
    df_3=df.iloc[index_ix:,column_ix:]
    
    df_1_sum=df_1.apply(sum,axis=1) # row sum
    df_2_sum=df_2.apply(sum,axis=0) # column sum
    df_3_sum=df_3.apply(sum,axis=0).sum() # all sum
    
    df_new_top=df.iloc[:index_ix,:column_ix]
    df_new_top['>='+str(column_loc)] = df_1_sum
 
    df_new_bottom=pd.DataFrame({'>='+str(index_loc):df_2_sum}).T
    df_new=pd.concat([df_new_top,df_new_bottom],axis=0)
    
    df_new.ix[-1,-1]=df_3_sum
    
    print df_new
    return df_new

# get data frame largest & sub-largest value for truncated bar plot
def get_df_trucated_val(df,fold=10,ratio=0.1):
    val_ls=np.sort(df.values,axis=None)[::-1]
    print val_ls
    
    status = 0

    for n,val in enumerate(val_ls[0:-1]):
        if val >= fold*val_ls[n+1]:
	    status = 1
            #print "one huge gap: %s[%s] <-> %s[%s]"%(val,n,val_ls[n+1],n+1)
            return status,val_ls[n+1]*(1+ratio),val*(1-ratio)

    if status == 0:
        return status,0,0


###  A B C
###  1 2 3
###  4 5 6

### A 1
### A 4
### B 2
### B 5
### C 3
### C 6
def df_melt(df,col_str='sample',val_str='val'):
    cols = df.columns
    dt = df.unstack(fill_value=cols,level=-1)
    ll = list(dt.index.get_level_values(level=0))  # A A B B C C
    return pd.DataFrame({col_str:ll,val_str:dt.values})

def df_to_md(df, savefn=None):
    printFuncRun('df_to_md')
    printFuncArgs()
    cols_keep = ['index']+list(df.columns)
    df['index'] = df.index
    cols = df.columns
    df2 = pd.DataFrame([['---']*len(cols)], columns=cols)
    df3 = pd.concat([df2,df])
    df3.replace('_','\_',inplace=True)
    df3[cols_keep].to_csv(savefn, sep='|', index=False)
    printFuncRun('df_to_md')

##################################################
### sam file realted
def check_sam_cigar_pattern(input_f="/Share/home/zhangqf/gongjing/paris-2016-05/result/SRR2814766_NM4/Aligned.out.sam",eg=0):
    print "check cigar pattern of %s"%(input_f)
    cigar_pattern_count_dict=collections.defaultdict(int)
    cigar_pattern_eg_dict={}
    with open(input_f,'r') as IN:
        for n,line in enumerate(IN):
            line=line.rstrip()
            if line and not line.startswith("@"):
                if n%1000000 == 0: print "process line %d"%(n)
                arr=line.split()
                #cigarstring=re.sub("[0-9]","",arr[5]) # need compile, slow
                cigarstring=arr[5].translate(None, digits) # 2'48
                #cigarstring=filter(lambda x: x.isalpha(), arr[5]) # 2'6
                MD_tag=arr[16]
                #print arr[3],cigarstring
                cigar_pattern_count_dict[cigarstring] += 1
                if eg:
			if not cigar_pattern_eg_dict.has_key(cigarstring):
			    cigar_pattern_eg_dict[cigarstring] = [len(MD_tag),line]
			elif len(MD_tag) > cigar_pattern_eg_dict[cigarstring][0]:
			    cigar_pattern_eg_dict[cigarstring] = [len(MD_tag),line]
			else:
			    pass
    #print cigar_feature_ls
    #cigar_feature_count_dict=collections.Counter(cigar_feature_ls)
    print cigar_pattern_count_dict
    if eg:
	    #with open(input_f.replace('sam','cigarPattern.sam'),'w') as OUT:
	    with open(rreplace(input_f,'sam','cigarPattern.sam'),'w') as OUT:
		for key,val in cigar_pattern_eg_dict.items():
		    print >>OUT,val[1]
	    print cigar_pattern_eg_dict

    df=pd.DataFrame.from_dict(cigar_pattern_count_dict,orient="index")
    print df
    df.columns=['count']
    ls=[]
    row_num=[]
    N_max=max([i.count('N') for i in df.index])
    for j in xrange(N_max+1):   
        for n,i in enumerate(df.index):
        #print i,i.count('N'),j
            if i.count('N') == j:
                ls.append(i)
                row_num.append(n)  
    df2=df.iloc[row_num,:]
    y_max=df2['count'].max()
    fig=plt.figure()
    df2.plot(kind='barh',figsize=(16,14),legend=False,color=['g'],fontsize=14,width=0.5,position=0.5)
    for i, label in enumerate(list(df2.index)):
        score = df2.ix[label]['count']
        y_pos=score*1.2
        plt.gca().annotate(str(score), (1, i),color='b')
    plt.xlabel("count")
    plt.ylabel("cigar string type")
    plt.title(input_f.split("/")[-1]+'\n'+'mapped cigar string pattern')
    #plt.savefig(input_f.replace('sam','cigarPattern.sam.bar.png'))
    plt.savefig(rreplace(input_f,'sam','cigarPattern.sam.bar.png'))
 
def sam_get_uniq_mapping(input_sam="/Share/home/zhangqf/gongjing/paris-2016-05/result/SRR2814766/Aligned.out.sam"):
    print "[sam_get_uniq_mapping]",timestamp()
    sam_keep_alignment_first_n_hits(input_sam=input_sam,HIn=1,NHn=1)
    print "[sam_get_uniq_mapping]",timestamp()

def sam_keep_alignment_first_n_hits(input_sam="/Share/home/zhangqf/gongjing/paris-2016-05/result/SRR2814766/Aligned.out.sam",HIn=1,NHn=100,writeout=1):
    assert input_sam.endswith("sam"), "[sam error] input_sam should be sam file: %s"%(input_sam)
    #out_sam=input_sam.replace("sam","keep"+str(HIn)+'.sam')
    out_sam=rreplace(input_sam,"sam","keep"+str(HIn)+'.sam')
    OUT=open(out_sam,'w')
    print "[sam_keep_alignment_first_n_hits]",timestamp()
    print "  - input_sam: %s"%(input_sam)
    print "  - out_sam: %s"%(out_sam)
    print "  - NH max: %s"%(NHn)
    print "  - HI max: %s"%(HIn)
    all_map_line=0
    uniq_reads_num=0
    gapped_reads=[]
    multi_reads_num=collections.defaultdict(int)
    with open(input_sam,'r') as IN:
        for n,line in enumerate(IN):
            if n%1000000 == 0: print "  * process line: %s"%(n)
            if line.startswith("@") or line.startswith("#"): 
                print >>OUT,line.strip()
            else:
                all_map_line += 1
                arr=line.split("\t")
                NH,HI=map(lambda x:int(x.split(":")[-1]),arr[11:13])
                cigar=arr[5]
                multi_reads_num[NH] += 1
		if "N" in cigar: gapped_reads.append(arr[0])
                if writeout and NH <= NHn and HI <= HIn:
                    print >>OUT,line.strip()
    OUT.close()
    print "[sam_keep_alignment_first_n_hits]",timestamp()
    print
    uniq_reads_num=multi_reads_num[1]
    multi_reads_num_refine={i:j/i for i,j in multi_reads_num.items() if i>=2}
    return out_sam,all_map_line,uniq_reads_num,sum(multi_reads_num_refine.values()),len(set(gapped_reads))

def sam2bam(input_sam="/Share/home/zhangqf/gongjing/paris-2016-05/result/SRR2814766/Aligned.out.sam",ref_fa=""):
    #out_bam=input_sam.replace("sam",'bam')
    out_bam=rreplace(input_sam,"sam",'bam')
    print "[sam2bam]",timestamp()
    print "  - input_sam: %s"%(input_sam)
    print "  - out_bam: %s"%(out_bam)
    subprocess.call(["samtools view -b %s > %s"%(input_sam,out_bam)],shell=True)
    print "[sam2bam]",timestamp()
    print
    return out_bam

def sortbam(input_bam="/Share/home/zhangqf/gongjing/paris-2016-05/result/SRR2814766/Aligned.out.bam"):
    #out_sorted_bam=input_bam.replace("bam","sorted.bam")
    out_sorted_bam=rreplace(input_bam,"bam","sorted.bam")
    print "[sortbam]",timestamp()
    print "  - input_bam: %s"%(input_bam)
    print "  - out_sorted_bam: %s"%(out_sorted_bam)
    subprocess.call(["samtools sort %s -o %s"%(input_bam,out_sorted_bam)],shell=True)
    print "[sortbam]",timestamp()
    print
    return out_sorted_bam

def indexbam(input_bam="/Share/home/zhangqf/gongjing/paris-2016-05/result/SRR2814766/Aligned.out.bam"):
    print "[indexbam]",timestamp()
    print "  - input_bam: %s"%(input_bam)
    subprocess.call(["samtools index %s"%(input_bam)],shell=True)
    print "[indexbam]",timestamp()
    print 
    return input_bam+'.bai' # bam index suffix by samtools: .bai

def sam2sortedbam(input_sam="/Share/home/zhangqf/gongjing/paris-2016-05/result/SRR2814700/Aligned.out.keep1.sam"):
    bam=sam2bam(input_sam=input_sam)
    sortedbam=sortbam(input_bam=bam)
    bamindex=indexbam(input_bam=sortedbam)
    return sortedbam

def multiBamCorrelationDeeptool(bam_file_ls,binsize,union_out_npz):
    project_pwd,scripts_pwd=get_project_scripts_pwd()
    bam_file_ls=[i if i.startswith("/Share") else project_pwd+'/'+i for i in bam_file_ls] 
    union_out_npz=add_file_pwd(union_out_npz)
    assert union_out_npz.endswith('npz'),"union_out_npz should be .npz: %s"%(union_out_npz)
    union_out_raw=union_out_npz.replace('.npz','.tab')
    union_out_raw=rreplace(union_out_npz,'.npz','.tab')
    save_png=rreplace(union_out_npz,'.npz','.pearson_coor.png')
    save_corr_matrix=rreplace(union_out_npz,'.npz','.pearson_coor_matrix.txt')
    title='/'.join(union_out_npz.split("/")[-2:])
    sample_labels=' '.join(['/'.join(i.split("/")[-2:]) for i in bam_file_ls])
    print "[multiBamCorrelationDeeptool]", timestamp()
    print "  - bam_file_ls:",bam_file_ls
    print "  - binsize: %s"%(binsize)
    print "  - union_out_npz: %s"%(union_out_npz)
    print "  - union_out_raw: %s"%(union_out_raw)
    print "  - title: %s"%(title)
    print "  - save_png: %s"%(save_png)
    print "  - save_corr_matrix: %s"%(save_corr_matrix)
    #subprocess.call(["multiBamSummary bins --bamfiles %s -out %s --outRawCounts %s --binSize=%s"%(' '.join(bam_file_ls),union_out_npz,union_out_raw,binsize)],shell=True)
    subprocess.call(['''plotCorrelation -in %s --corMethod pearson --plotTitle %s --whatToPlot scatterplot -o %s --outFileCorMatrix %s --labels %s'''%(union_out_npz,title,save_png,save_corr_matrix,sample_labels)],shell=True)
    print "[multiBamCorrelationDeeptool]", timestamp()
    print

# bam need sorted
def bamStatsBySamtoolsPlot(input_bam="result/SRR2814700/Aligned.out.sorted.bam"): 
    input_bam=add_file_pwd(input_bam)
    input_bam_stats=input_bam+'.stats'
    print "[bamStatsBySamtoolsPlot]",timestamp()
    print "  - input_bam: %s"%(input_bam)
    print "  - out_bam_stats: %s"%(input_bam_stats)
    subprocess.call(["samtools stats %s > %s"%(input_bam,input_bam_stats)],shell=True)
    input_bam_stats_plot_prefix=plotSamtoolsStats(input_bam_stats=input_bam_stats)
    print "[bamStatsBySamtoolsPlot]",timestamp()
    print
    return input_bam_stats,input_bam_stats_plot_prefix

def samtoolsStatsParse(input_stats=""):
    input_stats=add_file_pwd(input_stats)
    raw_total = 0
    read_len={}
    with open(input_stats,'r') as STAT:
        for line in STAT:
           if "raw total sequence" in line:
                raw_total=line.strip().split(":")[1].strip()
           if line.startswith("RL"):
               arr=line.split("\t")
               read_len[arr[1]] = int(arr[2])
    return int(raw_total),read_len

# may be run on mgt
def plotSamtoolsStats(input_bam_stats=""):
    input_bam_stats=add_file_pwd(input_bam_stats)
    input_bam_stats_plot=input_bam_stats+'.plot'+'/'
    plotBamstats="~/usr/samtools/bin/plot-bamstats"
    print "[plotSamtoolsStats]",timestamp()
    print "  - input_bam_stats: %s"%(input_bam_stats)
    print "  - out_plot_path: %s"%(input_bam_stats_plot)
    try:
        subprocess.call(["%s -p %s %s"%(plotBamstats,input_bam_stats_plot,input_bam_stats)],shell=True)
    except:
        print "  - [error] not plot success: %s"%(plotBamstats)
    print "[plotSamtoolsStats]",timestamp()
    print
    return input_bam_stats_plot

def plotSamtoolsStats2(sample_ls,bam_stats_ls):
    sample_bam_stats_fs=[sample+'/'+bam_stats for sample in sample_ls for bam_stats in bam_stats_ls]
    sample_bam_stats_pwds=[add_file_pwd('result/'+i) for i in sample_bam_stats_fs]
    read_len_ls=[]
    read_df_ls=[]
    for sample_bam_stats_f,sample_bam_stats_pwd in zip(sample_bam_stats_fs,sample_bam_stats_pwds):
        raw_total,read_len = samtoolsStatsParse(sample_bam_stats_pwd)
        read_len_ls.extend(read_len.keys())
    read_len_set=set(read_len_ls)
    for sample_bam_stats_f,sample_bam_stats_pwd in zip(sample_bam_stats_fs,sample_bam_stats_pwds):
        raw_total,read_len = samtoolsStatsParse(sample_bam_stats_pwd)
        read_len_ratio={}
        for i in read_len_set:
            if read_len.has_key(i): 
                read_len_ratio[i] = read_len[i]/float(raw_total)
            else:
                read_len_ratio[i] = 0
        read_df=pd.DataFrame.from_dict(read_len_ratio,orient='index') 
        read_df.sort_index(inplace=True)
        read_df.columns=[sample_bam_stats_f]
        read_df_ls.append(read_df)
        print sample_bam_stats_f
        print raw_total
        print read_len
        print read_len_ratio
        print read_df
    read_df_all=pd.concat(read_df_ls,axis=1)
    print read_df_all
    fig,ax=plt.subplots(figsize=(12,8))
    colors=['r','b','g','k','c','m']
    read_df_all.sort_index(axis=1,inplace=True)
    read_df_all['len']=read_df_all.index
    read_df_all['len']=read_df_all['len'].astype(float)
    read_df_all.sort_values('len',inplace=True)
    for n,col in enumerate(read_df_all.columns[0:-1]):
        print "plot: %s"%(col)
        label=col.split("/")[0]+' '+col.split("/")[1].split(".")[2]
        read_df_all.plot(kind='line',x='len',y=col,color=colors[n],ax=ax,label=label)
        read_df_all.plot(kind='scatter',x='len',y=col,color=colors[n],ax=ax)
    plt.ylabel('frequency')
    plt.ylim(0,)
    xmin=min(read_df_all['len'])-1
    xmax=max(read_df_all['len'])+1
    plt.xlim(xmin,xmax)
    plt.xlabel("mapped reads length")
    plt.legend(loc='upper left')
    savefn='result'+'/'+'plots'+'/'+'map_len_dist_'+'_'.join(sample_ls)+'.png'
    savefn=add_file_pwd(savefn)
    fig.savefig(savefn)
    plt.close()
    

def plotSamtoolsStatsCustome(sample_ls):
    gap_bam_ls=['Aligned.out.gapped.sorted.bam']
    nongap_bam_ls=['Aligned.out.nongapped.sorted.bam']
    gap_nongap_bed_ls=['gap_1_nongap_0.bed','gap_1_nongap_1.bed','gap_0_nongap_1.bed']
    sample_dict_gap={}
    gap_stats_ls=[]
    for s in sample_ls:
        for gap_bam in gap_bam_ls:
             for gap_nongap_bed in gap_nongap_bed_ls[0:2]:
                 stats_f=add_file_pwd('result'+'/'+s+'/'+gap_bam+'_'+gap_nongap_bed+'.bam.stats')
 	         gap_stats_ls.append(stats_f)
                 raw_total=samtoolsStatsParse(stats_f)
                 if not sample_dict_gap.has_key(s+':gap'): 
                     sample_dict_gap[s+':gap']=[raw_total]
		 else:
                     sample_dict_gap[s+':gap'].append(raw_total)
    nongap_stats_ls=[]
    sample_dict_nongap={}
    for s in sample_ls:
        for nongap_bam in nongap_bam_ls:
             for gap_nongap_bed in gap_nongap_bed_ls[1:]:
                 stats_f=add_file_pwd('result'+'/'+s+'/'+nongap_bam+'_'+gap_nongap_bed+'.bam.stats')
 	         nongap_stats_ls.append(stats_f)
                 raw_total=samtoolsStatsParse(stats_f)
                 if not sample_dict_nongap.has_key(s+':nongap'): 
                     sample_dict_nongap[s+':nongap']=[raw_total]
		 else:
                     sample_dict_nongap[s+':nongap'].append(raw_total)
    print sample_dict_gap
    print sample_dict_nongap
    out_f='result'+'/'+'plots'+'/'+'gap_nongap.'+'_'.join(sample_ls)+'.png'
    out_f=add_file_pwd(out_f)
    fig,ax=plt.subplots(nrows=1,ncols=2,figsize=(12,6))
    df1=pd.DataFrame.from_dict(sample_dict_gap,orient='index')
    df1.columns=['only_in_gap','both_in_gap_nongap']
    df2=pd.DataFrame.from_dict(sample_dict_nongap,orient='index')
    df2.columns=['both_in_gap_nongap','only_in_nongap']    
    colors = ['b', 'g', 'r','k','c','m']
    df1=df1.sort_index(axis=0)
    df2=df2.sort_index(axis=0)
    df1.plot(kind='barh',ax=ax[0],stacked=True,color=colors[0:2])
    ax[0].ticklabel_format(style='sci', axis='x', scilimits=(0,0))  # scientific x/y ticks
    df2.plot(kind='barh',ax=ax[1],stacked=True,color=colors[1:3])
    fig.suptitle("reads abundance in gapped & nongapped overlaped region")
    ax[1].yaxis.tick_right() # set x/y axis tick left/right/top/bottom
    fig.tight_layout()
    fig.savefig(out_f)

def starMapLogPlot(sample_ls):
    log_ls=[add_file_pwd('result'+'/'+s+'/'+'Log.final.out') for s in sample_ls]
    out_f='result'+'/'+'plots'+'/'+'starLog.'+'_'.join(sample_ls)+'.png'
    out_f=add_file_pwd(out_f)
    printFuncRun('starMapLogPlot')
    printFuncArgs()
    map_info_dict={}
    for s,log in zip(sample_ls,log_ls):
        with open(log,'r') as LOG:
            num_index=[5,6,9,24,29]
            content=LOG.readlines()
            content=getElementsFromList(content,num_index)
            readsNum,readsLen,uniq_ratio,multi_ratio,short_ratio = [i.split("\t")[1].strip() for i in content]
            map_info_dict[s+'\n'+readsNum+' '+readsLen]=[float(i.strip("%")) for i in [uniq_ratio,multi_ratio,short_ratio]]
    print map_info_dict
    print
    dict_bar_plot_by_dataframe_mutlti(d=map_info_dict,columns=['uniq_map','multi_map','unmap_short'],y_label="sample mapped ratio by category",x_label="sample",title_str="STAR map info",save_fn=out_f,orient='h',stacked=1,fig_size_x=12,fig_size_y=12)
    return map_info_dict

def samMapPlot(sample_ls):
    map_sam_ls=[add_file_pwd('result'+'/'+s+'/'+'Aligned.out.sam') for s in sample_ls]
    map_chimeric_ls=[add_file_pwd('result'+'/'+s+'/'+'Chimeric.out.sam') for s in sample_ls]
    out_f='result'+'/'+'plots'+'/'+'starMapDis.'+'_'.join(sample_ls)+'.png'
    out_f2='result'+'/'+'plots'+'/'+'starMapDis.'+'_'.join(sample_ls)+'.2.png'
    out_f=add_file_pwd(out_f)
    out_f2=add_file_pwd(out_f2)
    printFuncRun('samMapPlot')
    printFuncArgs()
    map_dict=collections.defaultdict(list)
    map_dict_2=collections.defaultdict(list)
    for sample,map_sam,map_chimeric in zip(sample_ls,map_sam_ls,map_chimeric_ls):
        keep_f,all_line,uniq_map,multi_map,gapped_reads = sam_keep_alignment_first_n_hits(input_sam=map_sam,HIn=1,NHn=100)
        #map_dict[sample+' '+'uniq'] = uniq_map
        #map_dict[sample+' '+'multi'] = multi_map
        #map_dict[sample+' '+'gapped'] = gapped_reads
        map_dict[sample].extend([uniq_map,multi_map,gapped_reads])
        map_dict_2[sample].extend([uniq_map+multi_map,gapped_reads])
        keep_f,all_line,uniq_map,multi_map,gapped_reads = sam_keep_alignment_first_n_hits(input_sam=map_chimeric,HIn=1,NHn=100)
        #map_dict[sample+' '+'chimeric'] = multi_map
        map_dict[sample].append(multi_map)
        map_dict_2[sample][1] += multi_map
    print map_dict
    print map_dict_2
    dict_bar_plot_by_dataframe_mutlti(d=map_dict,columns=['uniq_map','multi_map','gapped_reads','chimeric_reads'],y_label="sample mapped reads num by category",x_label="sample",title_str="STAR map info",save_fn=out_f,orient='h',stacked=1,fig_size_x=12,fig_size_y=12)
    dict_bar_plot_by_dataframe_mutlti(d=map_dict_2,columns=['uniq_map+multi_map','gapped_reads+chimeric_reads'],y_label="sample mapped reads num by category",x_label="sample",title_str="STAR map info",save_fn=out_f2,orient='h',stacked=1,fig_size_x=12,fig_size_y=12)
    return map_dict

# calculate genome wide single base coverage from bam_ls, calculator: bedtools,igvtools
# has checked, both bedtools & igvtools get the same result for a rnaseq mapped bam
# bedtools: all base for each chr; igvtools: only bases with coverage >0   
# genome wide single base coverage related
def genome_cov_base(bam_ls,calculator,origanism='human'):
    genome_fa='/Share/home/zhangqf/gongjing/paris-2016-05/result/sunlei/data/index/'+origanism+'/transcriptome/transcriptome.fa'
    genome_chr_size='/Share/home/zhangqf/gongjing/paris-2016-05/result/sunlei/data/index/'+origanism+'/transcriptome/chrNameLength.txt'

    printFuncRun('genome_cov_base')
    printFuncArgs() 
    calculator_avaliable_ls=['bedtools','igvtools']
    if calculator not in calculator_avaliable_ls:
        print "[error] genome wide coverage per base calculator not available for : %s"%(calculator)
        print "calculator_avaliable_ls:",calculator_avaliable_ls
        return
    for bam in bam_ls:
        bam_cov = bam+'.'+calculator+'.wig'
        if calculator == 'bedtools':
            subprocess.call(["bsub -q Z -o %s bedtools genomecov -ibam %s -g %s -split -d"%(bam_cov,bam,genome_chr_size)],shell=True)
        if calculator == 'igvtools':
            bam_bai=bam+'.bai'
            if not os.path.isfile(bam_bai):
                subprocess.call(["bsub -q Z samtools index %s"%(bam)],shell=True)
            subprocess.call(["bsub -q Z igvtools count -w 1 %s %s %s"%(bam,bam_cov,genome_fa)],shell=True)
            bedtools_genomecov_out_convert(cov=bam_cov,calculator='igvtools',chr_len=genome_chr_size) 
    printFuncRun('genome_cov_base') 
    return 

# for both bedtools,igvtools calculated base.cov, covert to a tab separated txt: chr   base_cov(space separated), and dict: key:chr,value:[covarage list]
# for compare with icshape, easy to retrieve coverage(count) for any postion on any chr(transcript)
def bedtools_genomecov_out_convert(cov,calculator='bedtools',chr_len=None):
    printFuncRun('bedtools_genomecov_out_convert')
    printFuncArgs()
    
    transcript,val_ls = '',[]
    cov_convert = cov+'.txt'
    
    if calculator == 'bedtools':
        transcript_dict=defaultdict(list)  # key: transcript val:[coverage list]
        with open(cov,'r') as COV:
            for n,line in enumerate(COV):  
                if not line or line.startswith('#'): continue
                transcript,pos,val = line.strip().split('\t')
                transcript_dict[transcript].append(val)
                if n%1000000 == 0: print " - process: %s"%(n)
                
    if calculator == 'igvtools':
        chr_len_dict = read_genome_chr_len(chr_len=chr_len)
        print "construct transcript_dict based on length",timestamp()
        transcript_dict = {chr:([0]*length) for chr,length in chr_len_dict.items()}
        print "construct transcript_dict based on length",timestamp()
        with open(cov,'r') as COV:
            for n,line in enumerate(COV):
                if not line or line.startswith(('#','track')): continue
                if line.startswith('variableStep'):
                    transcript = line.strip().split(' ')[1].split('=')[1]
                    continue
                pos,val = line.strip().split('\t')
                try:
                    transcript_dict[transcript][int(pos)-1] = val
                except:
                    print transcript,pos,len(transcript_dict[transcript])
                if n%1000000 == 0: print " - process: %s"%(n)
                    
    with open(cov_convert,'w') as COV_CONVERT:
        for n,(key,val) in enumerate(transcript_dict.items()):
           #if n <= 10: print key+'\t'+' '.join(map(str,val))
           if any(float(i)>0 for i in val):
               print >>COV_CONVERT,key+'\t'+' '.join([str(i).split('.')[0] for i in val])
           
    COV_CONVERT.close()
    printFuncRun('bedtools_genomecov_out_convert')
    return cov_convert,transcript_dict

def transCoor2genomeCoor(trans_dict, trans_name, site):
    site = int(site)
    if not trans_dict.has_key(trans_name):
        print "no such tx: %s, %s"%(trans_name, site)
        return False, False
    exon_tuple = trans_dict[trans_name]['exon'].split(',')
    exon_list = list()
    for mtuple in exon_tuple:
        exon_list.append( [int(i) for i in mtuple.split('-')] )
    strand = trans_dict[trans_name]['strand']
    accu = 0
    for i in range(len(exon_list)):
        range_of_tuple = abs( exon_list[i][1] - exon_list[i][0] ) + 1
        if accu + range_of_tuple + 1 > site:
            if strand == '+':
                return trans_dict[trans_name]['chr'],  exon_list[i][0]+site-accu-1
            else:
                return trans_dict[trans_name]['chr'],  exon_list[i][1]-(site-accu-1)
        accu += range_of_tuple
    print "tx query site > len: %s, %s"%(trans_name, site)
    return False, False

def loadGTFBed(file_name=None, organism='human'):
    printFuncRun('loadGTFBed')
    # /Share/home/zhangqf/lipan/DYNAMIC/GTF/Homo_sapiens.trans.bed
    if file_name is not None:
        pass
    elif organism == "human":
        file_name = '/Share/home/zhangqf5/gongjing/paris-2016-05/data/reference_data/GENCODE/human/gencode.v24.primary_assembly.annotation.bed'
    elif organism == 'mouse':
        file_name = '/Share/home/zhangqf5/gongjing/paris-2016-05/data/reference_data/GENCODE/mouse/gencode.vM9.primary_assembly.annotation.bed'
    else:
        return 
    printFuncArgs()
    H = open(file_name)
    line = H.readline()
    trans_dict = dict()
    while line: 
        if line.startswith('#'): line = H.readline(); continue 
        arr = line.strip().split()
        trans_dict[ arr[5] ] = { 'gene':arr[4], 'exon':arr[7], 'chr':arr[0], 'strand':arr[3] }
        line = H.readline()
    H.close()
    printFuncRun('loadGTFBed')
    return trans_dict

    # transBed = loadGTFBed('/Share/home/zhangqf/lipan/DYNAMIC/GTF/Homo_sapiens.trans.bed')
    # transCoor2genomeCoor(transBed, 'ENST00000365328', 133)
    # ('chr5', 139278905)


##################################################
### 
def plot_trun_site_mutation(input_f="/Share/home/zhangqf/gongjing/paris-2016-05/result/SRR2814710/Aligned.out.native.mutation.txt.ext.trun"):
    df=pd.read_csv(input_f,sep=" ",index_col=0)
    fig=plt.figure(2,figsize=(16,16),dpi=80)
    fig_size_x=12
    fig_size_y=8
    df2=df[(df.T != 0).any()].T.iloc[:,3:].sort_index(axis=0)
    df2.plot(kind='barh',stacked=True,figsize=(fig_size_x,fig_size_y),fontsize=14,width=0.8,position=0.5)
    print df
    print df2
    title_str='/'.join(input_f.split('/')[-2:])+'\n'+'truncation site mutation type ditribution'
    x_label="mutation type"
    y_label="count"
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    sum_row = {col: df2[col].sum() for col in df2}
    sum_all = sum(sum_row.values())
    title_str = title_str+' '+str(sum_all)
    plt.title(title_str)
    #fig.tight_layout()
    plt.gcf().subplots_adjust(bottom=0.15)
    file_str=input_f+'_bar.png'
    plt.savefig(file_str,)
    plt.close()

##################################################
### bar plot of a column of a dataframe, zoom in local region
def plot_frequency(df,col_str,sep_val,savefn,input_fn):
    fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(16,16))

    df[col_str].plot(kind='hist',bins=100,ax=ax)
    ax.grid(False)
    ax.set_title(input_fn+':'+col_str+','+str(df.shape[0]))
    ax.set_xlabel(col_str)

    ax1=fig.add_axes([0.35,0.55,0.5,0.3])
    df2=df.loc[df[col_str]<=sep_val,:]
    df2[col_str].plot(kind='hist',bins=50,ax=ax1)
    ax1.set_title(col_str+'<='+str(sep_val)+': '+str(df2.shape[0]))

    ax2=fig.add_axes([0.35,0.2,0.5,0.3])
    df3=df.loc[df[col_str]>sep_val,:]
    df3.loc[:,col_str].plot(kind='hist',bins=200,ax=ax2)
    ax2.set_title(col_str+'>'+str(sep_val)+': '+str(df3.shape[0]))

    #plt.show()

    fig.savefig(savefn)
    plt.close()
    
    print "======"
    print "plot bar graph:"
    print "    - col: ",col_str
    print "    - sep_val: ",sep_val
    print "    - savefn: ",savefn
    print
    
##################################################
### plot df two columns by min,25%,50%,75%,max seaborn jointplot
def plot_df_two_cols_describe_jointplot(df,x,y,fn):
        df[x]=df[x].astype(float)
        df[y]=df[y].astype(float)
        print "[plot_df_two_cols_describe_jointplot]"
        x_step_ls=df.describe().loc['min':'max',x]
	y_step_ls=df.describe().loc['min':'max',y]
	print "%s percentile: "%(x),x_step_ls
	print "%s percentile: "%(y),y_step_ls
        x_step_ls=refine_list(x_step_ls)
        y_step_ls=refine_list(y_step_ls)
	print "%s percentile[after refine]: "%(x),x_step_ls
	print "%s percentile[after refine]: "%(y),y_step_ls
        df_percentile_ls=[]
        plot_n_ls=[]

	for m,i in enumerate(x_step_ls):
	    for n,j in enumerate(y_step_ls):
		if m > 0 and n > 0:
		    print "%s:%s, %s:%s"%(x,m,y,n)
		    print "  - %s,%s"%(x_step_ls[m-1],i)
		    print "  - %s,%s"%(y_step_ls[n-1],j)
		    if x_step_ls[m-1] != i:
			tmp=df.loc[(df[x]>=x_step_ls[m-1]) & (df[x] < i),:]
		    else:
			tmp=df.loc[df[x]==i,:]
		    if y_step_ls[n-1] != j:
			tmp=tmp.loc[(tmp[y]>=y_step_ls[n-1]) & (tmp[y]<j),:]
		    else:
			tmp=tmp.loc[tmp[y]==j,:]
		    print "  - shape: ",tmp.shape
		    #print tmp.head()
                    d=collections.Counter(tmp['refBase']) # just for ivtools count
                    print "  - refbase count:",d
                    dd=["%s:%.2f"%(ii,jj/float(sum(d.values()))) for ii,jj in d.items()]
                    print "  - refbase frequency:",dd
		    n_point=str(tmp.shape[0])
		    x_str=x+':'+str(0.25*m-0.25)+'-'+str(0.25*m)
		    y_str=y+':'+str(0.25*n-0.25)+'-'+str(0.25*n)
		    try:
			a=sns.jointplot(x=x,y=y,data=tmp)
			#a.fig.text(0.75,0.85,'n='+n_point+'\n'+x_str+'\n'+y_str,color='b')
			a.fig.text(0.75,0.85,'n='+n_point,color='b',fontsize=16)
                        a.fig.suptitle(' '.join(dd))
			savefn=fn+'_'+x+'_'+y+'_'+str(m)+str(n)+'.png'
			a.savefig(savefn)
                        plt.close()
		    except:
			print "seaborn not plot success"
			fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(6,6))
			#fig.text(0.75,0.85,'n='+n_point+'\n'+x_str+'\n'+y_str,color='b')
			fig.text(0.75,0.85,'n='+n_point,color='b',fontsize=16)
			savefn=fn+'_'+x+'_'+y+'_'+str(m)+str(n)+'.png'
			#plt.show()
			fig.savefig(savefn)
			plt.close()
		    print
                    df_percentile_ls.append(tmp)
                    plot_n_ls.append(str(m)+str(n))

	subprocess.call(["montage %s{14,24,34,44,13,23,33,43,12,22,32,42,11,21,31,41}.png -geometry 360x360 %s"%(fn+'_'+x+'_'+y+'_',fn+'_'+x+'_'+y+'_'+'percentile_all.png')],shell=True)
        return df_percentile_ls,plot_n_ls



##################################################
### fasta file related
def prepend_str_before_f(fn, str=""):
    subprocess.call(["""sed -i '1,1s/^/%s/g' %s"""%(str.encode('string-escape'),fn)],shell=True)
    """
    with open(fn, 'w') as FN:
         FN.write(str)
         FN.flush()
         subprocess.call(["grep epoch full"], stdout=FN)
    """


def read_fa_into_dict(origanism='mouse',fa=None):
    """
    project_pwd,scripts_pwd=get_project_scripts_pwd()
    if origanism == "human":
        ref_fa="data/reference_data/GENCODE/human/GRCh38.primary_assembly.genome.chr.fa"
    else:
        ref_fa="data/reference_data/GENCODE/mouse/GRCm38.primary_assembly_chr.genome.fa"
    """
    if not fa is None:
       ref_fa = fa
    else:
       ref_fa=get_origanism_fa_path(origanism)
    print "[read_fa_into_dict]",timestamp()
    print "  - ref_fa: %s"%(ref_fa)
    ref_fa_dict=Fasta(ref_fa)
    print "[read_fa_into_dict]",timestamp()
    random_key = ref_fa_dict.keys()[0]
    print "random 10 keys:",random_key
    if len(random_key.strip().split()) > 1:
        print "whole > line as key, will refine"
        ref_fa_dict_refine = {key.split()[0]:val for key,val in ref_fa_dict.items()}
        random_key = ref_fa_dict_refine.keys()[0]
        print "random 10 keys after refine:",random_key
        return ref_fa_dict_refine
    print
    return ref_fa_dict

def get_origanism_fa_path(origanism='mouse'):
    project_pwd,scripts_pwd=get_project_scripts_pwd()
    if origanism == "human":
        ref_fa="data/reference_data/GENCODE/human/GRCh38.primary_assembly.genome.chr.fa"
    elif origanism == "mouse":
        ref_fa="data/reference_data/GENCODE/mouse/GRCm38.primary_assembly.genome.chr.fa"
    else:
        print "[error] origanism .fa not found: %s"%(origanism)
    ref_fa_path=project_pwd+'/'+ref_fa
    return ref_fa_path

def read_genome_chr_len(chr_len='/Share/home/zhangqf/gongjing/paris-2016-05/result/sunlei/data/index/mouse/transcriptome/chrNameLength.txt',origanism='mouse'):
    chr_len_dict={}
    with open(chr_len,'r') as CHR_LEN:
        for line in CHR_LEN:
            if not line or line.startswith('#'): continue
            chr,length = line.strip().split('\t')
            chr_len_dict[chr]=int(length)
    print "read into dict: %s"%(chr_len)
    return chr_len_dict



##################################################
### snp file related
def read_snp_into_dict(origanism='mouse',source='ucsc'):
    snp_f_dict={"mouse+ucsc":"data/ucsc_snp/MouseSnp142Common.txt",
                "mouse+dbsnp":"data/dbsnp/mouse_10090/mouse_10090_v146.bed",
                "human+ucsc":"data/ucsc_snp/HumanSnp142Common.txt",
                "human+dbsnp":"data/dbsnp/human_9606_b147_GRCh38p2/human_v147_all.bed"}
    snp_f=snp_f_dict[origanism+'+'+source]
    snp_f=add_file_pwd(snp_f)
    read_cols=[]
    if 'ucsc' in snp_f:
        read_cols=[1,2,3] # 0-based column num
    if 'dbsnp' in snp_f:
        read_cols=[0,1,2]
    snp_dict={} # key: chr+start+end
    print "[read_snp_into_dict]",timestamp()
    print "  - origanism: %s"%(origanism)
    print "  - source: %s"%(source)
    print "  - snp_f: %s"%(snp_f)
    with open(snp_f,'r') as SNP:
        for line in SNP:
            arr=line.strip().split()
            key=' '.join([arr[i].strip() for i in read_cols])
            snp_dict[key]=1
    print "  - len: %s"%(len(snp_dict))
    print "[read_snp_into_dict]",timestamp()
    print
    return snp_dict,len(snp_dict)


##################################################
### bed/fasta 
# get sequence from fasta in bed defined region, support extend 5'/3' region
def getFastaFromBed(input_bed="scripts/test_dataset/hum362.bed",origanism='human',extend_left=0,extend_right=0,outfmt='fa',strand=0):
    project_pwd,scripts_pwd=get_project_scripts_pwd()
    input_bed = add_file_pwd(input_bed)
    if origanism == 'human':
        ref_fa=project_pwd+'/'+'data/reference_data/GENCODE/human'+'/'+'GRCh38.primary_assembly.genome.chr.fa'
    elif origanism == 'mouse':
        ref_fa=project_pwd+'/'+'data/reference_data/GENCODE/mouse'+'/'+'GRCm38.primary_assembly.genome.chr.fa'
    else:
        print "[origanism error] unavailable origanism: %s"%(origanism)
    #output=input_bed.replace("bed",str(extend_left)+'_'+str(extend_right)+'_'+str(strand)+'.'+outfmt)
    output=rreplace(input_bed,"bed",str(extend_left)+'_'+str(extend_right)+'_'+str(strand)+'.'+outfmt)
    if extend_left != 0 or extend_right != 0:
        #bed_tmp=input_bed.replace("bed","tmp.bed")
        bed_tmp=rreplace(input_bed,"bed","tmp.bed")
        BED_TMP=open(bed_tmp,'w')
        with open(input_bed,'r') as INPUT_BED:
            for line in INPUT_BED:
                if not line.startswith("#"):
                    arr=line.strip().split()
                    try:
                        arr[1]=str(int(arr[1])-extend_left)
                        arr[2]=str(int(arr[2])+extend_right)
                        print >>BED_TMP,'\t'.join(arr)
                    except:
                        print line,arr[1],arr[2]
                        print traceback.print_exc()
        input_bed=bed_tmp
        BED_TMP.close() # keepInMind: close file handle before subprocess.call
    print "[getFastaFromBed]",timestamp()
    print "  - input_bed: %s"%(input_bed)
    print "  - origanism: %s"%(origanism)
    print "  - extend_left: %s"%(extend_left)
    print "  - extend_right: %s"%(extend_right)
    print "  - outfmt: %s"%(outfmt)
    print "  - output: %s"%(output)
    print "  - strand: %s"%(strand)
    print "  - ref_fa: %s"%(ref_fa)  
    if outfmt == "txt":
         if strand:
             subprocess.call(["bedtools getfasta -fo %s -tab -s -fi %s -bed %s"%(output,ref_fa,input_bed)],shell=True)        
         else:
             subprocess.call(["bedtools getfasta -fo %s -tab -fi %s -bed %s"%(output,ref_fa,input_bed)],shell=True)        
    if outfmt == 'fa':
         if strand:
             print "bedtools getfasta -fo %s -s -fi %s -bed %s"%(output,ref_fa,input_bed)
             subprocess.call(["bedtools getfasta -fo %s -s -fi %s -bed %s"%(output,ref_fa,input_bed)],shell=True)        
         else:
             subprocess.call(["bedtools getfasta -fo %s -fi %s -bed %s"%(output,ref_fa,input_bed)],shell=True)        
    print
    if extend_left != 0 or extend_right != 0: 
        subprocess.call(["rm %s"%(bed_tmp)],shell=True)
    print "[getFastaFromBed]",timestamp()
    print
    return output

def bedtoolsIntersect(input_bed1="result/xiongtuanlin_result/mou362.bed",input_bed2="data/supplementary_data/icshape_data/mES/icSHAPE.sim.bedgraph",intersect_bed="",mode="loj",strand=0):
    project_pwd,scripts_pwd=get_project_scripts_pwd()
    input_bed1 = add_file_pwd(input_bed1)
    input_bed2 = add_file_pwd(input_bed2)
    if intersect_bed == "":
        #intersect_bed=input_bed1.replace("bed","intersect.bed")
        intersect_bed=rreplace(input_bed1,"bed","intersect.bed")
    intersect_bed = add_file_pwd(intersect_bed)
    printFuncArgs()
    print "[bedtoolsIntersect]",timestamp()
    print "  - bed1: %s"%(input_bed1)
    print "  - bed2: %s"%(input_bed2)
    if mode=="loj":
        subprocess.call(['''bedtools intersect -a %s -b %s -wa -wb -loj > %s'''%(input_bed1,input_bed2,intersect_bed)],shell=True)
    if strand and mode=="loj":
        subprocess.call(['''bedtools intersect -a %s -b %s -wa -wb -loj -s > %s'''%(input_bed1,input_bed2,intersect_bed)],shell=True)
    print "  - out intersect bed: %s"%(intersect_bed)
    print "[bedtoolsIntersect]",timestamp()
    print
    return intersect_bed
 
def bedOverlapBySetOperation(input_bed1="result/xiongtuanlin_result/mou362.bed",input_bed2="data/supplementary_data/icshape_data/mES/icSHAPE.sim.bedgraph",intersect_bed_info_txt="",plotvenn=0,plotdist=0):
    input_bed1 = add_file_pwd(input_bed1)
    input_bed2 = add_file_pwd(input_bed2)
    print "[bedOverlapBySetOperation]",timestamp()
    print "  - bed1: %s"%(input_bed1)
    print "  - bed2: %s"%(input_bed2)
    bed1_dict=read_file_2_dict(input_f=input_bed1,key_col=[0,1,2],sep='\t',val_col=[18])
    bed2_dict=read_file_2_dict(input_f=input_bed2,key_col=[0,1,2],sep='\t',val_col=[3])
    bed1_set=set(bed1_dict.keys())
    bed2_set=set(bed2_dict.keys())
    bed1_1_bed2_1=bed1_set & bed2_set
    bed1_1_bed2_0=bed1_set - bed2_set
    bed1_0_bed2_1=bed2_set - bed1_set
    if plotvenn:
         venn3plot(mode='string',subsets_ls=(bed1_set,bed2_set),labels_ls=['/'.join(input_bed1.split("/")[-2:]),'/'.join(input_bed2.split("/")[-2:])],title_str="venn plot",save_fn=input_bed1+'_'+input_bed2.split("/")[-1]+'.overlap.venn.png')
    if plotdist:
         tmp=input_bed1+'_'+input_bed2.split("/")[-1]+'.overlap.tmp'
         TMP=open(tmp,'w')
         for i in bed1_1_bed2_0:
             print >>TMP,i+'\t'+bed1_dict[i]+'\t'+str(-1)+'\t'+'only_paris'
         for i in bed1_1_bed2_1:
             print >>TMP,i+'\t'+bed1_dict[i]+'\t'+bed2_dict[i]+'\t'+'both_paris_icshape'
         for i in bed1_0_bed2_1:
             print >>TMP,i+'\t'+str(-1)+'\t'+bed2_dict[i]+'\t'+'only_icshape'
         TMP.close()
         df=pd.read_csv(tmp,sep='\t',header=None)
         df.columns=['chr','start','end','paris_coverage','icshape_val','type']
         print df.head()
         print df.describe()
         paris_max_cutoff=500
         paris_min_cutoff=1
         df_paris=df[(df['type'] != 'only_icshape') & (df['paris_coverage']<paris_max_cutoff) & (df['paris_coverage']>paris_min_cutoff)]
         df_paris_2=df[(df['type'] != 'only_icshape') & (df['paris_coverage']<100) & (df['paris_coverage']>paris_min_cutoff)]
         tmp_log=open(tmp+'.log.csv','w')
        
         print >>tmp_log,"paris",'only'
         print >>tmp_log,df_paris.loc[df_paris['type']=='only_paris','paris_coverage'].describe()
         print >>tmp_log,"paris",'both'
         print >>tmp_log,df_paris.loc[df_paris['type']=='both_paris_icshape','paris_coverage'].describe()
         print >>tmp_log,"paris100",'only'
         print >>tmp_log,df_paris_2.loc[df_paris_2['type']=='only_paris','paris_coverage'].describe()
         print >>tmp_log,"paris100",'both'
         print >>tmp_log,df_paris_2.loc[df_paris_2['type']=='both_paris_icshape','paris_coverage'].describe()
         df_icshape=df[df['type'] != 'only_paris']
         print >>tmp_log,"icshape",'only'
         print >>tmp_log,df_icshape.loc[df_icshape['type']=='only_icshape','icshape_val'].describe()
         print >>tmp_log,"icshape",'both'
         print >>tmp_log,df_icshape.loc[df_icshape['type']=='both_paris_icshape','icshape_val'].describe()

         df_sns_violinplot(df=df_paris,col_str='paris_coverage',savefn=tmp+'.parisViolin.png',class_col="type",orient="")	
         df_sns_violinplot(df=df_paris_2,col_str='paris_coverage',savefn=tmp+'.parisViolin100.png',class_col="type",orient="")	
         df_sns_violinplot(df=df_icshape,col_str='icshape_val',savefn=tmp+'.icshapeViolin.png',class_col='type',orient="")
    print "%s & %s: %s"%(input_bed1,input_bed2,len(bed1_1_bed2_1))
    print "%s - %s: %s"%(input_bed1,input_bed2,len(bed1_1_bed2_0))
    print "%s - %s: %s"%(input_bed2,input_bed1,len(bed1_0_bed2_1))
    if intersect_bed_info_txt != "":
        intersect_bed_info_txt = add_file_pwd(intersect_bed_info_txt)
        with open(intersect_bed_info_txt,'w') as OL:
            print >> OL,len(bed1_1_bed2_0),len(bed1_1_bed2_1),len(bed1_0_bed2_1)
    print "[bedOverlapBySetOperation]",timestamp()
    print
    return bed1_1_bed2_0,bed1_1_bed2_1,bed1_0_bed2_1

def bedtoolsBedBamOverlap(input_bam="result/SRR2814766/Aligned.out.gapped.sorted.bam",input_bed="result/SRR2814766/gap_1_nongap_1.bed",out_bam=''):
    input_bam=add_file_pwd(input_bam)
    input_bed=add_file_pwd(input_bed)
    if out_bam == "":
        out_bam=input_bam+'_'+input_bed.split("/")[-1]+'.bam'
    out_bam=add_file_pwd(out_bam)
    print "[bedtoolsBedBamOverlap]",timestamp()
    print "  - input_bam: %s"%(input_bam)
    print "  - input_bed: %s"%(input_bed)
    print "  - out_bam: %s"%(out_bam)
    subprocess.call(["bedtools intersect -abam %s -b %s > %s"%(input_bam,input_bed,out_bam)],shell=True) 
    print "[bedtoolsBedBamOverlap]",timestamp()
    print
    return out_bam

def bedAnnotationByGTF(input_bed="result/SRR2814700/Aligned.out.sorted.bam.igvtoolsCallStrand.all.bed",anno_type='CDS:UTR:start_codon:stop_codon',origanism='mouse',anno_bed=""):
    input_bed=add_file_pwd(input_bed)
    if origanism == 'mouse':
        gtfprefix='data/reference_data/GENCODE/mouse' # gencode.vM9.primary_assembly.annotation.CDS.gtf
        ver='vM9'
    if origanism == 'human':
        gtfprefix='data/reference_data/GENCODE/human' # gencode.v24.primary_assembly.annotation.CDS.gtf
        ver='v24'
    if ':' in anno_type:
        anno_type=anno_type.split(":")
    else:
        anno_type=[anno_type] # convert to list
    if anno_bed == "":
        anno_bed=rreplace(input_bed,'bed','_'.join(anno_type)+'.bed')
    anno_bed=add_file_pwd(anno_bed)
    anno_gtfs=[gtfprefix+'/'+'gencode.'+ver+'.primary_assembly.annotation.'+type+'.uniq.gtf' for type in anno_type]
    anno_gtfs=[add_file_pwd(gtf) for gtf in anno_gtfs]
    print "[bedAnnotationByGTF]",timestamp()
    print "  - input_bed: %s"%(input_bed)
    print "  - anno_type: %s"%(anno_type)
    print "  - origanism: %s"%(origanism)
    print "  - anno_gtfs: ",anno_gtfs
    print "  - out anno_bed: %s"%(anno_bed)
    tmp_bed=input_bed
    for type,anno_gtf in zip(anno_type,anno_gtfs):
        tmp_out_bed=tmp_bed+'.'+type
        #subprocess.call(["bedtools intersect -a %s -b %s -bed -wa -wb -loj > %s"%(input_bed,' '.join(anno_gtfs),anno_bed)],shell=True)
        subprocess.call(["bedtools intersect -a %s -b %s -bed -wa -wb -loj > %s"%(tmp_bed,anno_gtf,tmp_out_bed)],shell=True)
        tmp_bed=tmp_out_bed
    subprocess.call(["cp %s %s"%(tmp_out_bed,anno_bed)],shell=True)
    print "[bedAnnotationByGTF]",timestamp()
    print
    return anno_bed

def bedAnnotationBySnpBed(input_bed="result/SRR2814700/Aligned.out.sorted.bam.igvtoolsCallStrand.all.bed",origanism='mouse',source='dbsnp',mode="loj",intersect_bed=""):
    input_bed = add_file_pwd(input_bed)
    if source == "ucsc":
        snpPrefix="data/ucsc_snp"
        snpBed=snpPrefix+'/'+origanism+'Snp142Common.bed'
    if source == 'dbsnp':
        snpPrefix="data/dbsnp"
        snpBed=snpPrefix+'/'+'dbsnp.'+origanism+'.bed'
    snpBed=add_file_pwd(snpBed)
    if intersect_bed == "":
        #intersect_bed=input_bed+'.'+source+'.bed'
        intersect_bed=rreplace(input_bed,'bed',source+'.bed')
    intersect_bed=add_file_pwd(intersect_bed)
    print "[bedAnnotationBySnpBed]",timestamp()
    print "  - input_bed: %s"%(input_bed)
    print "  - origanism: %s"%(origanism)
    print "  - source: %s"%(source)
    print "  - annotated snp file: %s"%(snpBed)
    if mode=="loj":
        subprocess.call(['''bedtools intersect -a %s -b %s -wa -wb -loj > %s'''%(input_bed,snpBed,intersect_bed)],shell=True)
    print "  - out annotated bed: %s"%(intersect_bed)
    print "[bedAnnotationBySnpBed]",timestamp()
    print
    return intersect_bed

def bedAnnotationByIcshape(input_bed="result/SRR2814700/Aligned.out.sorted.bam.igvtoolsCallStrand.all.bed",origanism='mouse',mode="loj",intersect_bed=""):
    input_bed=add_file_pwd(input_bed)
    if origanism == "mouse":
        icshape_bed="data/supplementary_data/icshape_data/mES/icSHAPE.sim.bedgraph"
    if origanism == "human":
        icshape_bed="data/supplementary_data/icshape_data/human/icSHAPE.sim.bedgraph"
    icshape_bed=add_file_pwd(icshape_bed)
    if intersect_bed == "":
        intersect_bed=input_bed+'_icshape'+'.bed'
        intersect_bed=rreplace(input_bed,'bed','icshape'+'.bed')
    intersect_bed=add_file_pwd(intersect_bed)
    printFuncRun('bedAnnotationByIcshape')
    printFuncArgs()
    if mode=="loj":
        subprocess.call(['''bedtools intersect -a %s -b %s -wa -wb -loj > %s'''%(input_bed,icshape_bed,intersect_bed)],shell=True)
    printFuncRun('bedAnnotationByIcshape')
    print
    return intersect_bed 

def gencode_biotype_load(txt='/Share/home/zhangqf5/gongjing/paris-2016-05/data/reference_data/GENCODE/gencode_biotype_table.txt'):
    printFuncRun('gencode_biotype_load')
    printFuncArgs()
    biotype_dict = nested_dict()

    with open(txt, 'r') as TXT:
        for line in TXT:
            line = line.rstrip()
            if not line: continue
            if line.startswith('#relation'): 
                arr = line.split()
                biotype_dict['relation']['header'] = arr[1].split(',')
                continue
            if line.startswith('#description'):
                arr = line.split()
                biotype_dict['description']['header'] = arr[1].split(',')
                continue
            if "|" in line: # description
                arr = line.split('|')
                biotype_dict['description'][arr[0]] = arr[1]
            elif "," in line: # relation
                arr = line.split(',')
                biotype_dict[arr[0]]['biotype'] = arr[0]
                biotype_dict[arr[0]]['type_num'] = arr[1]
                biotype_dict[arr[0]]['abstract_biotype'] = arr[2]
                biotype_dict[arr[0]]['type_description'] = biotype_dict['description'][arr[1]]
            else:
                print "un-parsed lines:",line
    for i,j in biotype_dict.items():
        if i != "relation" and i != "description":
            biotype_dict[i]['type_description'] = biotype_dict['description'][j['type_num']]
    print_dict(biotype_dict)
    printFuncRun('gencode_biotype_load')
    return biotype_dict.to_dict()


##################################################
### format conversion
def bigWigToBedGraph(bw,savefn=None):
    printFuncRun('bigWigToBedGraph')
    printFuncArgs()
    if savefn is None:
        savefn = bw.replace('.bw','.bg')
        savefn = savefn.replace('.bigwig','.bg')
    subprocess.call(["bigWigToBedGraph %s %s"%(bw,savefn)],shell=True)
    printFuncRun('bigWigToBedGraph')
    return savefn

def read_liftover_chain():
    human_chain_dir = '/Share/home/zhangqf5/gongjing/paris-2016-05/data/reference_data/GENCODE/human'
    human_chain_ls = ['hg18ToHg19','hg18ToHg38','hg19ToHg18','hg19ToHg38','hg38ToHg19']
    
    mouse_chain_dir = '/Share/home/zhangqf5/gongjing/paris-2016-05/data/reference_data/GENCODE/mouse'
    mouse_chain_ls = ['mm9ToMm10']

    liftover_chain_dict = {}

    for human_chain in human_chain_ls:
        liftover_chain_dict[human_chain] = human_chain_dir+'/'+human_chain+'.over.chain'        
    for mouse_chain in mouse_chain_ls:
        liftover_chain_dict[mouse_chain] = mouse_chain_dir+'/'+mouse_chain+'.over.chain'
        
    return liftover_chain_dict

def liftOver(oldFile,mode='hg18ToHg38',newFile=None):
    printFuncRun('liftOver')
    liftover_chain_dict = read_liftover_chain()
    map_chain = liftover_chain_dict[mode]
    if newFile is None:
        newFile = oldFile+mode.split('To')[1]
    printFuncArgs()
    subprocess.call(["liftOver %s %s %s"%(oldFile,map_chain,newFile)],shell=True)
    printFuncRun('liftOver')
    return newFile

def crossmap(oldFile,mode='hg18ToHg38',newFile=None):
    printFuncRun('crossmap')
    liftover_chain_dict = read_liftover_chain()
    map_chain = liftover_chain_dict[mode]
    if newFile is None:
        newFile = '.'.join(oldFile.split('.')[0:-1])+'.'+mode.split('To')[1]+'.'+oldFile.split('.')[-1]
    printFuncArgs()
    if oldFile.endswith('bed') or oldFile.endswith('bedgraph') or oldFile.endswith('bg'):
        input_fn_type = 'bed'
    else:
        input_fn_type = 'bed'
    subprocess.call(["CrossMap.py %s %s %s %s"%(input_fn_type,map_chain,oldFile,newFile)],shell=True)
    printFuncRun('crossmap')
    return newFile

def bed_uniq_keep_longest_one(bed,col=6):
    printFuncRun('bed_uniq_keep_longest_one')
    printFuncArgs()
    bed_col_ls = []
    bed_col_dir = {}
    num_all = 0
    with open(bed,'r') as BED:
        for line in BED:
            if not line or line.startswith('#'): continue
            num_all += 1
            arr=line.strip().split('\t')
            if bed_col_dir.has_key(arr[col]):
                new_len = int(arr[2])-int(arr[1])
                if new_len > bed_col_dir[arr[col]]:
                    bed_col_dir[arr[col]] = line.strip()
                else:
                    pass
            else:
               bed_col_dir[arr[col]] = line.strip()
               bed_col_ls.append(arr[col])
    bed_uniq = bed.replace('bed','uniq.bed')
    with open(bed_uniq,'w') as BED_UNIQ:
        for bed_col in bed_col_ls:
            print >>BED_UNIQ,bed_col_dir[bed_col] 
    print "all bed line: %s"%(num_all)
    print "all bed line(uniq): %s"%(len(bed_col_ls))
    printFuncRun('bed_uniq_keep_longest_one')
    return bed_uniq


def read_bed(bed):
    printFuncRun('read_bed')
    printFuncArgs()
    bed_dict = nested_dict()
    cols = ['chr','start','end','id','score','strand']
    with open(bed,'r') as BED:
        for line in BED:
            line = line.strip()
            if not line or line.startswith('#'): continue
            chrom,start,end,id,score,strand = line.split('\t')
            for i,j in zip(cols,[chrom,start,end,id,score,strand]):
                bed_dict[id][i] = j 
    print_dict_item0(bed_dict)
    printFuncRun('read_bed')
    return bed_dict



##################################################
### vcf related
def vcf_pattern_comparison(sample_ls,vcf_ls):
    printFuncRun('vcf_pattern_comparison')
    printFuncArgs()

    project_pwd,scripts_pwd=get_project_scripts_pwd()
    sample_vcf=['result/'+i+'/'+j for i in sample_ls for j in vcf_ls]
    sample_vcf_path=[add_file_pwd(vcf) for vcf in sample_vcf]
    sample_vcf_path_info={}
    for vcf in sample_vcf_path:
        arr=vcf.split("/")
        sample=arr[-3]
        caller=''
        if 'GATK' in vcf: 
	    caller='GATK'
        elif 'VarScan' in vcf: 
	    caller='VarScan'
        else:
            caller='samtools'
        sample_vcf_path_info[vcf]=[sample,caller]

    ref_fa_dict=read_fa_into_dict()

    callers=['GATK','samtools','VarScan']
    nrow=len(sample_ls)
    ncol=len(vcf_ls)
    ax_rows=nrow+1+1+1+1+1+1 #
    ax_cols=ncol+1+1
    fig,axs=plt.subplots(ax_rows,ax_cols,figsize=(8*ax_cols,8*ax_rows))

    sample_all_df=pd.DataFrame(columns=['ref','alt','count','sample','caller'])

    sample_all_df_list_dict={} # {sample:[info]}
    sample_all_df_all_cols={} # {sample:df}
    sample_all_df_list_dict_info={} 
    sample_all_shape0_dict=collections.defaultdict(list) # {sample:[shape0...]}
    sample_all_overlap_shape0_dict=collections.defaultdict(list) # {sample:[shape0...]}
    sample_all_for_box_dict={}

    #for n,(vcf,ax) in enumerate(zip(sample_vcf_path,axs.flat[0:nrow*ncol])):
    for n,vcf in enumerate(sample_vcf_path):
        print "sample & caller: ",sample_vcf_path_info[vcf]

        ax_nrow=n/ncol
        if sample_vcf_path_info[vcf][1] == "GATK":
            ax_ncol=0
        elif sample_vcf_path_info[vcf][1] == "samtools":
            ax_ncol=1
        elif sample_vcf_path_info[vcf][1] == "VarScan":
            ax_ncol=2
        else:
            print "[error] unknow caller: %s"%(sample_vcf_path_info[vcf][1])
        ax=axs[ax_nrow][ax_ncol]
        print "ax: [%s,%s]"%(ax_nrow,ax_ncol)

        bed=rreplace(vcf,'vcf','bed')
        if os.path.isfile(bed):
            print "[existed] %s"%(bed)
        else:
            print "[Not exist,will convert] %s"%(vcf)
            subprocess.call(["vcf2bed < %s > %s"%(vcf,bed)],shell=True)
        bed_dbsnp=rreplace(bed,'bed','dbsnp.bed')
        dbsnp_columns=['dbsnp_chr','dbsnp_start','dbsnp_end','dbsnp_id','dbsnp_val','dbsnp_strand']
        if not os.path.isfile(bed_dbsnp):
            bed_dbsnp=bedAnnotationBySnpBed(input_bed=bed,origanism='mouse',source='dbsnp',mode="loj",intersect_bed="")
        else:
            print "[existed] %s"%(bed_dbsnp)
        bed_dbsnp_ucsc=rreplace(bed_dbsnp,'bed','ucsc.bed')
        ucsc_columns=['ucsc_chr','ucsc_start','ucsc_end','ucsc_strand','ucsc_id','ucsc_ref','ucsc_alt']
        if not os.path.isfile(bed_dbsnp_ucsc):
            bed_dbsnp_ucsc=bedAnnotationBySnpBed(input_bed=bed_dbsnp,origanism='mouse',source='ucsc',mode="loj",intersect_bed="")
	else:
            print "[existed] %s"%(bed_dbsnp_ucsc)

        bed=bed_dbsnp_ucsc
        df=pd.read_csv(bed,header=None,sep="\t")
        vcf_bed_cols=['chr','start','end','none','quality','ref','alt','status','info','genotype','genotype_val']
	vcf_bed_cols.extend(dbsnp_columns)
        vcf_bed_cols.extend(ucsc_columns)
        df.columns=vcf_bed_cols

        df['dbsnp_status']=df['dbsnp_chr'].map(lambda x:int(x!='.'))
        df['ucsc_status']=df['ucsc_chr'].map(lambda x:int(x!='.'))
        df['snp_status']=df.apply(lambda row:row['dbsnp_status']+row['ucsc_status'],axis=1)
        df=df.loc[df['ref'].map(lambda x:len(x)==1) & df['alt'].map(lambda x:len(x)==1),:] 
        df_snp=df[df['snp_status'] > 0].shape[0] # entry has >=1 snp
        df=df[df['snp_status'] == 0] # entry has 0 snp
        print df.info()

        sample_all_shape0_dict[sample_vcf_path_info[vcf][0]].append(str(df.shape[0])+'\n'+'('+str(df_snp)+')') # key: sample val:[num_no_snp,num_snp]

        df_list=df[['chr','start','end','ref','alt','snp_status']].to_csv(None,header=False,index=False,sep=' ').split("\n")[0:-1]
        sample_all_df_list_dict[' '.join(sample_vcf_path_info[vcf])]=df_list # key: sample+caller val: mutsite info
        sample_all_df_all_cols[' '.join(sample_vcf_path_info[vcf])]=df.copy()
        #sns.countplot(x='ref',hue='alt',data=df,ax=ax)
        count2=df.groupby(['ref','alt']).count().unstack().fillna(0)
        legend_status = 1 if n==0 else 0
        print "legend_status:",legend_status
        #count2['chr'].plot(kind='bar',stacked=True,ax=ax,title=' '.join(sample_vcf_path_info[vcf])+' '+str(df.shape[0]),legend=legend_status)
        count2['chr'].plot(kind='bar',stacked=True,ax=ax,title=' '.join(sample_vcf_path_info[vcf])+' '+str(df.shape[0])+'('+str(df_snp)+')',legend=legend_status)
        
        if legend_status:
            ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),ncol=4)
            handles,labels=ax.get_legend_handles_labels() # store the label

        count=df.groupby(['ref','alt']).size() # for converient to convert to data frame
        print "count group by ref+alt:"
        print count
        count_df=pd.DataFrame({'count' :count}).reset_index()
        count_df['sample']=sample_vcf_path_info[vcf][0] # add sample label
        count_df['caller']=sample_vcf_path_info[vcf][1] # add caller label
        count_df.sort_values(by=['ref','alt'],axis=0,inplace=True)
        print count_df
        sample_all_df=sample_all_df.append(count_df)
        #count.plot(kind='bar',ax=ax)

    lambda1=lambda x:x['A']['T']+x['T']['A']+x['G']['T']+x['T']['G']

    for n,caller in enumerate(['GATK','samtools','VarScan']):
        sample=sample_ls[2:]
        sample0=sample[0]+' '+caller
        sample1=sample[1]+' '+caller
        sample2=sample[2]+' '+caller
        venn3_ax=axs[ax_nrow+1][n]
        venn3([set(sample_all_df_list_dict[sample0]),set(sample_all_df_list_dict[sample1]),set(sample_all_df_list_dict[sample2])],(sample0.split()[0],sample1.split()[0],sample2.split()[0]),ax=venn3_ax)
        #overlap3=set(sample_all_df_list_dict[sample0]) & set(sample_all_df_list_dict[sample1]) & set(sample_all_df_list_dict[sample2]) # in 3 sets

        overlap3=[j for i in [sample_all_df_list_dict[sample0],sample_all_df_list_dict[sample1],sample_all_df_list_dict[sample2]] for j in i] 
        ol_n=3
        overlap3=[i for i,j in collections.Counter(overlap3).items() if j >= ol_n] # in >=2 sets
        overlap3_ls=[] # [[],[],,,,] 
        for i in overlap3: overlap3_ls.append(i.split(' '))
        overlap3_df=pd.DataFrame(overlap3_ls,columns=['chr','start','end','ref','alt','snp_status'])
        overlap3_df_count=overlap3_df.groupby(['ref','alt']).count().unstack().fillna(0)
        print "[overlap3_df_count] count group by ref+alt: with ol_n=%s"%(ol_n)
        print overlap3_df_count
        overlap3_df_ax=axs[ax_nrow+2][n]
        overlap3_df_count['chr'].plot(kind='bar',stacked=True,ax=overlap3_df_ax,title='overlap: '+caller+' '+str(overlap3_df.shape[0]),legend=legend_status)
       

        """
        if caller == "GATK": # get overlap mutsite DP/AF/AC
                df_tm=pd.merge()
		overlap3_df['DP']=overlap3_df['info'].map(lambda x:float({i.split('=')[0]:i.split('=')[1] for i in x.split(";")}['DP']))
		overlap3_df['AF']=overlap3_df['info'].map(lambda x:float({i.split('=')[0]:i.split('=')[1] for i in x.split(";")}['AF']))
		overlap3_df['AC']=overlap3_df['info'].map(lambda x:float({i.split('=')[0]:i.split('=')[1] for i in x.split(";")}['AC']))
        	overlap3_df.plot(kind='scatter',x='DP',y='AF',c='AC',ax=axs[ax_nrow+2][n+1],title='overlap: '+caller+' '+str(overlap3_df.shape[0])+' '+'allele depth vs frequency',legend=1)
        """
        #overlap3_df_count_snp=overlap3_df.groupby(['ref','snp_status']).count().unstack().fillna(0) # plot snp
        #overlap3_df_count_snp['chr'].plot(kind='bar',stacked=True,ax=axs[ax_nrow+3][n],title='overlap: '+caller+' with snp '+str(overlap3_df.shape[0]),legend=1,colormap='jet')

        ax_mutsite_l_r=axs[ax_nrow+4][n]
        neighbor_l=1
        neighbor_r=-1
        neighbor=5
        overlap3_df['seq']=overlap3_df.apply(lambda x:ref_fa_dict[x['chr']][(int(x['start'])-neighbor):(int(x['start'])+neighbor+1)],axis=1)
        #overlap3_df['seq']=overlap3_df.apply(get_row_seq,ref_fa_dict=ref_fa_dict,neighbor=neighbor,axis=1)
        neighbor=1
        overlap3_df['seq_1']=overlap3_df.apply(lambda x:ref_fa_dict[x['chr']][(int(x['start'])-neighbor):(int(x['start'])+neighbor+1)],axis=1)
        #overlap3_df['seq_1']=overlap3_df.apply(get_row_seq,ref_fa_dict=ref_fa_dict,neighbor=neighbor,axis=1)
        """
        overlap3_df['seq_l']=overlap3_df.apply(lambda x:ref_fa_dict[x['chr']][(int(x['start'])-neighbor_l):(int(x['start'])+neighbor_r+1)],axis=1)
        overlap3_df_count=overlap3_df.groupby(['ref','seq_l']).count().unstack().fillna(0)
        overlap3_df_count=overlap3_df_count.applymap(lambda x:-x)
        overlap3_df_count_per=overlap3_df_count['chr'].apply(lambda x:[-float(i)/sum(x) for i in x],axis=1) # plot bar with percentage
        #overlap3_df_count_per.plot(kind='barh',stacked=1,ax=ax_mutsite_l_r,title='overlap seq around MutSite: '+caller+' '+str(overlap3_df.shape[0]),legend=0)
        ax_add1=ax_add_sub_ax(fig=None,ax=ax_mutsite_l_r,width="50%",height="50%",loc=1,axis_bgcol='white',xticks=None,yticks=None)
        count_dict=overlap3_df_count['chr'].to_dict()
        AU=map(lambda1,[count_dict])[0]
        #plt.text(0.5, 0.95, 'MutSite left 1 base'+'\n'+str(AU)+' %.2f'%(AU/float(overlap3_df.shape[0])), fontsize=20, bbox=dict(facecolor='white', alpha=0.7),transform=ax_mutsite_l_r.transAxes,horizontalalignment='center',verticalalignment='center')
        print "[overlap3_df_count] left 1"
        print overlap3_df_count['chr']
        print overlap3_df_count_per
        print count_dict
        neighbor_l=-1
        neighbor_r=1
        overlap3_df['seq_r']=overlap3_df.apply(lambda x:ref_fa_dict[x['chr']][(int(x['start'])-neighbor_l):(int(x['start'])+neighbor_r+1)],axis=1)
        overlap3_df_count=overlap3_df.groupby(['ref','seq_r']).count().unstack().fillna(0)
        #overlap3_df_count=overlap3_df_count.applymap(lambda x:-x)
        overlap3_df_count_per=overlap3_df_count['chr'].apply(lambda x:[float(i)/sum(x) for i in x],axis=1)
        overlap3_df_count_per.plot(kind='barh',stacked=1,ax=ax_mutsite_l_r,legend=0)
        count_dict=overlap3_df_count['chr'].to_dict()
        AU=abs(map(lambda1,[count_dict])[0])
        #plt.text(0.5, 0.05, 'MutSite right 1 base'+'\n'+str(AU)+' %.2f'%(AU/float(overlap3_df.shape[0])), fontsize=20, bbox=dict(facecolor='white', alpha=0.7),transform=ax_mutsite_l_r.transAxes,horizontalalignment='center',verticalalignment='center')
        print "[overlap3_df_count] right 1"
        print overlap3_df_count['chr']
        print overlap3_df_count_per
        print count_dict
        #ax_mutsite_l_r.hlines(y=0,xmin=-100,xmax=100,linestyles='--',colors='black') 
        ax_mutsite_l_r.vlines(x=0,ymin=-100,ymax=100,linestyles='--',colors='black') 

        print overlap3_df['seq'][0:5]
        print overlap3_df['seq_l'][0:5]
        print overlap3_df['seq_r'][0:5]
        print overlap3_df['ref'][0:5]
        """
        for nn,N in enumerate(['C','A','T','G']):
            seq_ls=overlap3_df[overlap3_df['ref']==N]['seq_1']
            fa=list2fa(seq_ls=seq_ls,savefn='result/plots/overlap'+str(ol_n)+'.'+caller+'.'+N+'.fa')
            weblogo_png=plotFaWeblogo(input_fa=fa,weblogo_f_title="")
            ax_add=ax_add_sub_ax(fig=None,ax=ax_mutsite_l_r,width="50%",height="50%",loc=nn+1,axis_bgcol='white',xticks=None,yticks=None)
            ax_show_img(fig=None,ax=ax_add,img_fn=weblogo_png)

        ax_mutsite_l_r_n=axs[ax_nrow+3][n]
        
        #ls_weblogo_plot(ls_ls=overlap3_df['seq'],pos_labels=range(-neighbor,neighbor+1),ax=ax_mutsite_l_r_n,title='overlap seq content: '+caller+' ')
        fa=list2fa(seq_ls=overlap3_df['seq'],savefn='result/plots/overlap3.fa')
        weblogo_png=plotFaWeblogo(input_fa=fa,weblogo_f_title="")        
        ax_show_img(fig=None,ax=ax_mutsite_l_r_n,img_fn=weblogo_png)


        if n >= 0:
                for m,sample01 in enumerate(sample_ls[0:2]):
			sample99=sample01+' '+caller
			sample99_ls=[]
			for i in sample_all_df_list_dict[sample99]: sample99_ls.append(i.split(' '))
			sample99_df=pd.DataFrame(sample99_ls,columns=['chr','start','end','ref','alt','snp_status'])
			#sample99_df_count=sample99_df.groupby(['ref','alt']).count().unstack().fillna(0)
			#print sample99_df_count
			#sample99_df_count['chr'].plot(kind='bar',stacked=True,ax=axs[ax_nrow+5][n],title='sample99: '+caller+' '+str(sample99_df.shape[0]),legend=0)
                        ax_mutsite_append=axs[ax_nrow+5+m][n]
			neighbor_l=1
			neighbor_r=-1
			neighbor=5
			sample99_df['seq']=sample99_df.apply(lambda x:ref_fa_dict[x['chr']][(int(x['start'])-neighbor):(int(x['start'])+neighbor+1)],axis=1)
                        #sample99_df['seq']=sample99_df.apply(get_row_seq,ref_fa_dict=ref_fa_dict,neighbor=neighbor,axis=1) # vcf has no strand info
                        neighbor=1
			sample99_df['seq_1']=sample99_df.apply(lambda x:ref_fa_dict[x['chr']][(int(x['start'])-neighbor):(int(x['start'])+neighbor+1)],axis=1)
                        #sample99_df['seq_1']=sample99_df.apply(get_row_seq,ref_fa_dict=ref_fa_dict,neighbor=neighbor,axis=1)
			for nn,N in enumerate(['C','A','T','G']):
			    seq_ls=sample99_df[sample99_df['ref']==N]['seq_1']
			    fa=list2fa(seq_ls=seq_ls,savefn='result/plots/'+sample99.split()[0]+'.'+caller+'.'+N+'.fa')
			    weblogo_png=plotFaWeblogo(input_fa=fa,weblogo_f_title="")
			    ax_add=ax_add_sub_ax(fig=None,ax=ax_mutsite_append,width="50%",height="50%",loc=nn+1,axis_bgcol='white',xticks=None,yticks=None)
			    ax_show_img(fig=None,ax=ax_add,img_fn=weblogo_png)

                        if caller == "GATK":
				ax_mutsite_append_hist=axs[ax_nrow+5+m][n+3]
				df_tmp=sample_all_df_all_cols[sample99]
				df_tmp['DP']=df_tmp['info'].map(lambda x:float({i.split('=')[0]:i.split('=')[1] for i in x.split(";")}['DP']))
				df_tmp['DP'].plot(kind='box',ax=ax_mutsite_append_hist,title='DP density',legend=sample99+' '+str(df_tmp.shape[0]),vert=False)
				ax_mutsite_append_hist.set_xlim(0,25)

                                sample_all_for_box_dict[sample99]=df_tmp['DP']
                        """
			sample99_df['seq_l']=sample99_df.apply(lambda x:ref_fa_dict[x['chr']][(int(x['start'])-neighbor_l):(int(x['start'])+neighbor_r+1)],axis=1)
			sample99_df_count=sample99_df.groupby(['ref','seq_l']).count().unstack().fillna(0)
			sample99_df_count=sample99_df_count.applymap(lambda x:-x)
			sample99_df_count_per=sample99_df_count['chr'].apply(lambda x:[-float(i)/sum(x) for i in x],axis=1)
			sample99_df_count_per.plot(kind='barh',stacked=1,ax=ax_mutsite_append,title=sample01+' seq around MutSite: '+caller+' '+str(sample99_df.shape[0]),legend=0)
                        #ax_add_sub_ax(fig=None,ax=ax_mutsite_append,width="50%",height="50%",loc=1,axis_bgcol='white',xticks=None,yticks=None)
                        count_dict=sample99_df_count['chr'].to_dict()
                        AU=map(lambda1,[count_dict])[0]
			#plt.text(0.5, 0.95, 'MutSite left 1 base'+'\n'+str(AU)+' %.2f'%(AU/float(sample99_df.shape[0])), fontsize=20, bbox=dict(facecolor='white', alpha=0.7),transform=ax_mutsite_append.transAxes,horizontalalignment='center',verticalalignment='center')
			print "[sample99_df_count] left 1"
                        print sample99_df_count['chr']
			print count_dict
			neighbor_l=-1
			neighbor_r=1
			sample99_df['seq_r']=sample99_df.apply(lambda x:ref_fa_dict[x['chr']][(int(x['start'])-neighbor_l):(int(x['start'])+neighbor_r+1)],axis=1)
			sample99_df_count=sample99_df.groupby(['ref','seq_r']).count().unstack().fillna(0)
			#sample99_df_count=sample99_df_count.applymap(lambda x:-x)
			sample99_df_count_per=sample99_df_count['chr'].apply(lambda x:[float(i)/sum(x) for i in x],axis=1)
			sample99_df_count_per.plot(kind='barh',stacked=1,ax=ax_mutsite_append,legend=0)
                        count_dict=sample99_df_count['chr'].to_dict()
                        AU=abs(map(lambda1,[count_dict])[0])
			#plt.text(0.5, 0.05, 'MutSite right 1 base'+'\n'+str(AU)+' %.2f'%(AU/float(sample99_df.shape[0])), fontsize=20, bbox=dict(facecolor='white', alpha=0.7),transform=ax_mutsite_append.transAxes,horizontalalignment='center',verticalalignment='center')
			print "[sample99_df_count] right 1"
                        print sample99_df_count['chr']
			print count_dict
			#ax_mutsite_append.hlines(y=0,xmin=-100,xmax=100,linestyles='--',colors='black') 
			ax_mutsite_append.vlines(x=0,ymin=-100,ymax=100,linestyles='--',colors='black') 
                        """

	if caller == "GATK":
		ax_mutsite_append_hist=axs[ax_nrow+3][n+3]
                for sample_tmp in sample_ls:
                    df_tmp=sample_all_df_all_cols[sample_tmp+' '+caller]
		    df_tmp['DP']=df_tmp['info'].map(lambda x:float({i.split('=')[0]:i.split('=')[1] for i in x.split(";")}['DP']))
	            #df_tmp['DP'].plot(kind='box',ax=ax_mutsite_append_hist,title='DP density',legend=sample99+' '+str(df_tmp.shape[0]),vert=False)
                    sample_all_for_box_dict[sample_tmp+' '+caller]=df_tmp['DP']
                df_tmp_box=pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in sample_all_for_box_dict.iteritems() ]))
                df_tmp_box.plot(kind='box',ax=ax_mutsite_append_hist,vert=False,title='depth at called mutation site')
	        ax_mutsite_append_hist.set_xlim(0,25)

        overlap3_bk=[j for i in [sample_all_df_list_dict[sample0],sample_all_df_list_dict[sample1],sample_all_df_list_dict[sample2]] for j in i] 
        for k in xrange(len(sample)):
            n=len([i for i,j in collections.Counter(overlap3_bk).items() if j == k+1])
            sample_all_overlap_shape0_dict[caller].append(n) # key: caller val: [1set,2set,3set]
        print sample_all_overlap_shape0_dict

    for n,sample in enumerate(sample_ls):
        call0=sample+' '+'GATK'
        call1=sample+' '+'samtools'
        call2=sample+' '+'VarScan'
        venn3_ax_sample=axs[n][ncol]
        venn3([set(sample_all_df_list_dict[call0]),set(sample_all_df_list_dict[call1]),set(sample_all_df_list_dict[call2])],('GATK','samtools','VarScan'),ax=venn3_ax_sample)

        overlap3=set(sample_all_df_list_dict[call0]) & set(sample_all_df_list_dict[call1]) & set(sample_all_df_list_dict[call2]) 
        overlap3_ls=[]
        for i in overlap3: overlap3_ls.append(i.split(' '))
        overlap3_df=pd.DataFrame(overlap3_ls,columns=['chr','start','end','ref','alt','snp_status'])
        overlap3_df_count=overlap3_df.groupby(['ref','alt']).count().unstack().fillna(0)
        ax_overlap_caller=axs[n][ncol+1]
        overlap3_df_count['chr'].plot(kind='bar',stacked=True,ax=ax_overlap_caller,title='overlap: '+sample+' '+str(overlap3_df.shape[0]),legend=legend_status)
    axs[-1, -1].axis('off')
    axs[-1,-1].legend(handles,labels,loc='center',fontsize='x-large') # reset plot labels on last axis

    sample_all_shape0_dict_df=pd.DataFrame(sample_all_shape0_dict,index=callers).T
    sample_all_overlap_shape0_dict_df=pd.DataFrame(sample_all_overlap_shape0_dict,index=['1 sample','2 sample','3 sample'])
    sample_all_overlap_shape0_dict_df=sample_all_overlap_shape0_dict_df[callers] 
    print sample_all_shape0_dict_df
    print sample_all_overlap_shape0_dict_df

    ax_df_static=axs[nrow,ncol]
    df_to_table_png(df=pd.concat([sample_all_shape0_dict_df,sample_all_overlap_shape0_dict_df],axis=0),fig=None,ax=ax_df_static,colLabels=callers,rowLabels=1,cellLoc='center',colwidth=0.17,savefn=None,scale_x=1.5,scale_y=4.2,fontsize=20)

    savefn='result/plots/'+'_'.join(sample_ls)+'_'+'GATK_samtools_varscan'+'.ge'+str(ol_n)+'.png'
    savefn=add_file_pwd(savefn)

    fig.tight_layout()
    fig.savefig(savefn)
    plt.close()

    print sample_all_df
    sample_all_df.to_csv(rreplace(savefn,'png','csv'),sep='\t',index=False)

    fig,ax=plt.subplots()
    #df_sns_facetgrid(df="",col="",col_wrap="",row="",hue="",plot_col="",plot_type="plt.hist",sharex=False,sharey=True,size=4,aspect=1,palette="",margin_titles=1,hue_order="",row_order="",col_order="",savefn="",titles_ls="",suptitle_str="")
    g=sns.FacetGrid(data=sample_all_df,row='caller',col='sample',hue='alt',size=6,margin_titles=True,sharey=False)
    g=(g.map(sns.barplot,'ref','count').add_legend())
    g.savefig(rreplace(savefn,'png','2.png')) 

    printFuncRun('vcf_pattern_comparison')
    print
    return

def get_row_seq(x,ref_fa_dict,neighbor=1):
    seq=ref_fa_dict[x['chr']][(int(x['start'])-neighbor):(int(x['start'])+neighbor+1)]
    base_comp={'A':'T','T':'A','C':'G','G':'C'}
    if x['strand'] == '-':
        seq=''.join([base_comp[i] for i in seq])[::-1]
    return seq

def read_hgnc_id_conversion(hgnc_txt='/Share/home/zhangqf/gongjing/paris-2016-05/data/HGNC/hgnc_complete_set.txt',key='symbol'):
    printFuncRun('read_hgnc_id_conversion')
    printFuncArgs()
    hgnc_id_conversion_dict = nested_dict()
    with open(hgnc_txt,'r') as HGNC_TXT:
        line_n = 0
        for line in HGNC_TXT:
            if not line or line.startswith('#'): continue
            if line.startswith('hgnc_id'):
                header_ls = line.strip('\n').split('\t')
                if key not in header_ls:
                    print "provide key not in headers: %s"%(key)
                    print "avalable headers:",header_ls
                    return
                key_idx = header_ls.index(key)
                print "  - %s idx: %s"%(key,key_idx)
                continue
          
            arr = line.strip('\n').split('\t')
            if len(arr) != len(header_ls):
                print line,arr
            key_val = arr[key_idx]
            for header,val in zip(header_ls,arr):
                hgnc_id_conversion_dict[key_val][header] = val
            line_n += 1
    print "read lines: %s, hgnc_id_conversion_dict: %s"%(line_n,len(hgnc_id_conversion_dict))
    print_dict_item0(hgnc_id_conversion_dict)
    printFuncRun('read_hgnc_id_conversion')
    return hgnc_id_conversion_dict

def read_gene2ensembl(gene2ensembl='/Share/home/zhangqf/gongjing/software/goatools-master/gene2ensembl',tax='9606',key='RNA_nucleotide_accession.version'):
    printFuncRun('read_gene2ensembl')
    printFuncArgs()
    gene2ensembl_dict = nested_dict()
    header_ls =['tax_id','GeneID','Ensembl_gene_identifier','RNA_nucleotide_accession.version','Ensembl_rna_identifier','protein_accession.version','Ensembl_protein_identifier']
    key_idx = header_ls.index(key)
    with open(gene2ensembl,'r') as TXT:
        line_n = 0
        for line in TXT:
            if not line or line.startswith('#'): continue
            arr = line.strip().split('\t')
            arr = [i.split('.')[0] for i in arr]
            tax_id = arr[0]
            key_val = arr[key_idx]
            for header,val in zip(header_ls,arr):
                gene2ensembl_dict[tax_id][key_val][header] = val
            line_n += 1
    print "line num: %s"%(line_n)
    tax_gene2ensembl_dict = gene2ensembl_dict[tax]
    print "%s entry: %s"%(tax,len(tax_gene2ensembl_dict))
    print_dict_item0(tax_gene2ensembl_dict)        
    printFuncRun('read_gene2ensembl')
    return tax_gene2ensembl_dict


def read_mygene(my_gene):
    mygene = nested_dict()
    
    with open(my_gene, 'r') as MYGENE:
        for line in MYGENE:
            line = line.strip("\n")
            if not line or line.startswith(('#','query')): continue
            gene_name, _id, _score, alias, ensembl, entrezgene, notfound, symbol = line.split('\t')
            if notfound == "True": continue
            if ensembl == "": continue
            #print line
            gene = json.loads(ensembl.replace('u','').replace('[','').replace(']','').split(',')[0].replace("'",'"'))
            if mygene.has_key(gene_name):
                _score = float(_score)
                if _score > mygene[gene_name]['score']:
                    mygene[gene_name]['score'] = _score
                    mygene[gene_name]['alias'] = _score
                    mygene[gene_name]['Gene ID'] = _score
                    mygene[gene_name]['notfound'] = _score
                    mygene[gene_name]['symbol'] = _score
            else:
                mygene[gene_name]['score'] = _score
                mygene[gene_name]['alias'] = alias
                mygene[gene_name]['Gene ID'] = gene['gene']
                mygene[gene_name]['notfound'] = notfound
                mygene[gene_name]['symbol'] = symbol
    print_dict_item0(mygene)
    return mygene


##################################################
### gtf file related
# read .gtf into dict, gtf_dict[feature][feature_id][desciption] , gtf_dict['gene']['ENSG0001']={'start':2,'end':9,....,'attr':{'gene_name':p53,'gene_type':protein_coding,...,},....,'transcript_id':[ENST001,ENST002]}
def read_gtf(gtf='/Share/home/zhangqf/gongjing/paris-2016-05/result/sunlei/data/index/human/genome/genome.gtf', savejson=0):
    printFuncRun('read_gtf')
    printFuncArgs()

    gtf_json = gtf+'.json'
    if os.path.isfile(gtf_json):
        printFuncRun('read from json: %s'%(gtf_json))
        json_in = open(gtf_json, 'rb')
        gtf_dict = json.load(json_in)
        json_in.close()
        printFuncRun('read from json: %s'%(gtf_json))
        return gtf_dict

    gtf_dict=defaultdict(dict)
    #gtf_dict = nested_dict()
    gtf_cols=['seq_name','source','feature','start','end','score','strand','frame','attr']
    tx_longest_dict = {}
    with open(gtf,'r') as GTF:
        n_line = 0
        for n,line in enumerate(GTF):
            if n%300000 == 0: print "read lines: %s"%(n),timestamp()
            if not line or line.startswith('#'): continue
            n_line += 1 
            arr = line.strip().split('\t')
            try:
                attr_dict = {i.strip().split(' ')[0]:i.strip().split(' ')[1] for i in arr[-1].replace('"','').split(';')[0:-1]}
                for i,j in attr_dict.items():
                    if i != "gene_name":
                        attr_dict[i] = j.split('.')[0]
            except:
                print line
            entry_id = attr_dict[arr[2]+'_id'] if attr_dict.has_key(arr[2]+'_id') else attr_dict['gene_id']
            if not gtf_dict[arr[2]].has_key(entry_id):
                gtf_dict[arr[2]][entry_id] = {'transcript_id':[]}
            for i in xrange(8):
                gtf_dict[arr[2]][entry_id][gtf_cols[i]] = arr[i]
            
            gtf_dict[arr[2]][entry_id]['attr'] = attr_dict
            try:
                gtf_dict[arr[2]][entry_id]['transcript_id'].append(attr_dict['transcript_id'])
            except:
                pass
            if arr[2] == "transcript":
                tx_len = int(arr[4]) - int(arr[3])
                if not tx_longest_dict.has_key(attr_dict['gene_id']):
                    tx_longest_dict[attr_dict['gene_id']] = [attr_dict['transcript_id'],tx_len]
                else:
                    if tx_len > tx_longest_dict[attr_dict['gene_id']][1]:
                        tx_longest_dict[attr_dict['gene_id']] = [attr_dict['transcript_id'],tx_len]
    for transcript_id,transcript_val in gtf_dict['transcript'].items():
        gene_id = transcript_val['attr']['gene_id']
        gtf_dict['gene'][gene_id]['transcript_id'].append(transcript_id)
    print "total entry lines: %s"%(n_line)
    feature_ls = gtf_dict.keys()
    feature_ls_num = 0
    feature_ls_transcript_num = 0
    feature_ls_uniq_transcript_num = 0
    for feature in feature_ls:  
        feature_num = len(gtf_dict[feature])
        feature_ls_num += feature_num
        feature_transcript_num = sum([len(j['transcript_id']) for i,j in gtf_dict[feature].items()])
        feature_uniq_transcript_num = sum([len(set(j['transcript_id'])) for i,j in gtf_dict[feature].items()])
        feature_ls_transcript_num += feature_transcript_num
        feature_ls_uniq_transcript_num += feature_uniq_transcript_num
        print "%s: %s"%(feature,feature_num),
    print "total: %s"%(feature_ls_num) 
    print "total(region based): %s"%(feature_ls_transcript_num)
    print "total(transcript based): %s"%(feature_ls_uniq_transcript_num)
    gtf_dict['tx_longest'] = {}
    for i,j in tx_longest_dict.items():
        tx = j[0]
        gtf_dict['tx_longest'][tx] = gtf_dict['transcript'][tx]
    feature_ls = gtf_dict.keys()
    for feature in feature_ls:
        print feature, feature_ls
        if gtf_dict[feature] == {}:
            print "empty"
            continue
        feature_idx0 = gtf_dict[feature].keys()[0]
        print feature,feature_idx0,gtf_dict[feature][feature_idx0]
    #for i,j in gtf_dict['gene'].items():
    #    gtf_dict['gene'][i.split('.')[0]] = j
    #for i,j in gtf_dict['transcript'].items():
    #    gtf_dict['transcript'][i.split('.')[0]] = j
    #    gtf_dict['transcript'][i.split('.')[0]]['len'] = int(j['end']) - int(j['start'])
    printFuncRun('read_gtf')

    print
    if savejson:
        printFuncRun('save gtf_dict to json')
        json_fn = gtf+'.json'
        json_out = open(json_fn,'wb')
        json.dump(gtf_dict, json_out)
        json_out.close()
        printFuncRun('save gtf_dict to json')
    return gtf_dict

def read_gtf_general(gtf='/Share/home/zhangqf5/gongjing/DNA-RNA-Protein-interaction-correlation-12-18/data/id_conversion/miRBase/mmu.gff3', attr_sep='='):
    printFuncRun('read_gtf_general')
    printFuncArgs()
    gtf_dict = nested_dict()
    gtf_cols=['seq_name','source','feature','start','end','score','strand','frame','attr']
    with open(gtf,'r') as GTF:
        for line in GTF:
            line = line.rstrip()
            if not line or line.startswith('#'): continue
            arr = line.strip().split('\t')
            try:
                attr_dict = {i.strip().split(attr_sep)[0]:i.strip().split(attr_sep)[1].split('.')[0] for i in arr[-1].replace('"','').split(';')[0:]}
            except:
                print line
            entry_id = attr_dict['ID']
            for i in xrange(8):
                gtf_dict[arr[2]][entry_id][gtf_cols[i]] = arr[i]
            gtf_dict[arr[2]][entry_id]['attr'] = attr_dict
    feature_ls = gtf_dict.keys()
    for feature in feature_ls:
        print feature
        print_dict_item0(gtf_dict[feature])
                           
    printFuncRun('read_gtf_general')
    return gtf_dict

def read_miRbase_ensembl_txt(txt='/Share/home/zhangqf/gongjing/DNA-RNA-Protein-interaction-correlation-12-18/data/id_conversion/human_mouse_miRbase_ensembl.txt',key='miRBase Accession(s)'):
    printFuncRun('read_miRbase_ensembl_txt')
    printFuncArgs()
    miRbase_ensembl_dict = nested_dict()
    with open(txt,'r') as TXT:
        for line in TXT:
            if not line or line.startswith(('#','Gene ID',"Gene stable ID")):
                header_ls = line.strip().split('\t')
                key_idx = header_ls.index(key)
                continue
            arr = line.strip().split('\t')
            key_val = arr[key_idx]
            for i,j in zip(header_ls,arr):
                miRbase_ensembl_dict[key_val][i] = j
    print_dict_item0(miRbase_ensembl_dict)
    printFuncRun('read_miRbase_emsembl_txt')
    return miRbase_ensembl_dict

def gffutils_create_db(gtf_fn='',db_fn=None,force_bool=True,keep_order_bool=True,merge_strategy_method='merge',sort_attribute_values_bool=True):
    printFuncRun('gffutils_create_db')
    printFuncArgs()
    if db_fn is None:
        print "db_fn is None"
        db_fn = gtf_fn+'.db'
        print "will write to db: %s"%(db_fn)
    db = gffutils.create_db(fn=gtf_fn,dbfn=db_fn,force=force_bool,keep_order=keep_order_bool,merge_strategy=merge_strategy_method, sort_attribute_values=sort_attribute_values_bool)
    printFuncRun('gffutils_create_db')


def read_icshape(icshape_out='/Share/home/zhangqf/gongjing/paris-2016-05/result/sunlei/data/icshape/mES/transcript/mes_cy_vivo_transcript.icshape.out',target_bed=None,target_bed_dict=None,save2json=None):
    printFuncRun('read_icshape')
    printFuncArgs()

    if target_bed:
        icshape_target_dict={}
        with open(target_bed,'r') as BED:
            for line in BED:
                if not line or line.startswith('#'): continue
                arr=line.strip().split('\t')
                icshape_target_dict[arr[0]]={'start':int(arr[1]),'end':int(arr[2])}
        print "icshape_target_dict from %s: %s"%(target_bed,len(icshape_target_dict))
        icshape_target_cal_gini_n = 0

    if target_bed_dict:
        icshape_target_dict=target_bed_dict
        print "icshape_target_dict from %s: %s"%(target_bed,len(icshape_target_dict))
        icshape_target_cal_gini_n = 0

    icshape_dict={}
    icshape_out_json = icshape_out+'.jsonss'
    if os.path.isfile(icshape_out_json):
        print "read from %s"%(icshape_out_json)
        output = open(icshape_out_json,'rb')
        icshape_dict = json.load(output)
        output.close()
    else:
        with open(icshape_out,'r') as ICSHAPE:
            print "read from %s"%(icshape_out)
            for line in ICSHAPE:
                if not line or line.startswith('#'): continue
                arr=line.strip().split('\t')
                transcript,length,rpkm,base_reactivity = arr[0].strip(),arr[1],arr[2],arr[3:]
                icshape_dict[transcript]={'icshape-length':length,'icshape-RPKM':rpkm,'base_reactivity_str':','.join(base_reactivity),'base_reactivity_ls':base_reactivity}
                icshape_dict[transcript]['gini']=gini([float(i) for i in base_reactivity if i != "NULL"])

                if (target_bed or target_bed_dict) and icshape_target_dict.has_key(transcript):
                    icshape_target_dict[transcript]={'icshape-length':length,'icshape-RPKM':rpkm,'base_reactivity_str':','.join(base_reactivity),'base_reactivity_ls':base_reactivity}
                    icshape_target_dict[transcript]['gini']=gini([float(i) for i in base_reactivity[icshape_target_dict[transcript]['start']:icshape_target_dict[transcript]['end']] if i != "NULL"])
                    icshape_target_cal_gini_n += 1

    if target_bed:
        print "icshape_target_cal_gini_n: %s"%(icshape_target_cal_gini_n)
        return icshape_dict,icshape_target_dict,icshape_out
    if target_bed_dict:
        print "icshape_target_cal_gini_n: %s"%(icshape_target_cal_gini_n)
        return icshape_dict,icshape_target_dict,icshape_out

    if not save2json is None:
        print "save icshape_dict to json"
        output = open(icshape_out+'.json', 'wb')
        json.dump(icshape_dict, output)
        output.close()

    printFuncRun('read_icshape')
    return icshape_dict,icshape_out

def gini(list_of_values,mode='gini',null_pct=1):
    if len(list_of_values) == 0: return -1
    if list_of_values.count('NULL')/float(len(list_of_values)) > null_pct: return -1
    list_of_values = [i for i in list_of_values if i != 'NULL']
    if len(list_of_values) == 0: return -1
    if type(list_of_values[0]) is str:
        list_of_values = map(float,list_of_values)
    if mode == 'mean_reactivity':
        return np.mean(list_of_values)
    if mode == 'gini':
        if sum(list_of_values) == 0: return 0.67
        sorted_list = sorted(list_of_values)
        height, area = 0, 0
        for value in sorted_list:
            height += value
            area += height - value / 2.
        fair_area = height * len(list_of_values) / 2.
        return (fair_area - area) / fair_area

##################################################
### public data related
def proteinatlas_coexpression(gene_id1='ENSG00000141510',gene_id2='ENSG00000146648',mode='celline',proteinatlas_dir='/Share/home/zhangqf5/gongjing/DNA-RNA-Protein-interaction-correlation-12-18/data/proteinatlas',df=None,only_corr_val=0):
    printFuncRun('proteinatlas_coexpression')
    printFuncArgs()
    if mode not in ['celline','tissue']:
        print "[error] mode not found: %s"%(mode)
        sys.exit()

    savefn_id = proteinatlas_dir+'/'+'exp_corr'+'/'+gene_id1+'_'+gene_id2+'.corr.png'
    if os.path.isfile(savefn_id):
        print "the query ids exist: %s"%(savefn_id)
        printFuncRun('proteinatlas_coexpression')
        sys.exit()

    tpm_csv = proteinatlas_dir+'/'+'rna_'+mode+'.csv'
    if df is None:
        df = pd.read_csv(tpm_csv,sep=',',header=0)
    
    df_gene1 = df[df['Gene']==gene_id1]
    df_gene2 = df[df['Gene']==gene_id2]
    print "%s: "%(gene_id1),df_gene1.shape
    #print  df_gene1.head()
    print "%s: "%(gene_id2),df_gene2.shape 

    if df_gene1.shape[0] == 0 or df_gene2.shape[0] == 0:
        print "cannot plot due to NOT abailable data"
        sys.exit()

    gene_name1 = df_gene1.iloc[0,1]
    gene_name2 = df_gene2.iloc[0,1]
    df_gene = pd.merge(df_gene1,df_gene2,how='outer',on='Sample')
    df_gene.rename(columns={'Value_x':gene_name1,'Value_y':gene_name2},inplace=True)
    df_gene.replace(0,0.001,inplace=True)
    #print list(df_gene[gene_name1])
    #print list(df_gene[gene_name2])
    df_gene["log2(TPM) of %s"%(gene_name1)] = np.log2(df_gene[gene_name1])
    df_gene["log2(TPM) of %s"%(gene_name2)] = np.log2(df_gene[gene_name2])
    #print df_gene.head()
    
    title = '' if df_gene1.shape[0] != df_gene2.shape[0] else "%s\n%s"%(mode,df_gene1.shape[0])
    
    df_sns_jointplot(col_str_x="log2(TPM) of %s"%(gene_name1),col_str_y="log2(TPM) of %s"%(gene_name2),savefn=savefn_id,df=df_gene,list1='list1',list2='list2',xlim=None,ylim=None,x_y_lim_same=1,title_str=title,title_suptitle='right',use_scale_x_y_lim=1,color=None,xlabel=None,ylabel=None)
    rho, pval = stats.pearsonr(df_gene["log2(TPM) of %s"%(gene_name1)],df_gene["log2(TPM) of %s"%(gene_name2)])

    random_corr_txt = proteinatlas_dir+'/'+'rna_'+mode+'.random_corr.txt'
    df_random = pd.read_csv(random_corr_txt,sep='\t',header=0)
    print df_random.head()
    df_random.dropna(inplace=True)
    print df_random.head(),df_random.shape
    #df_sns_jointplot(col_str_x='Pearson',col_str_y='p_val',savefn=savefn_id.replace('png','random.png'),df=df_random,x_y_lim_same=0,title_str='') 
    fig, ax = plt.subplots()
    df_random.plot.scatter(x='Pearson',y='p_val',ax=ax)
    ax.scatter(rho,pval,s=50,color='red')
    plt.savefig(savefn_id.replace('png','random.png'))
    plt.close()

    printFuncRun('proteinatlas_coexpression')

def proteinatlas_random_corr(mode='celline',proteinatlas_dir='/Share/home/zhangqf5/gongjing/DNA-RNA-Protein-interaction-correlation-12-18/data/proteinatlas',shuffle_n=100000):
    printFuncRun('proteinatlas_random_corr')
    printFuncArgs()
    if mode not in ['celline','tissue']:
        print "[error] mode not found: %s"%(mode)
        sys.exit()
        
    tpm_csv = proteinatlas_dir+'/'+'rna_'+mode+'.csv'
    df = pd.read_csv(tpm_csv,sep=',',header=0)
    gene_ls = list(set(list(df['Gene'])))
    print "gene_ls: %s"%(len(gene_ls))
    np.random.seed(1000)

    gene_corr_pair_ls = [list(np.random.choice(gene_ls,2)) for i in xrange(shuffle_n)]
    """
    corr_val_ls = []
    p_val_ls = []
    for n,gene_pair in enumerate(gene_corr_pair_ls):
        print n,gene_pair
        if n % 5000 == 0: print "process: %s"%(n)
        df_gene1 = df[df['Gene']==gene_pair[0]] 
        df_gene2 = df[df['Gene']==gene_pair[1]]
        df_gene = pd.merge(df_gene1,df_gene2,how='outer',on='Sample').replace(0,0.001)
        rho, pval = stats.spearmanr(df_gene['Value_x'],df_gene['Value_y']) 
        corr_val_ls.append(rho)
        p_val_ls.append(pval)
        print >>TXT,'\t'.join([gene_pair[0],gene_pair[1],str(rho),str(pval)])
    df_sns_jointplot(col_str_x=None,col_str_y=None,savefn=tpm_csv.replace('csv','random_corr.png'),df=None,list1=corr_val_ls,list2=p_val_ls,xlim=None,ylim=None,x_y_lim_same=0,title_str="",title_suptitle='right',use_scale_x_y_lim=1,color=None,xlabel=None,ylabel=None)
    """


    pool = Pool(90)
    results = pool.map(multi_run_wrapper_for_corr,[(i,j,k) for i,j,k in zip([df]*shuffle_n,[l[0] for l in gene_corr_pair_ls],[l[1] for l in gene_corr_pair_ls])])
    pool.close()
    
    random_corr_txt = tpm_csv.replace('csv','random_corr.txt')
    TXT = open(random_corr_txt,'w')
    header_ls = ['id1','id2','Pearson','p_val']
    print >>TXT,'\t'.join(header_ls)
    for i,j in zip(gene_corr_pair_ls,results):
        print >>TXT,'\t'.join(i+map(str,j))

    TXT.close()
    printFuncRun('proteinatlas_random_corr')

def multi_run_wrapper_for_corr(args):
        return df_gene_id_corr(*args)

def df_gene_id_corr(df,id1,id2):
    df_gene1 = df[df['Gene']==id1]
    df_gene2 = df[df['Gene']==id2]
    df_gene = pd.merge(df_gene1,df_gene2,how='outer',on='Sample').replace(0,0.001)
    #rho, pval = stats.spearmanr(df_gene['Value_x'],df_gene['Value_y'])
    rho, pval = stats.pearsonr(df_gene['Value_x'],df_gene['Value_y'])
    return rho, pval

def proteinatlas_wgcna(mode='celline',proteinatlas_dir='/Share/home/zhangqf5/gongjing/DNA-RNA-Protein-interaction-correlation-12-18/data/proteinatlas',id_ls=None,process=2):
    printFuncRun('proteinatlas_wgcna')
    printFuncArgs()
    tpm_csv = proteinatlas_dir+'/'+'rna_'+mode+'.csv'
    tpm_pivot = tpm_csv.replace('csv','pivot.txt')

    if process == 1:
        df = pd.read_csv(tpm_csv,sep=',',header=0)
        print df.head()
        df_pivot = df.pivot(index='Gene', columns='Sample', values='Value')
        print df_pivot.head(20) 
        df_pivot.to_csv(tpm_pivot,sep='\t',header=True,index=True)
       
        return df_pivot

    if process == 2:
        if id_ls is None:
            id_ls = ['ENSG00000141510','ENSG00000146648']
            id_ls = ['ENSG00000114805','ENSG00000163629','ENSG00000143858','ENSG00000153707','ENSG00000176771','ENSG00000153266']
        df_pivot = pd.read_csv(tpm_pivot,header=0,sep='\t')
        print df_pivot.head()
        sample_mean = df_pivot.mean() # consider exp_val = 0
        sample_mean_df = sample_mean.to_frame(name='sample_mean')
        df_pivot_genes = df_pivot[df_pivot['Gene'].isin(id_ls)].T
        df_pivot_genes.columns = df_pivot_genes.loc['Gene',:]
        
        print sample_mean_df
        print df_pivot_genes
        df_merge = pd.merge(sample_mean_df,df_pivot_genes,right_index=True,left_index=True)
        print df_merge

        df_merge.replace(0,0.001,inplace=True)
        df_merge = df_merge.applymap(lambda x:np.log2(x))
        col_dict = {}
        id_vs_bg_dict = {} 
        for i in df_merge.columns[1:]:
            above_mean_num = df_merge[df_merge[i] >= df_merge['sample_mean']].shape[0]
            col_dict[i] = "%s (>=mean: %s)"%(i,above_mean_num)
            id_vs_bg_dict[i] = above_mean_num
        df_merge.rename(columns=col_dict,inplace=True)
 
        fig,ax=plt.subplots(figsize=(54,16))
        df_merge.plot(ax=ax,marker='o')
        #df_merge.plot(ax=ax)
        savefn = proteinatlas_dir+'/'+'exp_corr'+'/'+'background_'+'_'.join(id_ls)+'.png'
        ax_legend_set_text_line_color(ax)
        ax.set_xticklabels([])
        ax.set_xlabel("Samples(%s)"%(df_merge.shape[0]))
        ax.set_ylabel("Log2(TPM)")
        ax.set_title("Expression vs background in different %s"%(mode))
        plt.tight_layout()
        plt.savefig(savefn)
        plt.close()        
    
        return id_vs_bg_dict

    printFuncRun('proteinatlas_wgcna')

def goatools_run(ensembl_ls,organism='human',savefn=None,bg_type='human_all'):
    printFuncRun('goatools_run')
    printFuncArgs()
    gene2ensembl_dict = goatools_genelist_enrich.read_ncbi_gene2ensembl(organism=organism)
    print_dict_item0(gene2ensembl_dict)
    ncbi_ls = [gene2ensembl_dict[i]['GeneID'] for i in ensembl_ls if gene2ensembl_dict.has_key(i)]
    goea_results_sig,goea_results_sig_tsv_fn = goatools_genelist_enrich.run_GOEA(gene_ls=ncbi_ls,savefn=savefn,bg_type=bg_type)

    goatools_results_tsv_barplot(goea_results_sig_tsv_fn)

    printFuncRun('goatools_run')

def goatools_results_tsv_barplot(tsv):
    printFuncRun('goatools_results_tsv_barplot')
    printFuncArgs()
    if not os.path.isfile(tsv):
        print "[not a file] %s"%(tsv)
        return 
    df_tsv = pd.read_csv(tsv,sep='\t',header=0)
    df_tsv['-log10(p_fdr_bh)'] = -np.log10(df_tsv['p_fdr_bh'])

    df_tsv_head10 = pd.concat([df_tsv[df_tsv['NS']==i].iloc[0:10,:] for i in ['BP','CC','MF']],axis=0)
    print df_tsv_head10

    fig,ax = plt.subplots(figsize=(16,16))
    g = sns.barplot(x='-log10(p_fdr_bh)',y='name',hue='NS',data=df_tsv_head10,ax=ax)
    ax.set_xlabel('-log10(pval)')
    ax.set_ylabel('GO term name')
    ax.set_title('GO annottaion')
    fig.tight_layout()
    fig.savefig(tsv.replace('tsv','png'))
    #plt.tight_layout()
    #plt.savefig(goea_results_sig_tsv_fn+'.png')
    plt.close()
    printFuncRun('goatools_results_tsv_barplot')

##################################################
### math/calculation related
def choose(n, k):
    """
    A fast way to calculate binomial coefficients by Andrew Dalke (contrib).
    """
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in xrange(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0



##################################################
### statistic test related
# test based on poisson distribution
def sig_test_poisson(x=0.3,mu=0.001):
    p=stats.distributions.poisson.pmf(x,mu)
    #print "poisson test:",p # work correctly
    return p

# test based on negative binomial distribution
def sig_test_nbinom(x=100,n=50,p=0.3):
    p=stats.distributions.nbinom.pmf(x,n,p)
    #print "nbinom test:",p # work correctly
    return p

# two sample distribution test
def ks_2samp(x,y):
    p=stats.ks_2samp(x,y)[1]
    return p

# two sample rank test
def sig_spearman_corr(x,y):
    p=stats.spearmanr(x,y)[0]
    return p

def sig_pearson_corr(x,y):
    p=stats.pearsonr(x,y)[0]
    return p

##################################################
### multiprocessing
# 
def multi_run_wrapper(args): 
    return add(*args) 

def add(x,y):
    return x+y

# calling method
# pool=Pool(4)
# results = pool.map(multi_run_wrapper,[(1,2),(2,3),(3,4)])  # call on multi_run_wrapper rather than add
# print results

################################################## 
### sklearn related
# kmeans clustering
def kmeans_multi_run_wrapper(args): 
    return sklearn_kmeans_one(*args) 

def sklearn_kmeans(data=None,data_fn=None,cols_keep=None,cols_annotate=None,savefn=None,n_cluster_ls=[2,3,4,5,6,7,8,9],init_method_ls=['k-means++','random'],n_init_tries_ls=[20]):
    printFuncRun('sklearn_kmeans')
    printFuncArgs()
    if data is None and data_fn is None:
        print "data & data_fn are both None, please provide either of them"
    if data_fn is None:
        print "data_fn is None, will use data directly"
    if data is None:
        print "data is None, will read from data_fn"
        data = pd.read_csv(data_fn,sep='\t',header=0,index_col=0)
        if savefn is None:
            savefn = data_fn+'.kmeans.performace.png'
    if cols_annotate is not None:
       data_copy = data.copy()
    if not cols_keep is None:
       data = data[cols_keep]
    cols_non_annotate = [i for i in data.columns if i not in cols_annotate]
    data = data[cols_non_annotate]
    print data.info()
 
    n_samples,n_features = data.shape
    print "row samles: %s, features: %s"%(n_samples,n_features)

    result_ls = []
    header_ls = ['method_name','time_elapse','n_cluster','n_init_tries','iner','silhouette_score']
    print('\t'.join(header_ls))

    para_set_ls = []
    para_set_full_ls = []
    for n_cluster in n_cluster_ls:
        for n_init_tries in n_init_tries_ls:
            for init_method in init_method_ls:
                para_set_ls.append((n_cluster,n_init_tries,init_method))
                para_set_full_ls.append((data,data_copy,data_fn,n_cluster,n_init_tries,init_method))

    if len(para_set_ls) < 15:
        print "total: %s,will process one by one"%(len(para_set_ls))
        for n_cluster,n_init_tries,init_method in para_set_ls:
            name,t1,n_cluster,n_init_tries,estimator_inertia,silhouette_score,df = sklearn_kmeans_one(data,data_copy,data_fn,n_cluster,init_method,n_init_tries)
            result_ls.append([name,t1,n_cluster,n_init_tries,estimator_inertia,silhouette_score])
    else:
        print "total: %s,will process in paralle"%(len(para_set_ls))
        print para_set_full_ls
        pool=Pool(4)
        result_ls = pool.map(kmeans_multi_run_wrapper,para_set_full_ls)

    result_df = pd.DataFrame(result_ls)
    result_df.columns = header_ls
    print(result_df)

    g = sns.FacetGrid(data=result_df,col='method_name',
                      row='n_init_tries',hue='method_name',sharex=False,sharey=True,margin_titles=True, size=5)
    g.map(sns.barplot,'n_cluster','silhouette_score')
    plt.savefig(savefn)
    plt.close()
    printFuncRun('sklearn_kmeans')

    return data_copy

"""
def sklearn_kmeans_one(data,data_copy,data_fn,n_cluster,init_method,n_init_tries):
                if init_method == 'PCA-based':
                    t0=time.time()
                    n_features = data.shape[1]
                    n_cluster = min(n_cluster,n_features)
                    pca = PCA(n_components=n_cluster).fit(data)
                    estimator = KMeans(init=pca.components_, n_clusters=n_cluster, n_init=n_init_tries)
                    estimator.fit(data)
                    cluster_label = estimator.labels_
                    df = data.copy()
                    df.index = cluster_label
                    df.sort_index(inplace=True)
                    sns_heatmap(df,df_fn=None,cols_keep=None,cols_drop=None,savefn=data_fn+'_'+str(n_cluster)+'_'+init_method+'_'+str(n_init_tries)+'.png')
                    
                    data_copy.index = cluster_label
                    data_copy.sort_index(inplace=True)
                    data_copy['cluster'] = data_copy.index
                    sns_heatmap_annotate(df=data_copy,cols_annotate=['TE_type','cluster'],savefn=data_fn+'_'+str(n_cluster)+'_'+init_method+'_'+str(n_init_tries)+'.anno.png')
                    silhouette_score = metrics.silhouette_score(data,estimator.labels_,metric='euclidean')
                    t1=time.time()-t0
                    name = "PCA-based"
                    print(' %9s %.2fs %s %s %i %.3f'%(name,t1,n_cluster,n_init_tries,estimator.inertia_,silhouette_score))
                    #result_ls.append([name,t1,n_cluster,n_init_tries,estimator.inertia_,silhouette_score])
                    return name,t1,n_cluster,n_init_tries,estimator.inertia_,silhouette_score

                t0=time.time()
                estimator = KMeans(n_clusters=n_cluster,init=init_method,n_init=n_init_tries,random_state=0)
                estimator.fit(data)
                cluster_label = estimator.labels_
                df = data.copy()
                df['cluster'] = cluster_label
                df['id'] = df.index
                df.index = df['cluster']
                df.sort_index(inplace=True)
                df.to_csv(data_fn+'_'+str(n_cluster)+'_'+init_method+'_'+str(n_init_tries)+'.label',sep='\t',header=True,index=False)
               
                del df['cluster']
                del df['id']
                sns_heatmap(df,df_fn=None,cols_keep=None,cols_drop=None,savefn=data_fn+'_'+str(n_cluster)+'_'+init_method+'_'+str(n_init_tries)+'.png')
                print df.head(2)

                data_copy.index = cluster_label
                data_copy.sort_index(inplace=True)
                data_copy['cluster'] = data_copy.index
                sns_heatmap_annotate(df=data_copy,cols_annotate=['TE_type','cluster'],savefn=data_fn+'_'+str(n_cluster)+'_'+init_method+'_'+str(n_init_tries)+'.anno.png')

                silhouette_score = metrics.silhouette_score(data,estimator.labels_,metric='euclidean')
                t1=time.time()-t0
                name = init_method
                print(' %9s %.2fs %s %s %i %.3f'%(name,t1,n_cluster,n_init_tries,estimator.inertia_,silhouette_score))
                #result_ls.append([name,t1,n_cluster,n_init_tries,estimator.inertia_,silhouette_score])
                return name,t1,n_cluster,n_init_tries,estimator.inertia_,silhouette_score,df
"""

def sklearn_kmeans_one(data,data_copy,data_fn,n_cluster,init_method,n_init_tries):
                if init_method == 'PCA-based':
                    t0=time.time()
                    n_features = data.shape[1]
                    n_cluster = min(n_cluster,n_features)
                    pca = PCA(n_components=n_cluster).fit(data)
                    estimator = KMeans(init=pca.components_, n_clusters=n_cluster, n_init=n_init_tries)
                    estimator.fit(data)
                    cluster_label = estimator.labels_
                    df = data.copy()
                    df.index = cluster_label
                    df.sort_index(inplace=True)
                    sns_heatmap(df,df_fn=None,cols_keep=None,cols_drop=None,savefn=data_fn+'_'+str(n_cluster)+'_'+init_method+'_'+str(n_init_tries)+'.pdf')
                    
                    data_copy.index = cluster_label
                    data_copy.sort_index(inplace=True)
                    data_copy['cluster'] = data_copy.index
                    sns_heatmap_annotate(df=data_copy,cols_annotate=['TE_type','cluster'],savefn=data_fn+'_'+str(n_cluster)+'_'+init_method+'_'+str(n_init_tries)+'.anno.png')
                    silhouette_score = metrics.silhouette_score(data,estimator.labels_,metric='euclidean')
                    t1=time.time()-t0
                    name = "PCA-based"
                    print(' %9s %.2fs %s %s %i %.3f'%(name,t1,n_cluster,n_init_tries,estimator.inertia_,silhouette_score))
                    #result_ls.append([name,t1,n_cluster,n_init_tries,estimator.inertia_,silhouette_score])
                    return name,t1,n_cluster,n_init_tries,estimator.inertia_,silhouette_score

                t0=time.time()
                estimator = KMeans(n_clusters=n_cluster,init=init_method,n_init=n_init_tries,random_state=0)
                estimator.fit(data)
                cluster_label = estimator.labels_
                df = data.copy()
                df['cluster'] = cluster_label
                df['id'] = df.index
                df.index = df['cluster']
                df.sort_index(inplace=True)
                df.to_csv(data_fn+'_'+str(n_cluster)+'_'+init_method+'_'+str(n_init_tries)+'.label',sep='\t',header=True,index=False)
               
                del df['cluster']
                del df['id']
                #sns_heatmap(df,df_fn=None,cols_keep=None,cols_drop=None,savefn=data_fn+'_'+str(n_cluster)+'_'+init_method+'_'+str(n_init_tries)+'.pdf',figsize_x=24,figsize_y=48)
                print df.head(2)

                data_copy.index = cluster_label
                data_copy.sort_index(inplace=True)
                data_copy['cluster'] = data_copy.index
                sns_heatmap_annotate(df=data_copy,cols_annotate=['cluster'],savefn=data_fn+'_'+str(n_cluster)+'_'+init_method+'_'+str(n_init_tries)+'.anno.pdf')

                silhouette_score = metrics.silhouette_score(data,estimator.labels_,metric='euclidean')
                t1=time.time()-t0
                name = init_method
                print(' %9s %.2fs %s %s %i %.3f'%(name,t1,n_cluster,n_init_tries,estimator.inertia_,silhouette_score))
                #result_ls.append([name,t1,n_cluster,n_init_tries,estimator.inertia_,silhouette_score])
                return name,t1,n_cluster,n_init_tries,estimator.inertia_,silhouette_score,df

##################################################
# print test
def printt():
    print "test from gj.py"



##################################################
### main function
def main():
    #check_sam_cigar_pattern(input_f="/Share/home/zhangqf/gongjing/paris-2016-05/result/SRR2814761/Aligned.out.sam")
    #check_sam_cigar_pattern(input_f="/Share/home/zhangqf/gongjing/paris-2016-05/result/mES_3rep_union/mES_3repUnion_DG_rmNGDGempty.sam")
    #check_sam_cigar_pattern(input_f="/Share/home/zhangqf/gongjing/paris-2016-05/result/mES_3rep_union/SRR2814766_DG_rmNGDGempty.sam")
    #check_sam_cigar_pattern(input_f="/Share/home/zhangqf/gongjing/paris-2016-05/result/SRR2814711/Aligned.out.sam")
    #check_sam_cigar_pattern(input_f="/Share/home/zhangqf/gongjing/paris-2016-05/result/SRR2814710/Aligned.out.sam")
    #plot_trun_site_mutation()
    #dict_bar_plot_by_dataframe()
    #getFastaFromBed()
    #getFastaFromBed(outfmt='txt')
    #getFastaFromBed(strand=1)
    #getFastaFromBed(strand=1,extend_left=30,extend_right=50)
    #sig_test_poisson()
    #sig_test_nbinom()
    #sam_keep_alignment_first_n_hits()
    #sam2sortedbam()
    #multiBamCorrelationDeeptool(bam_file_ls=["result/SRR2814766/Aligned.out.nongapped.sorted.bam","result/SRR2814767/Aligned.out.nongapped.sorted.bam","result/SRR2814768/Aligned.out.nongapped.sorted.bam","result/SRR2814766/Aligned.out.gapped.sorted.bam","result/SRR2814767/Aligned.out.gapped.sorted.bam","result/SRR2814768/Aligned.out.gapped.sorted.bam"],binsize=1000,union_out_npz="result/mES_3rep_union/paris_gapped_vs_nongapped_66_67_68.npz")
    #bedOverlapBySetOperation()
    #read_snp_into_dict('mouse','ucsc')
    #dict_bar_plot_by_dataframe_mutlti()
    #bedOverlapBySetOperation(input_bed1="result/SRR2814766/Aligned.out.gapped.sorted.bam.igvtoolsCallStrand.all.bed",input_bed2="result/SRR2814766/Aligned.out.nongapped.sorted.bam.igvtoolsCallStrand.all.bed",intersect_bed_info_txt="result/SRR2814766/gapped_vs_ungapped_overlap.txt")
    #bedOverlapBySetOperation(input_bed1="result/SRR2814767/Aligned.out.gapped.sorted.bam.igvtoolsCallStrand.all.bed",input_bed2="result/SRR2814767/Aligned.out.nongapped.sorted.bam.igvtoolsCallStrand.all.bed",intersect_bed_info_txt="result/SRR2814767/gapped_vs_ungapped_overlap.txt")
    #bedOverlapBySetOperation(input_bed1="result/SRR2814768/Aligned.out.gapped.sorted.bam.igvtoolsCallStrand.all.bed",input_bed2="result/SRR2814768/Aligned.out.nongapped.sorted.bam.igvtoolsCallStrand.all.bed",intersect_bed_info_txt="result/SRR2814768/gapped_vs_ungapped_overlap.txt")
    #bedtoolsBedBamOverlap()
    #bedtoolsBedBamOverlap(input_bam="result/SRR2814766/Aligned.out.gapped.sorted.bam",input_bed="result/SRR2814766/gap_1_nongap_0.bed",out_bam='')
    #bedtoolsBedBamOverlap(input_bam="result/SRR2814766/Aligned.out.gapped.sorted.bam",input_bed="result/SRR2814766/gap_0_nongap_1.bed",out_bam='')
    #bamStatsBySamtoolsPlot()
    #bedAnnotationByGTF()
    #bedAnnotationBySnpBed()
    #starMapLogPlot(sample_ls=['SRR2814761','SRR2814762','SRR2814763','SRR2814764','SRR2814765','SRR2814766','SRR2814767','SRR2814768','hela_S1_10','hela_S1_20','hela_S1_2_1','hela_S1_2_2'])
    #samMapPlot(sample_ls=['SRR2814761','SRR2814762','SRR2814763','SRR2814764','SRR2814765','SRR2814766','SRR2814767','SRR2814768','hela_S1_10','hela_S1_20','hela_S1_2_1','hela_S1_2_2'])
    #samMapPlot(sample_ls=['SRR2814700'])
    #plotSamtoolsStatsCustome(sample_ls=['SRR2814766','SRR2814767','SRR2814768'])
    #bedOverlapBySetOperation(input_bed1="result/SRR2814700/Aligned.out.sorted.bam.igvtoolsCallStrand.all.bed",input_bed2="data/supplementary_data/icshape_data/mES/icSHAPE.sim.bedgraph",intersect_bed_info_txt="",plotvenn=1,plotdist=1)
    #bedOverlapBySetOperation(input_bed1="result/SRR2814766/Aligned.out.sorted.bam.igvtoolsCallStrand.all.bed",input_bed2="data/supplementary_data/icshape_data/mES/icSHAPE.sim.bedgraph",intersect_bed_info_txt="",plotvenn=1,plotdist=1)
    #plotSamtoolsStats2(sample_ls=['SRR2814766','SRR2814767','SRR2814768'],bam_stats_ls=['Aligned.out.gapped.sorted.bam.stats','Aligned.out.nongapped.sorted.bam.stats'])
    #plotSamtoolsStats2(sample_ls=['SRR2814763','SRR2814764','SRR2814765'],bam_stats_ls=['Aligned.out.gapped.sorted.bam.stats','Aligned.out.nongapped.sorted.bam.stats'])
    #plotSamtoolsStats2(sample_ls=['SRR2814761','SRR2814762'],bam_stats_ls=['Aligned.out.gapped.sorted.bam.stats','Aligned.out.nongapped.sorted.bam.stats'])
    #vcf_pattern_comparison(sample_ls=['SRR2814711','SRR2814799','SRR2814766','SRR2814767','SRR2814768'],vcf_ls=['GATK/filter.output.vcf','SamtoolsBcftoolsCallMutation/Aligned.out.bam.vcf','SamtoolsBcftoolsCallMutation/Aligned.out.bam.VarScan.vcf'])
    #df_to_table_png(df=pd.DataFrame({'A':[1,2,3],'B':[4,5,6]}),fig=None,ax=None,colLabels=None,rowLabels=None,cellLoc='center',colwidth=0.17,savefn='df_to_table_png_test.png')
    #list2fa(seq_ls=['atcgt','ggtaa','gtcag','ggtag'],savefn='scripts/test_list2fa.fa')
    #plotFaWeblogo(input_fa="scripts/test_list2fa.fa",weblogo_f_title="")
    #ax_show_img(fig=None,ax=None,img_fn='scripts/test_list2fa.weblogo.png',savefn='test_from_img.png')
    #list_list_box_plot(list_list=[[1,2,3,4,5],[11,12,3,4,10,5,3,1,2],[20,11,33,22],range(3,30),np.arange(10,20,0.01)],list_label=['A','B','C','D','E','F'],fig=None,ax=None,savefn='test_list_bar.png',title="test",xlabel=None,ylabel=None)
    printt()
    sklearn_kmeans(data=None,data_fn='/Share/home/zhangqf/gongjing/paris-2016-05/result/sunlei/data/Riboseq_GSE73136_HEK293_2015_NatureMethod/transcript_feature_combine.gini.combine100.txt',cols_keep=['UTR5','CDS_H','CDS_M','CDS_T','UTR3'],savefn=None)
    #sns_heatmap(df=None,df_fn='/Share/home/zhangqf/gongjing/paris-2016-05/result/sunlei/data/Riboseq_GSE73136_HEK293_2015_NatureMethod/transcript_feature_combine.gini.combine100.txt',cols_keep=['UTR5','CDS_H','CDS_M','CDS_T','UTR3'],savefn='/Share/home/zhangqf/gongjing/paris-2016-05/result/sunlei/data/Riboseq_GSE73136_HEK293_2015_NatureMethod/transcript_feature_combine.gini.combine100.txt.heatmap.pdf')

##################################################
### final recall
if __name__ == "__main__":
    main()
