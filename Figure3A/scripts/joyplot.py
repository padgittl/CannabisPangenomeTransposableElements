import sys, re, os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from joypy import joyplot
import pandas as pd


########
# MAIN #
########

def createJoyplot(csvFile,teKeyword):
    values = []
    with open(csvFile,'r') as F:
        for line in F:
            line = line.strip().split(',')
            if 'Names' not in line:
                for i in line[:-1]:
                    values.append(float(i))
    maxValue = max(values)
    minValue = min(values)
    
    colorDict = {'CACTA':'#000000','Mutator':'#E69F00','PIF_Harbinger':'#0072B2','Tc1_Mariner':'#009E73','hAT':'#F0E442','helitron':'#CC79A7', 'Ty1-LTR':'#542788',
                 'unknown-LTR':'#998ec3', 'Ty3-LTR':'#b35806', 'Total':'gray'}
    colorList =	[]
    df = pd.read_csv(csvFile, nrows=1)
    for i in list(df.columns):
        if i != 'Names':
            if i in colorDict:
                colorID = colorDict[i]
                colorList.append(colorID)
            else:
                print(i,"not in colorDict")
                sys.exit()
    if len(colorList) == 1:
        colorList = colorList[0]
    if 'dna' in teKeyword:
        minValue = -0.25
        linewidth=0.0
    else:
        linewidth=0.25
    
    data = pd.read_csv(csvFile)
    #plt.rcParams.update({'font.size': 4})
    plt.rcParams.update({'font.size': 6})

    # overlap=1
    # figsize=(3.5,5)
    # dna and ltr
    # fig, axes = joyplot(data, by="Names", linewidth=linewidth, linecolor='white', legend=True, overlap=1.5, color=colorList, figsize=(5,4), x_range=[minValue,maxValue], title="Percent of each scaffolded assembly\ncovered by transposable elements", ylim='own', loc="lower center")
    fig, axes = joyplot(data, by="Names", linewidth=linewidth, linecolor='white', legend=True, overlap=1.5, color=colorList, figsize=(3.5,5), x_range=[minValue,maxValue], title="Percent of each scaffolded assembly\ncovered by transposable elements", ylim='own', loc="best")
    # axes.set_xlabel('Percent of each scaffolded assembly covered by transposable elements')
    # title="Percent of each scaffolded assembly covered by transposable elements")
    #plt.tight_layout()
    plt.savefig(teKeyword + '_joyplot.png',dpi=600)
    plt.savefig(teKeyword + '_joyplot.svg')
    plt.close()
    

usage = "Usage: " + sys.argv[0] + " <csv file> <TE keyword>\n"
if len(sys.argv) != 3:
    print(usage)
    sys.exit()

csvFile = sys.argv[1]
teKeyword = sys.argv[2]

createJoyplot(csvFile,teKeyword)
