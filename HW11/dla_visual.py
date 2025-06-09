import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
#从pyplot导入MultipleLocator类，这个类用于设置刻度间隔

file_data = "./dla_3_3_50000.csv"
file_pic = "./dla_3_3_50000.png"
title = "$DLA,K=3,ALPHA=3,N=50000$"

def dla_visual(filepath,t):
    x = np.loadtxt(filepath,dtype=np.float32,delimiter=",",usecols=(0),encoding='utf-8')
    y = np.loadtxt(filepath,dtype=np.float32,delimiter=",",usecols=(1),encoding='utf-8')
    plt.figure(figsize=(20, 20))
    plt.title(title)
    plt.xlabel("$X$")
    plt.ylabel("$Y$")
    plt.scatter(x,y,s=1)
    plt.savefig(file_pic)
    plt.show()
    
dla_visual(file_data,100000)