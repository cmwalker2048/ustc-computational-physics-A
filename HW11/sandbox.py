import numpy as np  #数据处理
import matplotlib.pyplot as plt #作图
from scipy.linalg import lstsq  #最小二乘法拟合

# 文件地址
#filedata = "./dla_sand_3_3_50000.csv"
filedata = "./dbm_sand_3_5_200.csv"
# 作图标题
#title = "$DLA,K=3,ALPHA = 3,N=50000$"
title = "$DBM,K=3,\eta=5,L=200$"
# 储存图片路径
#picpath = "./dla_sand_3_3_50000.png"
picpath = "./dbm_sand_3_5_200.png"

# SandBox数据拟合函数
def sandbox(filepath):
    # 导入数据
    r = np.loadtxt(filepath,dtype=np.float32,delimiter=",",usecols=(0),encoding='utf-8')
    n = np.loadtxt(filepath,dtype=np.float32,delimiter=",",usecols=(1),encoding='utf-8')
    # 做ln函数的处理
    R = np.log(r)
    N = np.log(n)
    
    #'''
    # 以下部分为提取直线部分
    max = 0
    length = 0;
    thruster = 0.9;
    for i in range(len(n)):
        if(i>=1):
            k = (N[i]-N[0])/(R[i]-R[0])
            if(np.abs(k)>max):
                max = k
            # 若斜率改变1-thruster，则认为不再是直线
            if(np.abs(k/max)<thruster):
                length = i-1
                break
        length = i
        pass
    
    # 改变数组长度（其范围为直线部分）
    R_f = R[0:length:1]
    N_f = N[0:length:1]
    #'''
    
    '''
    mem = []
    down =0;
    up =0;
    flag = 0;
    thruster = 0.99;
    for i in range(len(n)):
        k = (N[i]-N[i-1])/(R[i]-R[i-1])
        if(i>=2 and (flag==0)):
            if(np.abs(k-mem[-1])<thruster):
                flag=1;
                down=i;
        if(i>=2 and flag):
            if(np.abs(k-mem[-1])>thruster):
                up=i;
                break;
        if(i!=0):
            mem.append(k);
        if(i==(len(n)-1)):
            up=i;
            
    R_f = R[down:up:1]
    N_f = N[down:up:1]
    '''
    
    # 引用Scipy库，进行最小二乘法拟合
    A = np.vstack([R_f**0,R_f**1])
    sol,_r,rank,s = lstsq(A.T,N_f)
    
    # 获得拟合数据
    m = 100
    x = np.linspace(int(R_f[0]),int(R_f[-1]),m)
    y = sol[0] + sol[1] * x
    
    # 作图部分，请忽略！
    string1 = '$K=%0.4f$'%(sol[1])
    plt.plot(R,N,'o',label="Orignal Data")
    plt.plot(x,y,label='Fitted Data')
    plt.title(title)
    plt.xlabel("$\log{R}$")
    plt.ylabel("$\log{N}$")
    plt.legend(loc=4)
    plt.text(R_f[0],N_f[-1],string1,fontsize=15,verticalalignment="top",horizontalalignment="left")
    print('%0.4f'%(sol[1]))
    plt.savefig(picpath)
    plt.show()

# 调用函数
sandbox(filedata)