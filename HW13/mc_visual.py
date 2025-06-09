import numpy as np
import matplotlib.pyplot as plt

# 注：此文件为作图文件，助教可以不用检查！

# 一些图片的标题
title_1 = "$\Delta I /I$-$\log_{10}{\gamma}$ relationship $\left( \gamma \in \left[1,1.25\\right] \\right)$"
title_2 = "$\Delta I /I$-$\log_{10}{N}$ relationship"
title_3 = "$\Delta I /I$-$\log_{10}{\gamma}$ relationship $\left( \gamma \in \left[1.6,2.5\\right], x_0 = 1000 \\right)$"
title_4 = "$\Delta I /I$-$\log_{10}{\gamma}$ relationship $\left( \gamma \in \left[1.6,2.8\\right], k = 0.4 \\right)$"

def gamma(filepath):
    gamma = np.loadtxt(filepath,dtype=np.float32,delimiter=",",usecols=(0),encoding='utf-8')
    I = np.loadtxt(filepath,dtype=np.float32,delimiter=",",usecols=(1),encoding='utf-8')
    delta = np.loadtxt(filepath,dtype=np.float32,delimiter=",",usecols=(2),encoding='utf-8')
    
    x = np.log10(gamma)
    #plt.figure(figsize=(20, 20))
    plt.title(title_1)
    plt.xlabel("$\log_{10}{\gamma}$")
    plt.ylabel("$\Delta I / I$")
    plt.plot(x,delta,label='Curve')
    plt.scatter(x,delta,c='red',s=6,label='Data')
    plt.legend(loc=3,framealpha=1, shadow=True)
    #plt.savefig(file_pic)
    plt.show()
    
def efficient(filepath):
    gamma = np.loadtxt(filepath,dtype=np.float32,delimiter=",",usecols=(0),encoding='utf-8')
    I = np.loadtxt(filepath,dtype=np.float32,delimiter=",",usecols=(1),encoding='utf-8')
    delta = np.loadtxt(filepath,dtype=np.float32,delimiter=",",usecols=(2),encoding='utf-8')
    
    x = np.log10(gamma)
    #plt.figure(figsize=(20, 20))
    plt.title(title_2)
    plt.xlabel("$\log_{10}{N}$")
    plt.ylabel("$\Delta I / I$")
    plt.plot(x,delta,label='Curve')
    plt.scatter(x,delta,c='red',s=6,label='Data')
    plt.legend(loc=3,framealpha=1, shadow=True)
    #plt.savefig(file_pic)
    plt.show()
    pass
    
def x(filepath):
    gamma = np.loadtxt(filepath,dtype=np.float32,delimiter=",",usecols=(0),encoding='utf-8')
    I = np.loadtxt(filepath,dtype=np.float32,delimiter=",",usecols=(1),encoding='utf-8')
    delta = np.loadtxt(filepath,dtype=np.float32,delimiter=",",usecols=(2),encoding='utf-8')
    
    x = np.log10(gamma)
    #plt.figure(figsize=(20, 20))
    plt.title(title_3)
    plt.xlabel("$\log_{10}{\gamma}$")
    plt.ylabel("$\Delta I / I$")
    plt.plot(x,delta,label='Curve')
    plt.scatter(x,delta,c='red',s=6,label='Data')
    plt.legend(loc=1,framealpha=1, shadow=True)
    #plt.savefig(file_pic)
    plt.show()
    
def k(filepath):
    gamma = np.loadtxt(filepath,dtype=np.float32,delimiter=",",usecols=(0),encoding='utf-8')
    I = np.loadtxt(filepath,dtype=np.float32,delimiter=",",usecols=(1),encoding='utf-8')
    delta = np.loadtxt(filepath,dtype=np.float32,delimiter=",",usecols=(2),encoding='utf-8')
    
    x = np.log10(gamma)
    #plt.figure(figsize=(20, 20))
    plt.title(title_4)
    plt.xlabel("$\log_{10}{\gamma}$")
    plt.ylabel("$\Delta I / I$")
    plt.plot(x,delta,label='Curve')
    plt.scatter(x,delta,c='red',s=6,label='Data')
    plt.legend(loc=1,framealpha=1, shadow=True)
    #plt.savefig(file_pic)
    plt.show()

#gamma("./k_0025.csv")
#efficient("./temp.csv")
#x("./x_0_1000_s.csv")
k("./k_04_s.csv")