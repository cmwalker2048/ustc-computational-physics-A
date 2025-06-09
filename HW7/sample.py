import numpy as np # 数据处理
import matplotlib.pyplot as plt # 作图
import csv # 读取数据的工具

# 实验数据归一化函数
def p(datapath,outpath):
    x0 = np.loadtxt(datapath,dtype=np.float32,delimiter=",",usecols=(0),encoding='utf-8')
    y0 = np.loadtxt(datapath,dtype=np.float32,delimiter=",",usecols=(1),encoding='utf-8')
    sum = 0 #求和寄存器
    for i in range(0,len(y0)):
        sum+=y0[i]
    y = y0/sum # 更新 y = y0/sum(y0)
    # 写入文件循环
    with open(outpath,'w') as out: # 参数为:写入文件路径,写入模式
        for j in range(0,len(y0)):
            out.write("%.4f,%.4f\n"%(x0[j],y[j])) #写入文件

# 直接抽样函数
    # filepath随机数文件地址;datapath归一化p(x)地址;outpath输出文件地址;
def direct(filepath,datapath,outpath,n):# n为抽样的点数
    # 将CSV文件中的随机数导入数组orign中
    orign = np.loadtxt(filepath,dtype=np.float32,delimiter=",",usecols=(0),encoding='utf-8')
    # 归一化文件地址
    x0 = np.loadtxt(datapath,dtype=np.float32,delimiter=",",usecols=(0),encoding='utf-8')
    y0 = np.loadtxt(datapath,dtype=np.float32,delimiter=",",usecols=(1),encoding='utf-8')
    # 用于储存结果
    x = []
    # 抽样过程
    for j in range(0,n):
        sum =0
        for k in range(0,len(y0)):
            sum += y0[k]
            if (orign[j]<sum):
                x.append(x0[k])
                break
    
    # 写入文件循环
    with open(outpath,'w') as out: # 参数为:写入文件路径,写入模式
        for j in range(0,n):
            out.write("%.4f\n"%(x[j])) #写入文件
        
    
    
# 简单分布主函数 (最终未使用，请忽略！)
    # filepath随机数文件地址;datapath为数据地址;outpath输出文件地址;
def simple(filepath,datapath,outpath,n):# n为抽样个数
    # 将CSV文件中的随机数导入数组orign中
    orign = np.loadtxt(filepath,dtype=np.float32,delimiter=",",usecols=(0),encoding='utf-8')
    # 原始实验数据
    x0 = np.loadtxt(datapath,dtype=np.float32,delimiter=",",usecols=(0),encoding='utf-8')
    y0 = np.loadtxt(datapath,dtype=np.float32,delimiter=",",usecols=(1),encoding='utf-8')
    # 简单分布参数
    a = 2900
    b = 3013
    M = 38000
    xi_1 = orign*(b-a)+a
    eta = orign
    # 用于计点数的寄存器
    counter = 0
    # 用于舍去点的计数器
    reject = 0
    # 用于遍历随机数的数
    j = 0
    
    x = [] # 符合要求的x坐标
    
    while(counter <= n):
        # 求出xi_1在X轴上的区间范围
        i = 0
        while(x0[i] < xi_1[j]):
            i+=1
        u = eta[j+1]*M
        v = max(y0[i-1],y0[i])
        if (u<v):
            x.append(xi_1[j])
            counter+=1
        else:
            reject+=1
        j+=2
    
    # 写入文件循环
    with open(outpath,'w') as out: # 参数为:写入文件路径,写入模式
        for j in range(0,n-1):
            out.write("%.4f\n"%(x[j])) #写入文件
    
    print(counter/(counter+reject))
    '''
    plt.plot(x0,y0,'-')
    plt.hist(x,bins = 30,edgecolor='white',density=True)
    plt.xlabel('Energy($eV$)')
    plt.ylabel('Intensity(counts)')
    plt.title('Data')
    plt.show()
    '''


# 第二类舍选法主函数

    # filepath随机数文件地址;datapath为数据地址;outpath输出文件地址;
def trans(filepath,datapath,outpath,n):# n为抽样的点数
    # 将CSV文件中的随机数导入数组orign中
    orign = np.loadtxt(filepath,dtype=np.float32,delimiter=",",usecols=(0),encoding='utf-8')
    # 实验数据
    x0 = np.loadtxt(datapath,dtype=np.float32,delimiter=",",usecols=(0),encoding='utf-8')
    y0 = np.loadtxt(datapath,dtype=np.float32,delimiter=",",usecols=(1),encoding='utf-8')
    # F(x)的一些参数设置
    x_1 = 0.554028
    x_2 = 0.946955
    k_1 = 169.667
    b_1 = 2900
    k_2 = 25.45
    b_2 = 2979.9
    k_3 = 169.667
    b_3 = 2843.333
    M_1 = 5700
    M_2 = 38000
    # 分部的随机变量
    xi_1_1 = orign*k_1+b_1
    xi_1_2 = orign*k_2+b_2
    xi_1_3 = orign*k_3+b_3
    xi = orign
    eta = orign
        
    # 用于计点数的寄存器
    counter = 0
    # 用于舍去的点的计数
    reject = 0
    # 用于遍历随机数的数
    j = 0
    
    x = [] # 符合要求的x坐标
    
    # 舍选法主循环
    while(counter <= n):
        if(xi[j]<x_1):
            out = M_1
            sample = xi_1_1[j]
        elif(x_1<xi[j]<x_2):
            sample = xi_1_2[j]
            out = M_2
        else:
            sample = xi_1_3[j]
            out = M_1
        # 求出xi_1在X轴上的区间范围
        i = 0
        while(x0[i] < sample):
            i+=1
        u = eta[j+1]*out
        v = max(y0[i-1],y0[i])
        if (u<v):
            x.append(sample) # 记录合格的样本
            counter+=1
        else:
            reject+=1
        j+=2
    
    # 写入文件循环
    with open(outpath,'w') as out: # 参数为:写入文件路径,写入模式
        for j in range(0,n-1):
            out.write("%.4f\n"%(x[j])) #写入文件
    
    # 输出效率
    print(counter/(counter+reject))
    
    '''
    plt.plot(x0,y0,'-')
    plt.hist(x,bins = 30,edgecolor='white',density=0)
    plt.xlabel('Energy($eV$)')
    plt.ylabel('Intensity(counts)')
    plt.title('Data')
    plt.show()
    '''
# 作图函数，请忽略！
def temp_plot(filepath1,filepath2):
    #orign = np.loadtxt(filepath1,dtype=np.float32,delimiter=",",usecols=(0),encoding='utf-8')
    x = np.loadtxt(filepath1,dtype=np.float32,delimiter=",",usecols=(0),encoding='utf-8') 
    y = np.loadtxt(filepath1,dtype=np.float32,delimiter=",",usecols=(1),encoding='utf-8')
    data = np.loadtxt(filepath2,dtype=np.float32,delimiter=",",usecols=(0),encoding='utf-8')
    plt.plot(x,y,'-',label='original data')
    plt.xlabel('$Energy(eV)$')
    plt.ylabel('F(x)')
    plt.title("direct $10^2$")
    plt.hist(data[0:100:1],bins = 100,edgecolor='white',density=1,label='direct')
    #sns.kdeplot(orign)
    plt.legend(loc=2)
    plt.show()

# 利用直接抽样法
#direct('./16807.csv','./unitlize.csv','./direct.csv',10000000)
# 利用简单分布
#simple('./16807.csv','./data.csv','./simple.csv',1000000)
# 利用第二类舍选法
#trans('./16807.csv','./data.csv','./trans.csv',100000)
# 暂时作图
temp_plot('./unitlize.csv','./direct.csv')