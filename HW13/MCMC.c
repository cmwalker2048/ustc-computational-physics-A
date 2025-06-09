#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#define ALPHA 2             //常数 alpha = 2
#define BETA 1              //常数 beta = 1
#define ANSWER 2            //积分的答案 I_0 = 2
#define K 0.1               //热化系数 m = N * k
#define X_0 1               //初始点 x_0 = 1
#define G_MAX 100           //gamma 可取值的最大值
#define G_MIN 0.0001        //gamma 可取值的最小值
#define G_LEN 1000          //gamma于区间 [G_MIN,G_MAX]中的等分份数
#define LENTH 1000000       //Markov Chain 的链长

// ############# 以下均为输出用的输出文件路径 ################
#define FILE_ACCURATE_1 "./accurate_1.csv"
#define FILE_ACCURATE_2 "./accurate_2.csv"
#define FILE_ACCURATE_3 "./accurate_3.csv"
#define FILE_ACCURATE_4 "./accurate_4.csv"
#define FILE_X_0_1000 "./x_0_1000_s.csv"
#define FILE_X_0_100 "./x_0_100_s.csv"
#define FILE_EFFICIENT "./efficient.csv"
#define FILE_TEMP "./temp.csv"
#define FILE_K_0025 "./k_0025_s.csv"
#define FILE_K_04 "./k_04_s.csv"

// 随机数生成器，可以直接调用rn()使用，利用局部静态变量实现
float rn(){
    // 16807产生器的参数
    int A = 16807; //参数a
    int B = 0; //参数b
    long M = 2147483647; //参数m
    static long i = 23942907; //起始值i_0的大小(种子值)
    // Schrage法的参数
    long Q = 127773; //参数q（除数）
    int R = 2836; //参数r（余数）
    float x;
    while(1){
        x = (i*1.0) / M; //得到随机数
        i = A*(i % Q) - R*(i / Q);
        if(i < 0){
            i = i + M; //Schrage法的一个步骤
        }
        break;
    }
    return x;
}

// Metropolis-Hasting法中判别式 r 的表达式
double func(double x1,double x2,double gamma){
    // x2 为新值 ; x1 为旧值 ;
    return pow(x2/x1,(double)(ALPHA-1)) * exp(-(x2-x1)/(double)BETA) * exp((x2-x1)/(double)gamma);
}

//double型绝对值函数
double abs_d(double a){
    if(a>0){
        return a;
    }
    else{
        return -a;
    }
}

// Metropolis-Hasting方法主程序，包括了计算精度的部分
void mc_accurate(){
    int i,N,counter;    //一些计数器和 N值
    long j;             //长整数型计数器
    double *I;          //用于储存积分值的指针数组
    double *accurate;   //用于储存相对误差 delta I = Delta I / I的指针数组
    double *gamma;      //用于储存gamma值的指针数组
    double *x;          //用于储存Markov Chain的指针数组
    double step;        //等分gamma时用的等差数列的公差
    double x_t;         //产生新的x值的试探值 x_t
    double r;           //Metropolis-Hasting 法中的判别式 r
    double sum;         //计算积分 I 时的加法寄存器 sum
    double m;           //热化步数 m

    I = (double *)malloc(G_LEN*sizeof(double));         //初始化数组 I
    accurate = (double *)malloc(G_LEN*sizeof(double));  //初始化数组 accurate
    gamma = (double *)malloc(G_LEN*sizeof(double));     //初始化数组 gamma
    // 给gamma附上值，为[G_MIN,G_MAX]间的等差数列
    step = (G_MAX-G_MIN)*1.0/G_LEN; //产生公差
    for(i=0;i<G_LEN;i++){
        gamma[i] = G_MIN+step*i;
    }
    
    // 对于不同的gamma值，分别进行Metropolis法进行计算
    for(i=0;i<G_LEN;i++){

        N = 0;      //计数器的初始化
        counter=0;  //计数器的初始化
        sum=0;      //计数器的初始化

        x = (double *)malloc(LENTH*sizeof(double));     //初始化数组 x
        x[0] = X_0; //初始点 x_0 = 1 或其他值

        //Metropolis-Hasitng法主程序
        for(j=0;j<LENTH;j++){
            x_t = -log(rn()); // [0,infty]生成随机数
            r = func(x[j],x_t,gamma[i]); //计算判别式 r
            // 下面进行一些逻辑运算
            if(r>1){
                x[j+1]=x_t;
            }
            else{
                if(r>rn()){
                    x[j+1]=x_t;
                }
                else{
                    x[j+1]=x[j];
                }
            }
            counter+=1; //计数器 +1
        }

        // 重要抽样法计算积分值
        for(j=0;j<LENTH;j++){
            sum += pow(x[j]-ALPHA*BETA,2.0);
        }

        N = counter;    //获得 N 值，为总步数
        m = N*K;        //获得 m 值，为热化步数
        I[i]=(sum/(double)(N-m)); // 得到积分值
        accurate[i]=abs_d((I[i]-ANSWER)*1.0/I[i]); //计算相对误差 delta I

        printf("gamma No.%d\n",i); // 输出为第几个gamma值的循环
    }

    // 将结果写入文件中
    FILE *fp =NULL;
    fp = fopen(FILE_ACCURATE_1,"w");
    
    // 第一列为gamma值 ；第二列为 积分值 I；第三列为 相对误差值
    for(i=0;i<G_LEN;i++){
        fprintf(fp,"%3.4lf,%3.4lf,%3.4lf\n",gamma[i],I[i],accurate[i]);
    }
    fclose(fp);

}

// Metropolis-Hasting方法的主程序，包括了计算效率的部分
void mc_efficient(){
    int i,N,counter,index;  //一些计数器和 N值
    long j;                 //长整数型计数器
    double *x;              //用于储存Markov Chain的指针数组
    double *I;              //用于储存积分值的指针数组
    double *accurate;       //用于储存相对误差 delta I = Delta I / I的指针数组
    double step;            //等分gamma时用的等差数列的公差
    double x_t;             //产生新的x值的试探值 x_t
    double r;               //Metropolis-Hasting 法中的判别式 r
    double sum;             //计算积分 I 时的加法寄存器 sum
    double m;               //热化步数 m
    double gamma;           //gamma值
    double len;             //Markov Chain的链长

    gamma = 1.0568;         //gamma = 1.0568(此处积分相对误差最小)
    I = (double *)malloc(G_LEN*sizeof(double));         //初始化数组 I
    accurate = (double *)malloc(G_LEN*sizeof(double));  //初始化数组 accurate

    index = 8;  //N的最大值为 10^index = 10^8

    //对于不同数量级的N值，分别进行Metropolis法进行计算
    for(i=1;i<=index;i++){

        N = 0;          //计数器的初始化
        counter=0;      //计数器的初始化
        sum=0;          //计数器的初始化

        len = pow((double)10,(double)i);    //求出链长
        x = (double *)malloc((len)*sizeof(double)); //初始化数组 x
        x[0] = X_0;     //初始点 x_0 = 1 或其他值

        //Metropolis-Hasting法主程序
        for(j=0;j<len;j++){
            x_t = -log(rn());           // [0,infty]生成随机数
            r = func(x[j],x_t,gamma);   //计算判别式 r
            //下面进行一些逻辑运算
            if(r>1){
                x[j+1]=x_t;
            }
            else{
                if(r>rn()){
                    x[j+1]=x_t;
                }
                else{
                    x[j+1]=x[j];
                }
            }
            counter+=1; // 计数器 +1
        }

        // 重要抽样法计算积分值
        for(j=0;j<len;j++){
            sum += pow(x[j]-ALPHA*BETA,2.0);
        }

        N = counter;    //获得 N 值，为总步数
        m = N*K;        //获得 m 值，为热化步数
        I[i]=(sum/(double)(N-m)); // 得到积分值
        accurate[i]=abs_d((I[i]-ANSWER)*1.0/I[i]); //计算相对误差 delta I

        printf("N=10^(%d)\n",i); // 输出为第几个gamma值的循环
    }

    // 将结果写入文件中
    FILE *fp =NULL;
    fp = fopen(FILE_EFFICIENT,"w");

    // 第一列为 N值 ；第二列为 积分值 I；第三列为 相对误差值
    for(i=0;i<index;i++){
        fprintf(fp,"%lf,%3.7lf,%3.7lf\n",pow((double)10,(double)i),I[i],accurate[i]);
    }
    fclose(fp);
}

//调用程序进行的主函数
void main(){
    // Metropolis-Hasting方法主程序，包括了计算精度的部分
    mc_accurate();
    // Metropolis-Hasting方法主程序，包括了计算效率的部分
    //mc_efficient();
}
