#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#define PI 3.1415926    //圆周率
#define L 200          //正方形边界的长度
#define K 3             //种子边界的长度
#define HALF_L L/2      //正方形边界长度的一半
#define HALF_K K/2      //种子边界长度的一半
//#define ETA 1           //eta 和击穿有关的一个系数
#define PHI_0 1         //phi_0 种子的电势大小
#define ITER_NUM 100    //Jacobi迭代次数
#define ALPHA 0.01      //迭代运算时的电势的精度要求

#define SAND_R K/2   //sandbox计数法 起始长度
#define STEP K/2     //sandbox计数法 间隔步长
#define SAND_N (HALF_L)/(STEP)  //sandbox计数法 采样点总步长数
#define THRUSTER 0.01   //sandbox阈值判断

#define TOL 40000       //预期总模拟点数

#define FILE_DATA "./dbm_3_0_5_200.csv"
#define FILE_SAND "./dbm_sand_3_0_5_200.csv"

//随机数产生器
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
    return x;//返回随机数
}

    // 自定义double型的绝对值函数
double abs_d(double a){
    if(a>0){
        return a;
    }
    else{
        return -a;
    }
}

//用于检查格点有无被占有的函数，若被占有，则返回1.
int check(int **dot,int x,int y){
    if(dot[x][y]==1){ //dot[x][y]代表第x行j列的占有状态
        return 1;
    }
    else{
        return 0;
    }
}

//重置电势状态函数，将占有部分记为1,边界部分记为0，其余部分记为0.5
void refresh_phi(int **dot,double **phi){
    int i,j;
    for(i=1;i<L-1;i++){
        for(j=1;j<L-1;j++){
            if((check(dot,i,j))){
                phi[i][j] = PHI_0;  //占有部分
            }
            else{
                phi[i][j] = 0.5;    //其余部分
            }
        }
    }
    //边界部分
    for(i=0;i<L;i++){
        phi[i][L-1] = 0;
        phi[i][0] = 0;
        phi[0][i] = 0;
        phi[L-1][i] = 0;
    } 
}


//重置占有状态边缘的函数
void refresh_mark(int **mark){
    int i,j;
    for(i=0;i<L;i++){
        for(j=0;j<L;j++){
            mark[i][j]=0;   //全部记为0
        }
    }
}

//检查有无击穿的函数
void exam(int **mark,int *charge){
    int i,j;
    for(i=1;i<L-1;i++){
        //若占有状态边缘部分接近边界部分，则击穿
        if(mark[1][i]==1 || mark[L-2][i]==1 || mark[i][1]==1 || mark[i][L-2]==1){
            charge[0] = 1;
        }
    }
}

//Jacobi迭代部分
double iteration(int **dot,double **phi){
    int i,j;    //计数器
    double delta;   //求出前后变化误差
    double max;     //求出最大的误差，返回主函数，用于精度判断
    for(i=0;i<L;i++){
        for(j=0;j<L;j++){
            if(!(check(dot,i,j))){
                delta = phi[i][j] - (phi[i+1][j]+phi[i-1][j]+phi[i][j+1]+phi[i][j-1])/4;
                phi[i][j] = (phi[i+1][j]+phi[i-1][j]+phi[i][j+1]+phi[i][j-1])/4;
                if(abs_d(delta)>max){
                    max = abs_d(delta);
                }
            }
        }
    }
    return max;
}

//搜索函数，用于找出占有状态边缘的函数，若于边缘，则mark[i][j]记为1；否则为0
void search_2(int **dot,int **mark){
    int i,j;        //计数器
    int flag=0;     //判断标志
    int count=0;    //计数器
    //于格点内部的一般状态判断
    for(i=2;i<L-2;i++){
        for(j=2;j<L-2;j++){
            if(dot[i][j]==0){
                count=0;
                flag=0;
                if(dot[i+1][j]==1){
                    flag=1;
                    count++;
                }
                if(dot[i-1][j]==1){
                    flag=1;
                    count++;
                }
                if(dot[i][j-1]==1){
                    flag=1;
                    count++;
                }
                if(dot[i][j+1]==1){
                    flag=1;
                    count++;
                }
                if((count!=4) && flag){
                    mark[i][j]=1;
                }
            }
        }
    }
    //靠近边界时的特殊状态判断
    for(i=1;i<L-1;i++){
        if(dot[i][1]==0){
            if(dot[i][2]==1){
                mark[i][1]=1;
            }
        }
        if(dot[i][L-2]==0){
            if(dot[i][L-3]==1){
                mark[i][L-2]=1;
            }
        }
        if(dot[1][i]==0){
            if(dot[2][i]==1){
                mark[1][i]=1;
            }
        }
        if(dot[L-2][i]==0){
            if(dot[L-3][i]==1){
                mark[L-2][i]=1;
            }
        }
    }
}

//sandbox计数法
    //传入占有状态矩阵，和两个储存数组
void sandbox(int **dot,int *sandbox_r,int *sandbox_n){
    int i,j;        //计数器
    int r = SAND_R; //sandbox计数法 起始长度
    int step = STEP;    //sandbox计数法 间隔步长
    int len = 0;    //N截至长度，记录N不再变大时的步长
    long counter=0; //计数器
    double delta=0; //误差记录器
    int flag=0;
    int k=0;        //计数器
    while(r<=HALF_L){
        counter=0;
        //记录于r矩形内的点数
        for(i=HALF_L-r;i<HALF_L+r;i++){
            for(j=HALF_L-r;j<HALF_L+r;j++){
                if(dot[i][j]==1){
                    if(i==L-1 || i==0 || j==L-1 || j==0){
                        flag=1;
                        break;
                    }
                    counter++;
                }
            }
        }
        if(flag){
            len = k;
            break;
        }
        sandbox_r[k]=r;         //储存数据
        sandbox_n[k]=counter;   //储存数据
        r+=step;                //增加r
        ///*
        //判断结束条件
        if(k!=0){
            delta = sandbox_n[k]-sandbox_n[k-1];
            //若差值变化小的时候，则退出循环
            if(delta==0){
                len = k;
                break;
            }
        }
        k++;
        //*/
    }

    //写入文件
    FILE *fp =NULL;
    fp = fopen(FILE_SAND,"w");

    for(i=0;i<len;i++){
        fprintf(fp,"%d,%d\n",sandbox_r[i],sandbox_n[i]);
    }
    fclose(fp);
}

//主函数
void main(){
    float r;        //暂时储存随机数
    long row = L;   //行数
    long col = L;   //列数
    long i;         //计数器
    long j;         //计数器
    int counter=0;  //计数器
    long n=0;       //实际占有点数
    int **dot;      //dot[x][y]代表第x行y列的占有状态
    int **mark;     //mark[x][y]代表占有状态的边缘部分，1代表该部分，0则相反
    double **phi;   //phi[x][y]中储存该点的电势值
    int *charge;    //charge是一个是否击穿的标志
    int flag=0;     //判断器
    double delta;   //Jacobi迭代之后，电势的前后变化精度
    double sum = 0; //求和寄存器
    double sum_2 =0;//求和寄存器
    double eta =0.5;  //eta的值
    double phi_0 = 1;//种子的电势值
    double p=0;     //记录生长的概率值
    int *sandbox_r; //sandbox_r[x]用于记录R的值
    int *sandbox_n; //sandbox_n[x]用于记录N的值
    sandbox_r = (int *)malloc((SAND_N+1)*sizeof(int));  //给数组赋值长度
    sandbox_n = (int *)malloc((SAND_N+1)*sizeof(int));  //给数组赋值长度
    dot = (int **)malloc(row*(sizeof(int*)));   //给数组赋值长度
    mark = (int **)malloc(row*(sizeof(int*)));  //给数组赋值长度
    phi = (double **)malloc(row*(sizeof(double*))); //给数组赋值长度
    charge = (int *)malloc(1*sizeof(int));  //给数组赋值长度
    //给数组赋值长度
    for(i=0;i<row;i++){
        dot[i] = (int *)malloc(col*(sizeof(int)));
        mark[i] = (int *)malloc(col*(sizeof(int)));
        phi[i] = (double *)malloc(col*(sizeof(double))); 
    }

    //initialize
    printf("initialize\n");
    for(i=0;i<L;i++){
        for(j=0;j<L;j++){
            dot[i][j] = 0;
        }
    }
        //给起始种子的部位记为占有状态
    for(i=HALF_L-HALF_K;i<=HALF_L+HALF_K;i++){
        for(j=HALF_L-HALF_K;j<=HALF_L+HALF_K;j++){
            dot[i][j] = 1;
            n+=1;   //记录占有状态总数
        }
    }
    //给边缘电势为0部分赋值
    for(i=0;i<L;i++){
        dot[i][L-1] = 1;
        dot[i][0] = 1;
        dot[0][i] = 1;
        dot[L-1][i] = 1;
    }
    charge[0] = 0;//给判断击穿标志赋值

    //当还未击穿，且总点数小于预期点数的时候有的循环
    while(counter<TOL && !charge[0]){
        //每一次循环前，更新电势状态phi及占有部分的边缘状态mark
        refresh_phi(dot,phi);
        refresh_mark(mark);

        //进行Jacobi迭代
        for(i=0;i<ITER_NUM;i++){
            delta = iteration(dot,phi);
            if(delta< ALPHA){
                break;  //满足精度要求时退出迭代
            }
        }

        //搜索占有状态，找出边缘部分，并判断有无击穿
        search_2(dot,mark);
        exam(dot,charge);

        p = 0;      //概率累计求和器
        r = rn();   //生成一个随机数
        sum = 0;    //总速度之和
        sum_2 =0;   //累计速度求和器

        //若为边缘部分，则求出其速度，并计入总和sum
        for(i=0;i<L;i++){
            for(j=0;j<L;j++){
                if(mark[i][j]==1){
                    sum+= n*pow(abs_d(phi_0-phi[i][j]),eta);
                }
            }
        }

        //重新遍历，并求出累计的概率，当累计概率大于随机数时，则占有状态进行生长
        for(i=0;i<L-0;i++){
            flag=0;
            for(j=0;j<L;j++){
                //判断是否为边缘部分
                if(mark[i][j]==1){
                    sum_2 += n*pow(abs_d(phi_0-phi[i][j]),eta);
                    p = sum_2/sum;  //求出累计概率
                    //累计概率大于随机数时，占有状态生长
                    if(r<p){
                        dot[i][j]=1;
                        flag=1;
                        n+=1;
                        break;
                    }
                }
            }
            if(flag){
                break;  //若完成生长，则开始新的循环
            }
        }
        counter++;  //新增占有状态的计数器+1
        printf("order=%dth\n",n);   //输出现有的状态总数
    }
    printf("%d\n",charge[0]);   //输出是否击穿的标志，若为1，则击穿；若为0，则反之

    //写入文件
    FILE *fp =NULL;
    fp = fopen(FILE_DATA,"w");

    for(i=1;i<L-1;i++){
        for(j=1;j<L-1;j++){
            if(dot[i][j]==1){
                fprintf(fp,"%d,%d\n",i,j);
            }
        }
    }
    fclose(fp);

    //进行Sandbox计数
    printf("calculating...\n");
    sandbox(dot,sandbox_r,sandbox_n);
    printf("done!");
}