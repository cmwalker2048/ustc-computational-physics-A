#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include <string.h>

#define DELTA_R 0.2         //固定步长
#define PI 3.1415926        // pi 
#define LENTH 100000        //Markov Chain 的链长
#define LAMBDA 0.001        //阈值

// ############# 以下均为输出用的输出文件路径 ################
#define FILE_R "./m_rosenbluth.csv"
#define FILE_OUT "./out.txt"
#define FILE_TEMP "./out_temp.txt"
#define FILE_TEMP_LAMBDA "./out_temp_lambda.txt"
#define FILE_TEMP_DELTA_R "./out_temp_delta_r.txt"


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

// 系统哈密顿量表达式
double hamilton(double x,double y){
    return (-2*(x*x+y*y) + 0.5*(pow(x,4.0)+pow(y,4.0))+0.5*pow(x-y,4.0));
}

// Metropolis-Rosenbluth 方法中的判别式
double R_rosenbluth(double x_1,double y_1,double x_2,double y_2,double beta){
    return exp(-beta*(hamilton(x_2,y_2)-hamilton(x_1,y_1)));
}

// 生成一个0-2pi的随机数
double theta(){
    return 2*PI*rn();
}

// 生成一个随机的步长
double delta(){
    return (rn()-0.5)*DELTA_R;  //其中 DELTA_R 为一固定步长
}

// x^2 的求和器
double x_square(double **cor){
    long j;             //长整数型计数器
    double sum;         //寄存器

    sum = 0;

    for(j=0;j<LENTH;j++){
        sum += pow(cor[j][0],2.0);
    }

    return sum;
}

// y^2 的求和器
double y_square(double **cor){
    long j;             //长整数型计数器
    double sum;         //寄存器

    sum = 0;

    for(j=0;j<LENTH;j++){
        sum += pow(cor[j][1],2.0);
    }

    return sum;
}

// x^2 + y^2 的求和器
double x_y_square(double **cor){
    long j;             //长整数型计数器
    double sum;         //寄存器

    sum = 0;

    for(j=0;j<LENTH;j++){
        sum += (pow(cor[j][0],2.0)+pow(cor[j][1],2.0));
    }

    return sum;
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

// 热化过程判断：若完成热化，返回 1；否则返回 0
int thruster(double h1,double h2){
    // 判断 Hamilnion 的相对改变值有无超过阈值
    if(abs_d(h1-h2)/h1 < LAMBDA){
        return 1;
    }
    else{
        return 0;
    }
}


// Metropolis-Rosenbluth法 主函数
void m_rosenbluth(){
    int i,flag;    //一些计数器和逻辑判断器
    long m;             // 热化步数
    long j;             //长整数型计数器
    double **cor;       //用于储存坐标
    double beta[3];     //用于储存 beta值
    double x_t;         //用于储存 x的试探解
    double y_t;         //用于储存 y的试探解
    double Theta;       //用于寄存随机产生的角度值 theta
    double r;           //用于储存表达式 r
    char filepath[3][100] = {"./m_rosenbluth_1.csv","./m_rosenbluth_2.csv","./m_rosenbluth_3.csv"};
    
    // 以下三个是得到其他结果时的输出路径
    //char filepath[3][100] = {"./m_rosenbluth_1_temp.csv","./m_rosenbluth_2_temp.csv","./m_rosenbluth_3_temp.csv"};
    //char filepath[3][100] = {"./m_rosenbluth_1_temp_lambda.csv","./m_rosenbluth_2_temp_lambda.csv","./m_rosenbluth_3_temp_lambda.csv"};
    //char filepath[3][100] = {"./m_rosenbluth_1_temp_delta_r.csv","./m_rosenbluth_2_temp_delta_r.csv","./m_rosenbluth_3_temp_delta_r.csv"};

    beta[0] = 0.2;      // beta = 0.2
    beta[1] = 1;        // beta = 1
    beta[2] = 5;        // beta = 5
    m = 0;              // 热化步数
    flag = 0;           // 归零逻辑判断器

    // 建立二维数组 cor[][]
    cor = (double **)malloc(LENTH*sizeof(double *));
    for(i=0;i<LENTH;i++){
        cor[i] = (double *)malloc(2*sizeof(double));
    }

    // 初始点的获取
    ///*
    cor[0][0] = 10 - 20 * rn();
    cor[0][1] = 10 - 20 * rn();
    //*/
    // 一个另取的初始点
    /*
    cor[0][0] = 8;
    cor[0][1] = -8;
    */
    printf("Initialize done.\n");

    // 循环取遍不同的 beta 值
    for(i=0;i<3;i++){
        printf("beta = %lf\n",beta[i]);
        // 开始主循环
        for(j=0;j<LENTH;j++){
            //防止数组溢出
            if(j==LENTH-1){
                break;
            }
            // 产生新点
            Theta = theta(); //获得一个随机的角度
            x_t = cor[j][0] + delta()*cos(Theta);
            y_t = cor[j][1] + delta()*sin(Theta);

            // 获得 r 并进行逻辑判断
            r = R_rosenbluth(cor[j][0],cor[j][1],x_t,y_t,beta[i]);
            if(r>1){
                cor[j+1][0] = x_t;
                cor[j+1][1] = y_t;
            }
            else{
                if(rn()<r){
                    cor[j+1][0] = x_t;
                    cor[j+1][1] = y_t;
                }
                else{
                    cor[j+1][0] = cor[j][0];
                    cor[j+1][1] = cor[j][1];
                }
            }

            // 判断是否完成热化 并记录热化步数 m
            if(!flag){
                if(thruster(hamilton(cor[j][0],cor[j][1]),hamilton(cor[j+1][0],cor[j+1][1]))){
                    flag = 1;
                    m = j;
                };
            }

            printf("j=%ld\n",j);
        }

        // 输出结果
        printf("beta = %lf\n",beta[i]);
        printf("average x^2 = %4.4lf\n",(x_square(cor)/(LENTH-m)));
        printf("average y^2 = %4.4lf\n",(y_square(cor)/(LENTH-m)));
        printf("average x^2 + y^2 = %4.4lf\n",(x_y_square(cor)/(LENTH-m)));

        // 于文件中写入结果
        FILE *fp_out = NULL;
        fp_out = fopen(FILE_OUT,"a");

        fprintf(fp_out,"beta = %lf\n",beta[i]);
        fprintf(fp_out,"average x^2 = %4.4lf\n",(x_square(cor)/(LENTH-m)));
        fprintf(fp_out,"average y^2 = %4.4lf\n",(y_square(cor)/(LENTH-m)));
        fprintf(fp_out,"average x^2 + y^2 = %4.4lf\n",(x_y_square(cor)/(LENTH-m)));

        fclose(fp_out);

        // 将原始数据写入文件中
        FILE *fp =NULL;
        fp = fopen(filepath[i],"w");

        // 第一列为x值 ；第二列为 y；
        for(j=0;j<LENTH;j++){
            fprintf(fp,"%2.4lf,%2.4lf\n",cor[j][0],cor[j][1]);
        }
        fclose(fp);
    }

}

//调用程序进行的主函数
void main(){
    // Metropolis-Rosenbluth法 主函数
    m_rosenbluth();
}
