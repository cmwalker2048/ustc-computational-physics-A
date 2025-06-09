#include<math.h>
#include<stdio.h>

// PB21020660 何金铭
// HW1 —— 利用Schrage方法编写随机数子程序

void main(){
    // 16807产生器的参数
    int A = 16807; //参数a
    int B = 0; //参数b
    long M = 2147483647; //参数m
    long i = 23942907; //起始值i_0的大小(种子值)
    // Schrage法的参数
    long Q = 127773; //参数q（除数）
    int R = 2836; //参数r（余数）
    // 默认产生随机数的个数
    int n = 10000;
    // 运算中暂时储存的产生的随机数
    float x;

    //让用户输入想得到的随机数个数
    printf("Please input the number of Random Number you want:");
    scanf("%ld",&n);
    
    FILE *fp = NULL; //文件指针
    fp = fopen("./16807.csv","w"); // 打开文件

    // 主循环
    for(int j=1;j<=n;j=j+1){
        x = (i*1.0) / M; //得到随机数
        //利用Schrage方法取模
        i = A*(i % Q) - R*(i / Q);
        if(i < 0){
            i = i + M; //Schrage法的一个步骤
        }
        fprintf(fp,"%f\n",x); //于文件中写入随机数
    }
    fclose(fp); //关闭文件
}