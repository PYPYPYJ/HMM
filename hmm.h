#pragma once
#include "nrutil.h"

typedef struct
{
	int N;			//隐状态数
	int M;			//观测状态数
	double* pi;		//初始概率模型
	double** A;		//转移概率模型
	double** B;		//发射概率模型
} HMM;

//hmm.c
void ReadHMM(HMM* hmm, FILE* fp);		//从文件读取已存在的hmm模型
void InitHMM(HMM* hmm, int N, int M);		//随机初始化hmm模型
void freeHMM(HMM* hmm);
void printHMM(HMM* hmm);

//forward.c			前向算法算alpha矩阵
void Forward(HMM* hmm, long T, int* O, double** alpha);
double ForwardWithScale(HMM* hmm, long T, int* O, double** alpha, double* scale);
//backward.c		后向算法算beta矩阵
void Backward(HMM* hmm, long T, int* O, double** beta);
void BackwardWithScale(HMM* hmm, long T, int* O, double** beta, double* scale);
//baumwelch.c		baumwelch算法估算参数直至收敛
int BaumWelch(HMM* hmm, long T, int* O, double* logprobinit, double* logprobfinal);
void ComputeGamma(HMM* hmm, long T, double** alpha, double** beta, double** gamma);
void ComputeXi(HMM* hmm, long T, int* O, double** alpha, double** beta, double*** Xi);
void freeXi(double*** Xi, long T, int N);
void printABG(int N, long T, double** alpha, double** beta, double** gamma);

//hmmrand.c		设置和获取随机数
int getseed(void);
void setseed(int seed);
double getrandom(void);