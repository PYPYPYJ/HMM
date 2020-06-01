#pragma once
#include "nrutil.h"

typedef struct
{
	int N;			//��״̬��
	int M;			//�۲�״̬��
	double *pi;		//��ʼ����ģ��
	double **A;		//ת�Ƹ���ģ��
	double **B;		//�������ģ��
} HMM;

//hmm.c
void ReadHMM(HMM *hmm, FILE *fp);		//���ļ���ȡ�Ѵ��ڵ�hmmģ��
void InitHMM(HMM *hmm, int N, int M);		//�����ʼ��hmmģ��
void freeHMM(HMM *hmm);
void printHMM(HMM *hmm);
void ReadSequence(FILE *fp, int *T, int **O);	//���ļ���ȡ�Ѵ��ڵĹ۲�����

//forward.c			ǰ���㷨��alpha����
void Forward(HMM *hmm, int T, int *O, double **alpha);
double ForwardWithScale(HMM *hmm, int T, int *O, double **alpha, double *scale);
//backward.c		�����㷨��beta����
void Backward(HMM *hmm, int T, int *O, double **beta);
void BackwardWithScale(HMM *hmm, int T, int *O, double **beta, double *scale);
//baumwelch.c		baumwelch�㷨�������ֱ������
int BaumWelch(HMM *hmm, int T, int *O, double *logprobinit, double *logprobfinal);
void ComputeGamma(HMM *hmm, int T, double **alpha, double **beta, double **gamma);
void ComputeXi(HMM *hmm, int T, int *O, double **alpha, double **beta, double ***Xi);
void freeXi(double ***Xi, int T, int N);
void printABG(int N, int T, double **alpha, double **beta, double **gamma);

//hmmrand.c		���úͻ�ȡ�����
int getseed(void);
void setseed(int seed);
double getrandom(void);