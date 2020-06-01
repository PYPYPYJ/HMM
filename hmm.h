#pragma once
#include "nrutil.h"

typedef struct
{
	int N;			//��״̬��
	int M;			//�۲�״̬��
	double* pi;		//��ʼ����ģ��
	double** A;		//ת�Ƹ���ģ��
	double** B;		//�������ģ��
} HMM;

//hmm.c
void ReadHMM(HMM* hmm, FILE* fp);		//���ļ���ȡ�Ѵ��ڵ�hmmģ��
void InitHMM(HMM* hmm, int N, int M);		//�����ʼ��hmmģ��
void freeHMM(HMM* hmm);
void printHMM(HMM* hmm);

//forward.c			ǰ���㷨��alpha����
void Forward(HMM* hmm, long T, int* O, double** alpha);
double ForwardWithScale(HMM* hmm, long T, int* O, double** alpha, double* scale);
//backward.c		�����㷨��beta����
void Backward(HMM* hmm, long T, int* O, double** beta);
void BackwardWithScale(HMM* hmm, long T, int* O, double** beta, double* scale);
//baumwelch.c		baumwelch�㷨�������ֱ������
int BaumWelch(HMM* hmm, long T, int* O, double* logprobinit, double* logprobfinal);
void ComputeGamma(HMM* hmm, long T, double** alpha, double** beta, double** gamma);
void ComputeXi(HMM* hmm, long T, int* O, double** alpha, double** beta, double*** Xi);
void freeXi(double*** Xi, long T, int N);
void printABG(int N, long T, double** alpha, double** beta, double** gamma);

//hmmrand.c		���úͻ�ȡ�����
int getseed(void);
void setseed(int seed);
double getrandom(void);