#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#define VALUE_TYPE int
typedef struct
{
	VALUE_TYPE *e;
	int *columnindex;
	int *rpos; //���е�һ������Ԫ��λ�ñ�,length=row+1
	int mu, nu, tu; 
}RLSMatrix;
void multiSMatrix(RLSMatrix M,RLSMatrix N,RLSMatrix *Q)
{
	 if(M.nu!=N.mu)
	 {
	 	printf("Matrix %d %d cannot be multiplied\n",M.nu,N.mu);
	 	return;
	 }
	 Q->mu=M.mu;
	 Q->nu=N.nu; 
	 Q->tu=0; 
	 Q->rpos=(int*)malloc((1+M.mu)*sizeof(int));//�Ծ���Q��ʼ��
	 Q->e=(VALUE_TYPE*)malloc(sizeof(VALUE_TYPE)*Q->mu*Q->nu);
	 Q->columnindex=(int*)malloc(sizeof(int)*Q->mu*Q->nu);
	 if(M.tu*N.tu!=0)
	 { //����Q�Ƿ������
	 	VALUE_TYPE ctemp[N.nu];
		 for(int arow=0;arow<M.mu;arow++)
		 { //����M��ÿһ��
			
			for(int num=0;num<N.nu;num++)
				ctemp[num]=0; //��ǰ�и�Ԫ���ۼ�������
			Q->rpos[arow] = Q->tu;
			for(int p = M.rpos[arow]; p < M.rpos[arow+1]; p++) //�Ե�ǰ����ÿһ����Ԫ
			{
				 int brow=M.columnindex[p]; //�ҵ���ӦԪ�ھ���N�е��к�*/
				 
				 for (int q=N.rpos[brow];q<N.rpos[brow+1];q++)
				 {
					 int ccol=N.columnindex[q]; /*�˻�Ԫ���ھ���Q�е��к�*/
					 ctemp[ccol]+=M.e[p]*N.e[q];
				 }/*for q*/
			}//���Q�е�crow��=arow���еķ���Ԫ
			for(int ccol=0; ccol<Q->nu; ++ccol) //ѹ���洢���з���Ԫ
			if(ctemp[ccol])
			{
			Q->columnindex[Q->tu]=ccol;
				Q->e[Q->tu]=ctemp[ccol];
				 
				++Q->tu;
			}//if
		}//for arow
		Q->rpos[M.mu] = Q->tu;
	}//if
	return;
}
int main()
{
	RLSMatrix M;
	RLSMatrix N;
	RLSMatrix Q;
	M.mu=4;
	M.nu=3;
	M.tu=3;
	M.e=(VALUE_TYPE*)malloc(sizeof(VALUE_TYPE)*M.tu);
	M.columnindex=(int*)malloc(sizeof(int)*M.tu);
	//M.data[0].i=0;
	M.columnindex[0]=0;
	M.e[0]=1;
	//M.data[1].i=1;
	M.columnindex[1]=0;
	M.e[1]=1;
	//M.data[2].i=2;
	M.columnindex[2]=1;
	M.e[2]=1;
	M.rpos=(int *)malloc(sizeof(int)*(M.mu+1));
	M.rpos[0]=0;
	M.rpos[1]=1;
	M.rpos[2]=2;
	M.rpos[3]=3;
	M.rpos[4]=3;
	N.mu=3;
	N.nu=2;
	N.tu=2;
	N.e=(VALUE_TYPE*)malloc(sizeof(VALUE_TYPE)*N.tu);
	N.columnindex=(int*)malloc(sizeof(int)*N.tu);
	
	//N.data[0].i=1;
	N.columnindex[0]=0;
	N.e[0]=1;
	//N.data[1].i=2;
	N.columnindex[1]=1;
	N.e[1]=1;
	N.rpos=(int *)malloc(sizeof(int)*(M.mu+1));
	N.rpos[0]=0;
	N.rpos[1]=0;
	N.rpos[2]=1;
	N.rpos[3]=2;
	multiSMatrix(M,N,&Q);
	for(int i=0;i<Q.tu;i++)
	{
		printf("%d  %d\n",Q.columnindex[i],Q.e[i]);
	}
	printf("------------\n");
	for(int i=0;i<Q.mu;i++)
	{
		printf("%d  %d\n",Q.rpos[i],Q.rpos[i+1]);
	}
	RLSMatrix fM,ff;
	fM.mu=2;
	fM.nu=4;
	fM.tu=5;
	fM.e=(VALUE_TYPE*)malloc(sizeof(VALUE_TYPE)*fM.tu);
	fM.columnindex=(int*)malloc(sizeof(int)*fM.tu);
	//M.data[0].i=0;
	fM.columnindex[0]=0;
	fM.e[0]=1;
	//M.data[1].i=1;
	fM.columnindex[1]=1;
	fM.e[1]=1;
	fM.columnindex[2]=3;
	fM.e[2]=1;
	//M.data[1].i=1;
	fM.columnindex[3]=0;
	fM.e[3]=1;
	fM.columnindex[4]=2;
	fM.e[4]=1;
	fM.rpos=(int *)malloc(sizeof(int)*(fM.mu+1));
	fM.rpos[0]=0;
	fM.rpos[1]=3;
	fM.rpos[2]=3;
	multiSMatrix(Q,fM,&ff);
	printf("------------\n");
	for(int i=0;i<ff.tu;i++)
	{
		printf("%d  %d\n",ff.columnindex[i],ff.e[i]);
	}
	printf("------------\n");
	for(int i=0;i<ff.mu;i++)
	{
		printf("%d  %d\n",ff.rpos[i],ff.rpos[i+1]);
	}
}
