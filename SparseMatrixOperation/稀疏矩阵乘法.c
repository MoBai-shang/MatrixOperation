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
	int *rpos; //各行第一个非零元的位置表,length=row+1
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
	 Q->rpos=(int*)malloc((1+M.mu)*sizeof(int));//对矩阵Q初始化
	 Q->e=(VALUE_TYPE*)malloc(sizeof(VALUE_TYPE)*Q->mu*Q->nu);
	 Q->columnindex=(int*)malloc(sizeof(int)*Q->mu*Q->nu);
	 if(M.tu*N.tu!=0)
	 { //矩阵Q是非零矩阵
	 	VALUE_TYPE ctemp[N.nu];
		 for(int arow=0;arow<M.mu;arow++)
		 { //处理M的每一行
			
			for(int num=0;num<N.nu;num++)
				ctemp[num]=0; //当前行各元素累加器清零
			Q->rpos[arow] = Q->tu;
			for(int p = M.rpos[arow]; p < M.rpos[arow+1]; p++) //对当前行中每一非零元
			{
				 int brow=M.columnindex[p]; //找到对应元在矩阵N中的行号*/
				 
				 for (int q=N.rpos[brow];q<N.rpos[brow+1];q++)
				 {
					 int ccol=N.columnindex[q]; /*乘积元素在矩阵Q中的列号*/
					 ctemp[ccol]+=M.e[p]*N.e[q];
				 }/*for q*/
			}//求得Q中第crow（=arow）行的非零元
			for(int ccol=0; ccol<Q->nu; ++ccol) //压缩存储该行非零元
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
