#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>
#include <sys/time.h>
#define VALUE_TYPE int
#define nthreads 16
#define VALUE_TYPE int
typedef struct
{
	VALUE_TYPE *e;
	int *columnindex;
	int *rpos; //各行第一个非零元的位置表,length=row+1
	int mu, nu, tu; 
}RLSMatrix;
void scan(int *array,int n)
{
	int old,New;
	old=array[0];
	array[0]=0;
	for(int i=1;i<n;i++)
	{
		New =array[i];
		array[i]=array[i-1]+old;
		old=New;
	}
 } 
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
	 
	 if(M.tu*N.tu!=0)
	 { //矩阵Q是非零矩阵
	 	
	 	
	 	VALUE_TYPE *values;
	 	int *colindex,*counts;
	 	values=(VALUE_TYPE*)malloc(sizeof(VALUE_TYPE)*Q->mu*Q->nu);
	 	colindex=(int*)malloc(sizeof(int)*Q->mu*Q->nu);
	 	counts=(int*)malloc(sizeof(int)*(1+M.mu));
	 	#pragma omp parallel for
		 for(int arow=0;arow<M.mu;arow++)
		 { //处理M的每一行
		 VALUE_TYPE ctemp[N.nu];
		 int base=arow*N.nu;
		 int base0=base;
		 memset(ctemp, 0, sizeof(VALUE_TYPE) *N.nu);
			//for(int num=0;num<N.nu;num++)
			//	ctemp[num]=0; //当前行各元素累加器清零
			
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
			colindex[base]=ccol;
				values[base]=ctemp[ccol];
				base++;
			}//if
			counts[arow]=base-base0;
			
		}//for arow
		memcpy(Q->rpos, counts, M.mu*sizeof(int));
		scan(Q->rpos,M.mu+1);
		Q->tu=Q->rpos[M.mu];
		Q->e=(VALUE_TYPE*)malloc(sizeof(VALUE_TYPE)*Q->tu);
	 	Q->columnindex=(int*)malloc(sizeof(int)*Q->tu);
	 	#pragma omp parallel for
		for(int arow=0;arow<M.mu;arow++)
		{
			int base=arow*N.nu;
			int qbase=Q->rpos[arow];
			for(int mm=0;mm<counts[arow];mm++)
			{
				Q->columnindex[qbase]=colindex[base];
				Q->e[qbase]=values[base];
				base++;
				qbase++;
			}
		}
		free(values);
	free(counts);
	free(colindex);
	}//if
	
	return;
}
/*
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
	 //omp_set_num_threads(nthreads);
	 	VALUE_TYPE *ctemp;
		ctemp=(VALUE_TYPE*)malloc(sizeof(VALUE_TYPE)*M.nu*N.nu);
			memset(ctemp, 0, sizeof(VALUE_TYPE) *N.nu*M.nu);
		 for(int arow=0;arow<M.mu;arow++)
		 { //处理M的每一行
		 //printf("%d\n",arow);
		 int ctempnums=M.rpos[arow+1]-M.rpos[arow];
		 Q->rpos[arow] = Q->tu;
		 if(ctempnums>0)
		 {
		 	
			#pragma omp parallel for
			for(int p = M.rpos[arow]; p < M.rpos[arow+1]; p++) //对当前行中每一非零元
			{
				int sta=(p-M.rpos[arow])*N.nu;
				int en=sta+N.nu;
				
				 int brow=M.columnindex[p]; //找到对应元在矩阵N中的行号
				 
				 for (int q=N.rpos[brow];q<N.rpos[brow+1];q++)
				 {
					 int ccol=N.columnindex[q]; //乘积元素在矩阵Q中的列号
					 ctemp[sta+ccol]=M.e[p]*N.e[q];
				 }//for q
			}//求得Q中第crow（=arow）行的非零元
			#pragma omp parallel for
			for(int mii=0;mii<N.nu;mii++)
			{
				for(int mjj=1;mjj<ctempnums;mjj++)
				{
					ctemp[mii]+=ctemp[mii+N.nu*mjj];
				}
			}
			
			for(int ccol=0; ccol<Q->nu; ++ccol) //压缩存储该行非零元
			if(ctemp[ccol])
			{
				Q->columnindex[Q->tu]=ccol;
				Q->e[Q->tu]=ctemp[ccol];
				 
				++Q->tu;
			}//if
		 }
			
		}//for arow
		Q->rpos[M.mu] = Q->tu;
		free(ctemp);
	}//if
	return;
}*/
/*
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
	 omp_set_num_threads(nthreads);
	 	VALUE_TYPE *ctemp;
		ctemp=(VALUE_TYPE*)malloc(sizeof(VALUE_TYPE)*M.nu*N.nu);
			
		 for(int arow=0;arow<M.mu;arow++)
		 { //处理M的每一行
		 //printf("%d\n",arow);
		 int ctempnums=M.rpos[arow+1]-M.rpos[arow];
		 Q->rpos[arow] = Q->tu;
		 if(ctempnums>0)
		 {
		 	
			#pragma omp parallel for
			for(int p = M.rpos[arow]; p < M.rpos[arow+1]; p++) //对当前行中每一非零元
			{
				int sta=(p-M.rpos[arow])*N.nu;
				int en=sta+N.nu;
				
				 int brow=M.columnindex[p]; //找到对应元在矩阵N中的行号
				 
				 for (int q=N.rpos[brow];q<N.rpos[brow+1];q++)
				 {
					 int ccol=N.columnindex[q]; //乘积元素在矩阵Q中的列号
					 ctemp[sta+ccol]=M.e[p]*N.e[q];
				 }//for q
			}//求得Q中第crow（=arow）行的非零元
			#pragma omp parallel for
			for(int mii=0;mii<N.nu;mii++)
			{
				for(int mjj=1;mjj<ctempnums;mjj++)
				{
					ctemp[mii]+=ctemp[mii+N.nu*mjj];
				}
			}
			
			for(int ccol=0; ccol<Q->nu; ++ccol) //压缩存储该行非零元
			if(ctemp[ccol])
			{
				Q->columnindex[Q->tu]=ccol;
				Q->e[Q->tu]=ctemp[ccol];
				 
				++Q->tu;
			}//if
		 }
			
		}//for arow
		Q->rpos[M.mu] = Q->tu;
	}//if
	return;
}*/
/*
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
	 VALUE_TYPE *c;
	c= (VALUE_TYPE*)malloc(sizeof(VALUE_TYPE)*Q->mu*Q->nu);
	 Q->rpos=(int*)malloc((1+M.mu)*sizeof(int));//对矩阵Q初始化
	 Q->e=(VALUE_TYPE*)malloc(sizeof(VALUE_TYPE)*Q->mu*Q->nu);
	 Q->columnindex=(int*)malloc(sizeof(int)*Q->mu*Q->nu);
	 if(M.tu*N.tu!=0)
	 { //矩阵Q是非零矩阵
	 omp_set_num_threads(16);
	 	VALUE_TYPE ctemp[N.nu];
	 	#pragma omp parallel for
		 for(int arow=0;arow<M.mu;arow++)
		 { //处理M的每一行
			for(int num=0;num<N.nu;num++)
				ctemp[num]=0; //当前行各元素累加器清零
			Q->rpos[arow] = Q->tu;
			for(int p = M.rpos[arow]; p < M.rpos[arow+1]; p++) //对当前行中每一非零元
			{
				 int brow=M.columnindex[p]; //找到对应元在矩阵N中的行号
				 
				 for (int q=N.rpos[brow];q<N.rpos[brow+1];q++)
				 {
					 int ccol=N.columnindex[q]; //乘积元素在矩阵Q中的列号
					 ctemp[ccol]+=M.e[p]*N.e[q];
				 }//for q
			}//求得Q中第crow（=arow）行的非零元
			int base=arow*N.nu;
			for(int ccol=0; ccol<Q->nu; ++ccol) //压缩存储该行非零元
				c[base+ccol]=ctemp[ccol];
		}//for arow
		for(int arow=0;arow<M.mu;arow++)
		{
			Q->rpos[arow]=Q->tu;
			int base=arow*N.nu;
			
			for(int ccol=0; ccol<Q->nu; ++ccol) 
			{
				
				if(c[base+ccol])
				{
					Q->columnindex[Q->tu]=ccol;
					Q->e[Q->tu]=c[base+ccol];
					Q->tu++;
				}
			}
		}
		Q->rpos[M.mu] = Q->tu;
	}//if
	return;
}*/
/*
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
	 VALUE_TYPE *c;
	c= (VALUE_TYPE*)malloc(sizeof(VALUE_TYPE)*Q->mu*Q->nu);
	 Q->rpos=(int*)malloc((1+M.mu)*sizeof(int));//对矩阵Q初始化
	 Q->e=(VALUE_TYPE*)malloc(sizeof(VALUE_TYPE)*Q->mu*Q->nu);
	 Q->columnindex=(int*)malloc(sizeof(int)*Q->mu*Q->nu);
	 if(M.tu*N.tu!=0)
	 { //矩阵Q是非零矩阵
	 
	 	VALUE_TYPE ctemp[N.nu];
		 omp_set_num_threads(16);
	 	#pragma omp parallel for
		 for(int arow=0;arow<M.mu;arow++)
		 { //处理M的每一行
			for(int num=0;num<N.nu;num++)
				ctemp[num]=0; //当前行各元素累加器清零
			Q->rpos[arow] = Q->tu;
			for(int p = M.rpos[arow]; p < M.rpos[arow+1]; p++) //对当前行中每一非零元
			{
				 int brow=M.columnindex[p]; //找到对应元在矩阵N中的行号
				 for (int q=N.rpos[brow];q<N.rpos[brow+1];q++)
				 {
					 int ccol=N.columnindex[q]; //乘积元素在矩阵Q中的列号
					 ctemp[ccol]+=M.e[p]*N.e[q];
				 }//for q
			}//求得Q中第crow（=arow）行的非零元
			int base=arow*N.nu;
			for(int ccol=0; ccol<N.nu; ++ccol) //压缩存储该行非零元
				c[base+ccol]=ctemp[ccol];
		}//for arow
		for(int arow=0;arow<M.mu;arow++)
		{
			Q->rpos[arow]=Q->tu;
			int base=arow*N.nu;
			
			for(int ccol=0; ccol<N.nu; ++ccol) 
			{
				
				if(c[base+ccol])
				{
					Q->columnindex[Q->tu]=ccol;
					Q->e[Q->tu]=c[base+ccol];
					Q->tu++;
				}
			}
		}
		Q->rpos[M.mu] = Q->tu;
	}//if
	return;
}*/
int main()
{
	omp_set_num_threads(nthreads);
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
	RLSMatrix f1,f2;
	f1.mu=4;
	f1.nu=4;
	f1.tu=15;
	f1.e=(VALUE_TYPE*)malloc(sizeof(VALUE_TYPE)*f1.tu);
	f1.columnindex=(int*)malloc(sizeof(int)*f1.tu);
	//M.data[0].i=0;
	int a[]={0,1,2,3,0,1,2,3,0,1,2,3,0,1,2};
	int b[]={1,5,7,3,2,6,6,2,3,7,5,1,4,8,4};
	for(int i=0;i<f1.tu;i++)
	{
		f1.columnindex[i]=a[i];
		f1.e[i]=b[i];
	}
	
	f1.rpos=(int *)malloc(sizeof(int)*(fM.mu+1));
	f1.rpos[0]=0;
	f1.rpos[1]=4;
	f1.rpos[2]=8;
	f1.rpos[3]=12;
	f1.rpos[4]=15;
	multiSMatrix(ff,f1,&f2);
	printf("------------\n");
	for(int i=0;i<f2.tu;i++)
	{
		printf("%d  %d\n",f2.columnindex[i],f2.e[i]);
	}
	printf("------------\n");
	for(int i=0;i<f2.mu;i++)
	{
		printf("%d  %d\n",f2.rpos[i],f2.rpos[i+1]);
	}
}
