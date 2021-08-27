#include<stdio.h>
#include<malloc.h>
typedef struct//三元组定义
{
	int i;//行号 
	int j;//列号 
	float e;//元素值 
}Triple; 
typedef struct//三元组顺序表定义 
{
	int mu;//行数值 
	int nu;//列数值
	int tu;//非零元素个数
	Triple *data; 
}TSMatrix;
typedef struct//行逻辑链接的顺序表
{
Triple *data; 
int *rpos; //各行第一个非零元的位置表
int mu, nu, tu; 
}RLSMatrix;
void MatrixCompress(float **A,TSMatrix &A_C,int row,int col)//矩阵压缩 
{
	int i,j,sum=0;
	A_C.mu=row;
	A_C.nu=col;
	A_C.tu=0;
	for(i=0;i<row;i++)
		for(j=0;j<col;j++)
			if(A[i][j])
				sum++;
	A_C.data=(Triple*)malloc(sizeof(Triple)*sum);
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
			if(A[i][j])
			{
				A_C.data[A_C.tu].i=i;
				A_C.data[A_C.tu].j=j;
				A_C.data[A_C.tu].e=A[i][j];
				A_C.tu++;
			}
	}
	printf("STMatrix Compress succed\n");
}
void RLSMatrixCompress(float **A,RLSMatrix &A_C,int row,int col)//矩阵压缩 
{
	int i,j,sum=0;
	A_C.mu=row;
	A_C.nu=col;
	A_C.tu=0;
	for(i=0;i<row;i++)
		for(j=0;j<col;j++)
			if(A[i][j])
				sum++;
	A_C.data=(Triple*)malloc(sizeof(Triple)*sum);
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
			if(A[i][j])
			{
				A_C.data[A_C.tu].i=i;
				A_C.data[A_C.tu].j=j;
				A_C.data[A_C.tu].e=A[i][j];
				A_C.tu++;
			}
	}
	A_C.rpos=(int*)malloc(sizeof(int)*A_C.mu);	
	for(j=0;j<A_C.mu;j++)
		A_C.rpos[j]=-1;
	for(j=0,i=0;j<A_C.tu;j++)
	{
		if(A_C.data[j].i==i)
		{
			A_C.rpos[i]=j;
			i++;
		}
		else if(A_C.data[j].i>i)
		{
			A_C.rpos[i]=j;
			i++;
			j--;
		}
		if(i>=A_C.mu)
			break;
	}
	if(A_C.rpos[A_C.mu-1]<0)
		A_C.rpos[A_C.mu-1]=A_C.tu;
	//printf("RLSMatrix Compress succed\n");
}
void MatrixFastTranspose(TSMatrix M, TSMatrix &T)
{
	int j,col,q,p,*num,*copt;
	T.mu=M.nu;
	T.nu=M.mu;
	T.tu=M.tu;
	if(M.tu)
	{
		//num存储矩阵A中某列的非零元素的个数；
		num=(int*)malloc(sizeof(int)*M.nu);
		//cpot：初值表示矩阵A中某列的第一个非零元素在B中的位置。
		copt=(int*)malloc(sizeof(int)*M.nu);
		//num赋初值
		for(j=0;j<M.nu;j++)
			num[j]=0;
		//num赋值
		for(j=0;j<M.tu;j++)
			num[M.data[j].j]++;
		//copt赋初值 
		copt[0]=0;
		//cpot赋值
		for(j=1;j<M.nu;j++)
			copt[j]=copt[j-1]+num[j-1];	
		for(p=0;p<M.tu;p++)
		{
			col=M.data[p].j;//取得M中第p个元素的列值
			q=copt[col]; //第col列的第一个非零元在T中位置
			//将M存储到T中正确的位置
			T.data[q].i=M.data[p].j;
			T.data[q].j=M.data[p].i;
			T.data[q].e=M.data[p].e;
			copt[col]++;
		}
		printf("\nMatrix Fast Transpose succed\n");
	}
}
void MultiMatrix()//稀疏矩阵的乘法运算
{
	float *ctemp;
	RLSMatrix M,N,Q;
	int row,col,i,j;
	int arow,k,p,tp,brow,t,q,ccol,MAXSIZE;
	float **A,**B;
	printf("Input the row number and column number of matrix which you used\nrow,col:");
	scanf("%d,%d",&row,&col);
	printf("Input the matrix in the order with row fistly\n"); 
	A=(float**)malloc(sizeof(float)*row);
	for(i=0;i<row;i++)
	{
		A[i]=(float*)malloc(col*sizeof(float));
		for(j=0;j<col;j++)
			scanf("%5f",A[i]+j);
	}
	RLSMatrixCompress(A,M,row,col);
	printf("Input the row number and column number of matrix which you used\nrow,col:");
	scanf("%d,%d",&row,&col);
	printf("Input the matrix in the order with row fistly\n"); 
	B=(float**)malloc(sizeof(float)*row);
	for(i=0;i<row;i++)
	{
		B[i]=(float*)malloc(col*sizeof(float));
		for(j=0;j<col;j++)
			scanf("%5f",B[i]+j);
	}
	RLSMatrixCompress(B,N,row,col);
	MAXSIZE=M.mu*N.nu;
	ctemp=(float*)malloc(sizeof(float)*M.nu);
	Q.data=(Triple*)malloc(sizeof(Triple)*MAXSIZE);
	if(M.nu!=N.mu)
	 	printf("ERROR:M.nu!=N.mu\n");
	else
	{
		Q.rpos=(int*)malloc(sizeof(int)*M.mu);
		Q.mu=M.mu;Q.nu=N.nu; Q.tu=0; //对矩阵Q初始化
		if(M.tu*N.tu)//矩阵Q是非零矩阵
		{ 
			for(arow=0;arow<M.mu;arow++) //处理M的每一行
			{
				if(M.rpos[arow+1]-M.rpos[arow]<1)//当前行元素全为0 
					continue;
				else
				{
					//当前行各元素累加器清零
					for(k=0;k<M.nu;k++)
						ctemp[k]=0;
						
					Q.rpos[arow] = Q.tu;
					if(arow<M.mu-1) //设置矩阵M每行的上限标号,即每行非零元素的位置范围 
						tp = M.rpos[arow+1]; 
					else
						tp =M.tu;
					for(p = M.rpos[arow]; p < tp; p++) //对当前行中每一非零元
					{
						brow=M.data[p].j; //找到对应元在矩阵N中的行号*/
						if(brow<N.mu-1) //设置矩阵N每行的上限标号,即每行非零元素的位置范围
							t=N.rpos[brow+1];
						 else 
						 	t=N.tu;
						for (q=N.rpos[brow];q<t;q++)//对当前行中每一非零元
						{
							 ccol=N.data[q].j; /*乘积元素在矩阵Q中的列号*/
							 ctemp[ccol]+=M.data[p].e*N.data[q].e;
						}/*for q*/
					}//求得Q中第crow（=arow）行的非零元
					for(ccol=0; ccol<Q.nu; ++ccol) //压缩存储该行非零元
						if(ctemp[ccol])
						{
							if(++Q.tu>MAXSIZE)
							{
								printf("ERROR");
								break;
							}
							Q.data[Q.tu-1].i=arow;
							Q.data[Q.tu-1].j=ccol;	
							Q.data[Q.tu-1].e = ctemp[ccol];
						}//if
				}
				
			}//for arow
		}//if
	printf("Multi Matrix succed\nThis is answer:\n");
	for(i=0;i<Q.tu;i++)
	printf("row,col,value:%4d%4d%8.2f\n",Q.data[i].i+1,Q.data[i].j+1,Q.data[i].e);
	}
}

main()
{
	int row,col,i,j;
	float **A;
	TSMatrix A_C,T;
	RLSMatrix M,N,Q;
	int *num,*copt;
	MultiMatrix();
	printf("Input the row number and column number of matrix which you used\nrow,col:");
	scanf("%d,%d",&row,&col);
	printf("Input the matrix in the order with row fistly\n"); 
	A=(float**)malloc(sizeof(float)*row);
	for(i=0;i<row;i++)
	{
		A[i]=(float*)malloc(col*sizeof(float));
		for(j=0;j<col;j++)
			scanf("%5f",A[i]+j);
	}
	MatrixCompress(A,A_C,row,col);
	for(i=0;i<A_C.tu;i++)
		printf("%8.2f",A_C.data[i].e);
	MatrixFastTranspose(A_C,T);
	for(i=0;i<T.tu;i++)
		printf("%8.2f",T.data[i].e);
	printf("\n");
	
}
