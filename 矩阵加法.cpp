#include <stdlib.h>
#include<stdio.h>
#include <sys/time.h>
#include <omp.h>
#include<math.h>
#include <immintrin.h>
int main(int argc,char **argv)
{
	int n=atoi(argv[1]);
	int repeat=10;
	//create vectors
	int *a=(int *)malloc(sizeof(int)*n);
	int *b=(int *)malloc(sizeof(int)*n);
	int *c=(int*)malloc(sizeof(int)*n);
	for(int i=0;i<n;i++)
	{
		a[i]=1;
		b[i]=2;
	}
	struct timeval t1,t2;
	gettimeofday(&t1,NULL);
	for(int i=0;i<repeat;i++)
	{
		for(int j=0;j<n;j++)
		{
			c[i]=a[i]+b[i];
		}
	}
	gettimeofday(&t2,NULL);
	double time_serial=(t2.tv_sec-t1.tv_sec)*1000+(t2.tv_usec-t1.tv_usec)/1000;
	time_serial/=repeat;
	time_serial/=1000.0;
	printf("C serial takes %f sec\n",time_serial);
	gettimeofday(&t1,NULL);
	for(int i=0;i<repeat;i++)
	{
		#pragma omp parallel for
		for(int j=0;j<n;j++)
		{
			c[i]=a[i]+b[i];
			//c[i]=pow(sqrt(a[i]),sqrt(b[i]));
		}
	}
	gettimeofday(&t2,NULL);
	double time_omp=(t2.tv_sec-t1.tv_sec)*1000+(t2.tv_usec-t1.tv_usec)/1000;
	time_omp/=repeat;
	time_omp/=1000.0;
	printf("C-OpenMP serial takes %f sec",time_omp);
	int loop=n/8;
	gettimeofday(&t1,NULL);
	for(int i=0;i<repeat;i++)
	{
		#pragma omp parallel for
		for(int k=0;k<loop;k++)
		{
			_m256i aavx2=_mm256_loadu_si256((_m256i*)(&a[li*8]));
			_m256i bavx2=_mm256_loadu_si256((_m256i*)(&b[li*8]));
			_m256i cavx2=_mm256_add_epi32(aavx2,bavx2);
			_mm256_storeu_si256((_m256i*)(&c[li*8]),cavx2);
		}
		for(int j=loop*8;j<n;j++)
		{
			c[i]=a[i]+b[i];
			//c[i]=pow(sqrt(a[i]),sqrt(b[i]));
		}
	}
	gettimeofday(&t2,NULL);
	double time_avx2=(t2.tv_sec-t1.tv_sec)*1000+(t2.tv_usec-t1.tv_usec)/1000;
	time_avx2/=repeat;
	time_avx2/=1000.0;
	printf("C-Avx2 takes %f sec",time_avx2);
	free(a);
	free(b);
	free(c);
}
