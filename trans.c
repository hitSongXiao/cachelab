/* 
 * trans.c - Matrix transpose B = A^T
 *
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[N][M], int B[M][N]);
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 */ 
#include <stdio.h>
#include "cachelab.h"

int is_transpose(int M, int N, int A[N][M], int B[M][N]);

/* 
 * transpose_submit - This is the solution transpose function that you
 *     will be graded on for Part B of the assignment. Do not change
 *     the description string "Transpose submission", as the driver
 *     searches for that string to identify the transpose function to
 *     be graded. 
 */
char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[N][M], int B[M][N])
{
    int i,j,k,l;
    int b;
    int temp1,temp2,temp3,temp4,temp5,temp6;
    int T1;
    if(M == 64)
    {
        b=8;
        for(i = 0; i < N; i +=b)
            for(j = 0; j < M ;j +=b)
            {
                for(k=0;k<4;k++)
                {
                    for(l=0;l<4;l++)
                    {
                        if(i+k!=j+l)
                        {
                            B[j+l][i+k]=A[i+k][j+l];
                            B[j+l][i+k+4]=A[i+k][j+(l+1)%4+4];//错位，使B在变换时减少驱逐次数（减少抖动），B分块右侧的第j+l行存储的为A的第j+(l+1)%4+4列元素，若不偏移应存储第j+l+4列元素
                        }
                        else
                        {
                            T1=l;
                            temp1=A[i+k][j+l];
                            temp2=A[i+k][j+(l+1)%4+4];
                        }
                    }
                    if(i==j)
                    {
                        B[j+T1][i+k]=temp1;
                        B[j+T1][i+k+4]=temp2;
                    }
                }
                for(k=0;k<4;k++)
                {
                    temp1=A[i+4+k][j+4];
                    temp2=A[i+4+k][j+3];
                    B[j+4][i+k]=B[j+3][i+4+k];
                    B[j+4][i+k+4]=temp1;
                    B[j+3][i+4+k]=temp2;
                }
                for(k=0;k<3;k++)
                {
                    if(i==j)
                    {
                        temp4=A[i+k+5][j+5+k];
                        temp3=A[i+k+5][j+k];//若i==j，则在k+1==l时，则读A的第i+4+k行与写B的第j+5+k行之间会发生抖动，应将该过程分离，因此，先读出A的第i+4+k的元素，等待后续写入B
                        temp5=A[i+4+k][j+k];
                        temp6=A[i+4+k][j+k+5];      //若i==j，则在l==k时，读A的第i+4+k与写B的第j+k行之间会发生抖动
                    }
                    for (l = 0; l < 4; l++)
                    {
                        if(l==k&&i==j)
                        {
                            temp1=temp5;
                            temp2=temp6;
                        }
                        else if(k+1==l&&i==j)
                        {
                            temp1=temp3;
                            temp2=temp4;
                        }
                        else
                        {
                            temp1=A[i+4+l][j+k];
                            temp2=A[i+4+l][j+5+k];
                        }
                        B[j+5+k][i+l]=B[j+k][i+4+l];
                        B[j+5+k][i+4+l]=temp2;
                        B[j+k][i+4+l]=temp1;
                    }
                }
            }
    }
    else
    {
        if(M==32)
            b=8;//一个块32字节，可容纳8个整型数，因此块的长和宽为8
        else
            b=23;//根据将矩阵的每一行在cache的映射组数计算出来可以发现，第24行与第1行同时写到cache的同一组上，会发生抖动，因此，分块的边长应为23
        for(i=0;i<N;i+=b)
        {
            for(j=0;j<M;j+=b)
            {
                for(k=i;k<i+b&&k<N;k++)//当访问A[i][j]时，A[i][j]——A[i][j+7]的元素写入cache中，同样，访问B[j][i]时，B[j][i]——B[j][i+7]元素写入cache中，因此2个矩阵遍历8*8方块时，只会产生16次不命中（不考虑抖动）
                {
                    for(l=j;l<j+b&&l<M;l++)
                    {
                        if(k!=l)
                        {
                            B[l][k]=A[k][l];
                        }
                        else//由于cache的规格为32*32字节，A与B的大小为32*32*4字节，因此，A与B的同一行会映射到cache的同一块上，会造成抖动，读取A的k行时，cache被A的第k行写满，当读取B的第l行时，由于B==l，（且A的大小为cache的整数倍）因此A与B的地址映射到同一组中，一组只有一块，因此会覆盖写入的A，当A再次访问下一个元素时，会再次不命中
                        {
                            temp1=A[k][l];
                            T1=l;
                        }
                    }
                    if(i==j)
                    {
                        B[T1][T1]=temp1;//避免抖动的发生，在A访问1行（8个值）后在写入B的对角线元素，可避免同时访问A与B的同一行造成抖动（如果M=61,A大小为4087字节，与cache大小的4倍相差9字节，同样可能发生抖动）
                    }
                }
            }
        }
    }
}

/* 
 * You can define additional transpose functions below. We've defined
 * a simple one below to help you get started. 
 */ 

/* 
 * trans - A simple baseline transpose function, not optimized for the cache.
 */
char trans_desc[] = "Simple row-wise scan transpose";
void trans(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, tmp;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            tmp = A[i][j];
            B[j][i] = tmp;
        }
    }    

}

/*
 * registerFunctions - This function registers your transpose
 *     functions with the driver.  At runtime, the driver will
 *     evaluate each of the registered functions and summarize their
 *     performance. This is a handy way to experiment with different
 *     transpose strategies.
 */

void registerFunctions()
{
    /* Register your solution function */
    registerTransFunction(transpose_submit, transpose_submit_desc); 

    /* Register any additional transpose functions */
    registerTransFunction(trans, trans_desc); 

}

/* 
 * is_transpose - This helper function checks if B is the transpose of
 *     A. You can check the correctness of your transpose by calling
 *     it before returning from the transpose function.
 */
int is_transpose(int M, int N, int A[N][M], int B[M][N])
{
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; ++j) {
            if (A[i][j] != B[j][i]) {
                return 0;
            }
        }
    }
    return 1;
}

