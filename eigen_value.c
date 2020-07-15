#include<stdio.h>
#define delta 1e-100

int N;
double A[10000][10000],L[10000][10000],R[10000][10000],Ans[10000][10000];


void reset(double B[10000][10000],int flg){
    int i,j;
    for(i=0;i<N;++i){
        for(j=0;j<N;j++){
            if(i==j)B[i][j]=flg;
            else B[i][j]=0;
        }
    }
}

void LU(double D[10000][10000]){
    int i,j,k;
    reset(L,1);
    reset(R,0);

    //LR分解(LU分解)をする
    for(j=0;j<N;++j){
        for(i=0;i<N;++i){
            if(i<=j){
                //R成分
                R[i][j]=D[i][j];
                for(k=0;k<i;++k){
                    R[i][j]-=L[i][k]*R[k][j];
                }
            }else{
                //L成分
                L[i][j]=D[i][j];
                for(k=0;k<j;++k){
                    L[i][j]-=L[i][k]*R[k][j];
                }
                L[i][j]/=R[j][j];
            }
        }
    }
}

int main(){
    int i,j,k,flg=0,steps=0;
    //もろもろ標準入力する
    scanf("%d",&N);
    for(i=0;i<N;++i)for(j=0;j<N;++j)scanf("%lf",&A[i][j]);

    //下三角成分がすべてゼロ近傍になるまでループ
    while(1){
        //LU分解
        LU(A);

        //行列の掛け算
        for(i=0;i<N;i++){
            for(j=0;j<N;j++){
                A[i][j]=0;
                for(k=0;k<N;k++)A[i][j]+=R[i][k]*L[k][j];
            }
        }

        //収束したか判定
        flg=1;
        for(i=0;i<N;i++){
            for(j=0;j<i;j++){
                if(A[i][j]<-1*delta || A[i][j]>delta){
                    flg=0;
                    break;
                }
            }
        }
        if(flg==1)break;
    }

    printf("\n Eigen values\n");
    for(i=0;i<N;i++){
        printf("%d: %lf\n",i+1,A[i][i]);
    }
}