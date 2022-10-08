#include<stdio.h>
#include <stdlib.h>
#include<math.h>
int nclusters =0;
int *cluster_ass;
int LSH(int dim, int ndata, double *data,int m, double W, double **h, double *b,int *cluster_start, int *cluster_size, int **cluster_hashval);
int search_LSH(int dim, int ndata, double *data,int m, double W, double **h, double *b,int nclusters, int *cluster_start, int *cluster_size, int **cluster_hashval,double *query_pt, double *result_pt) ;

int LSH(int dim, int ndata, double *data,int m, double W, double **h, double *b,int *cluster_start, int *cluster_size, int **cluster_hashval){
    int *cluster_assign = malloc((ndata)*sizeof(int));
    //iterate all points
    for(int i=0;i<ndata*dim;i=i+dim){
        int hashval[m];
        //finding hashvalue to each point
        for(int j=0;j<m;j++){
            double res = (double)0;
            int p=0;
            for(int k=i;k<i+dim;k++){
                res = res + data[k]*h[j][p];
                p=p+1;
            }
            hashval[j] = (int)floor(((res-b[j])/W));
        }
        //adding the hashvalues to cluster_hashval if 
        if(nclusters == 0){
            cluster_hashval[nclusters] = (int*)malloc(m*sizeof(int));
            for(int c=0;c<m;c++){
                cluster_hashval[nclusters][c] = hashval[c];
            }
            cluster_assign[i/dim] = nclusters;
            nclusters = nclusters + 1;
        }
        else{
            int isSame = 0;
            for(int q=0;q<nclusters;q++){
                int temp=0;
                for(int s=0;s<m;s++){
                    if(hashval[s] == cluster_hashval[q][s]){
                        temp = temp+1;
                    }
                }
                if(temp == m){
                    isSame = 1;
                    cluster_assign[i/dim] = q;
                    break;
                }
            }
            if(!isSame){
                cluster_hashval[nclusters] = (int*)malloc(m*sizeof(int));
                for(int r=0;r<m;r++){
                    cluster_hashval[nclusters][r] = hashval[r];
                }
                cluster_assign[i/dim] = nclusters;
                nclusters = nclusters + 1;
            }
        }
    }
    cluster_ass = cluster_assign;
    return nclusters;
}
int search_LSH(int dim, int ndata, double *data,int m, double W, double **h, double *b,int nclusters, int *cluster_start, int *cluster_size, int **cluster_hashval,double *query_pt, double *result_pt){
    //calculating hashvalue of query point
    int resHashVal[m];
    for(int i=0;i<m;i++){
        double res = (double)0;
        for(int j=0;j<dim;j++){
            res = res + h[i][j]*query_pt[j];
        }
        resHashVal[i] = (int)floor(((res-b[i])/W));
    }
    //search the cluster
    int cluster;
    int isClusterFound = 0;
    for(int i=0;i<nclusters;i++){
        int temp=0;
        for(int j=0;j<m;j++){
            if(resHashVal[j] == cluster_hashval[i][j]){
                temp = temp+1;
            }
        }
        if(temp == m){
            cluster = i;
            isClusterFound = 1;
            break;
        }
    }
    //the closest points in the cluster are
    if(isClusterFound){
        printf("The query point belongs to the cluster %d ", cluster);
        printf("\nThe points in the cluster are:-\n");
        for(int i=0;i<ndata;i++){
            if(cluster_ass[i] == cluster){
                for(int j=i;j<i+dim;j++){
                    printf("%f ",data[j]);
                }
                printf("\n");
            }
        }
    }else{
        printf("\nNo cluster found");
    }
}
int main(){
    int ndata = 100;
    int dim = 16;
    int m = 5;
    double W = 0.5;
    double *b  = malloc((m)*sizeof(double));
    double *data = malloc((ndata*dim)*sizeof(double));
    double **h = (double**)malloc(m*sizeof(double*));
    int **cluster_hashval = (int**)malloc(ndata*sizeof(int*));
    int *cluster_size = malloc(ndata*sizeof(int));
    int *cluster_start = malloc(ndata*sizeof(int));
    double *query_pt = malloc(dim*sizeof(double));
    double *result_pt = malloc(dim*sizeof(double));
    double centroid[dim];
    for(int i=0;i<m;i++){
        h[i] = (double*)malloc(dim*sizeof(double));
    }
    for(int i=0; i< ndata*dim; i++)
    {
        data[i] = ((double)rand()/(RAND_MAX));
    }

    for(int i=0;i<dim;i++){
        query_pt[i] = ((double)rand()/(RAND_MAX));
    }
    //generate random vectors
    for(int i=0;i<m;i++){
        double sum =(double)0;
        for(int j=0;j<dim;j++){
            h[i][j] = 1.0 * rand() / (RAND_MAX / 2) - 1;
            sum = sum + pow(h[i][j],2);
        }
        for(int k=0;k<dim;k++){
            h[i][k] = h[i][k]/sqrt(sum);
        }
    }
    // //print random vectors
    // for(int i=0;i<m;i++){
    //     for(int j=0;j<dim;j++){
    //         printf("%f ",h[i][j]);
    //     }
    //     printf("\n");
    // }

    //find centroid
    for(int i=0;i<dim;i++){
        double res = (double)0;
        for(int j=0;j<ndata*dim;j=j+dim)
        {
            res= res + data[j];
        }
        centroid[i] = res/ndata;
    }
    //setting up array of b
    for(int i=0;i<m;i++){
        double res = (double)0;
        for(int j=0;j<dim;j++){
            res = res + h[i][j]*centroid[j];
        }
        b[i] = res;
    }
    //print b array
    // for(int i=0;i<m;i++){
    //     printf("%f ",b[i]);
    // }
    printf("The number of clusters are -----------%d\n",LSH(dim,ndata,data,m,W,h,b,cluster_start,cluster_start,cluster_hashval));
    search_LSH(dim,ndata,data,m,W,h,b,nclusters,cluster_start,cluster_size,cluster_hashval,query_pt,result_pt);
}