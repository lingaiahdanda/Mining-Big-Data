#include<stdio.h>
#include <stdlib.h>
#include<math.h>
#include<float.h>
int nclusters =0;
int *cluster_ass;
int LSH(int dim, int ndata, double *data,int m, double W, double **h, double *b,int *cluster_start, int *cluster_size, int **cluster_hashval);
int search_LSH(int dim, int ndata, double *data,int m, double W, double **h, double *b,int nclusters, int *cluster_start, int *cluster_size, int **cluster_hashval,double *query_pt, double *result_pt) ;
void closestPoint(int cluster,int ndata,double *query_pt,int dim,double *data);

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
//searching the points
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
        closestPoint(cluster,ndata,query_pt,dim,data);
    }else{
        double minSum = DBL_MAX;
        int closest_cluster;
        for(int i=0;i<nclusters;i++){
            double sum = (double)0;
            for(int j=0;j<m;j++){
                sum = sum + pow((resHashVal[j]-cluster_hashval[i][j]),2);
            }
            if(sum < minSum){
                minSum = sum;
                closest_cluster = i;
            }
        }
        printf("\nNo cluster found\nThe closest cluster is %d\n",closest_cluster);
        closestPoint(closest_cluster,ndata,query_pt,dim,data);
    }
}

//find closest point
void closestPoint(int cluster,int ndata,double *query_pt,int dim,double *data){
    double minSum = DBL_MAX;
        int closest_point;
        printf("The query point belongs to the cluster %d ", cluster);
        printf("\nThe points in the cluster are:-\n");
        for(int i=0;i<ndata;i++){
            if(cluster_ass[i] == cluster){
                double sum=(double)0;
                int iterator = 0;
                for(int j=i;j<i+dim;j++){
                    sum = sum + (pow((query_pt[iterator]-data[j]),2));
                    iterator = iterator+1;
                    printf("%f ",data[j]);
                }
                printf(" Distance from query pt-> %f ",sqrt(sum));
                if(sqrt(sum) < minSum){
                    minSum = sum;
                    closest_point = i;
                }
                printf("\n");
            }
        }
        printf("The closest point is %d and the point is :-\n",closest_point);
        for(int j=closest_point;j<closest_point+dim;j++){
            printf("%f ",data[j]);
        }
        printf("\n");
}
int main(){
    int ndata = 1000000;
    int dim = 16;
    int m = 5;
    double W = 0.3;
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