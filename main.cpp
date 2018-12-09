#include<bits/stdc++.h>
#include "munkres.h"
using namespace std;

int main()
{
    srand(time(NULL));
    int n,m; ///size of the graph n*m
    int instances=20; /// total number of instances

    cout<<"Enter n:";
    cin>>n;
    cout<<"Enter m:";
    cin>>m;
    ofstream output;
    output.open("data.csv"); ///output file. First column represents the cr, second and third column represents the total weight achieved by greedy and OPT respectively. Each row represents each instance.
    output<<"cr,Greedy,OPT"<<endl;
    for(int x=0; x<instances; x++)
    {
        int w=5000;
        int weight[500][500];
        int greedy[500][500];
        Matrix<double> matrix(n, n);
        double value;
        double max_value=-1;

        ///assigning initial weight 0 to n*n matrix
        for(int i=0; i<n; i++)
        {
            for(int j=0; j<n; j++)
            {
                weight[i][j]=0;
                greedy[i][j]=0;
            }
        }

        ///assigning random weight
        for(int i=0; i<m; i++)
        {
            for(int j=0; j<n; j++)
            {
                weight[i][j]=abs(rand()%(w+1));
            }
        }

        ///hungarian algorithm starts
        for ( int row = 0 ; row < n ; row++ )
        {
            for ( int col = 0 ; col < n ; col++ )
            {
                value = weight[row][col];
                if(value>max_value)
                    max_value=value;
                matrix(row,col) = value;
            }
        }
        for ( int row = 0 ; row < n ; row++ )
        {
            for ( int col = 0 ; col < n ; col++ )
            {
                matrix(row,col) = max_value-matrix(row,col);
            }
        }

        Munkres<double> opt;
        opt.solve(matrix);
        double result_optimal=0;
        for(int i=0; i<m; i++)
        {
            for(int j=0; j<n; j++)
            {
                if(matrix(i,j)==0)
                {
                    result_optimal+=weight[i][j];
                }
            }
        }


        ///greedy algorithm starts
        double result_greedy = 0;
        for(int i=0; i<m; i++)
        {
            int max_n = -1;
            int index = -1;
            for(int j=0; j<n; j++)
            {
                if(greedy[i][j]==0)
                {
                    if(weight[i][j]>max_n)
                    {
                        max_n = weight[i][j];
                        index = j;
                    }
                }
            }
            for(int a=0; a<m; a++)
            {
                greedy[a][index]=-1;
            }
            greedy[i][index]=1;
        }
        for(int i=0; i<m; i++)
        {
            for(int j=0; j<n; j++)
            {
                if(greedy[i][j]==1)
                {
                    result_greedy+=weight[i][j];
                }
            }
        }
        double cr = result_optimal/result_greedy;
        cout<<"CR : "<<cr<<endl;
        output<<cr<<","<<result_greedy<<","<<result_optimal<<endl;
    }

    output.close();

    return 0;

}
