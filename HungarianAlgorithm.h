#ifndef HUNGARIANALGORITHM_H_INCLUDED
#define HUNGARIANALGORITHM_H_INCLUDED
#include "BaseStation.h"
#include "munkres.h"
#include <bits/stdc++.h>

double **SINR_value= new double*[CELLULAR];
double **SINR_value_cpy = new double*[CELLULAR];

void calculateSINR_normal(BaseStation &base)
{
    int i,j;
    double sinr_c,sinr_d;
    double with_out,with;
    base.initializeSharedStatus();
    for(i=0; i<CELLULAR; i++)
    {
        for(j=0; j<CELLULAR; j++)
        {
            if(i<D2D)
            {
                with = 0;
                sinr_c = base.SINR_dl_c(j,i);
                sinr_d = base.SINR_dl_d(i);
                double x,y;
                if(base.constraintCheck(sinr_c,j,sinr_d,i))
                {
                    with += base.getSumRate(sinr_c);
                    with += base.getSumRate(sinr_d);
                }

                sinr_c = base.SINR_dl_c(j,-1);
                with_out = base.getSumRate(sinr_c)+base.getSumRate(0);
                SINR_value[i][j]=with;
                SINR_value_cpy[i][j]=with;
            }
            else
            {
                sinr_c = base.SINR_dl_c(j,-1);
                with_out = base.getSumRate(sinr_c)+base.getSumRate(0);
                SINR_value[i][j]=with_out;
                SINR_value_cpy[i][j]=with_out;
            }

        }
    }
}

void calculateSINR_graph_based(BaseStation &base)
{
    int i,j;
    double sinr_c,sinr_d;
    double with_out,with;
    base.initializeSharedStatus();
    for(i=0; i<CELLULAR; i++)
    {
        for(j=0; j<CELLULAR; j++)
        {
            if(i<D2D)
            {
                with = 0;
                sinr_c = base.SINR_dl_c(j,i);
                sinr_d = base.SINR_dl_d(i);
                double x,y;
                if(base.constraintCheck(sinr_c,j,sinr_d,i))
                {
                    x= base.getSumRate(sinr_c);
                    y= base.getSumRate(sinr_d);

                }
                double xx = base.SINR_dl_c(j,-1);
                double yy =  base.getSumRate(xx);
                SINR_value[i][j]=x+y-yy;
                SINR_value_cpy[i][j]=x+y-yy;
            }
            else
            {
                with_out = 0;
                SINR_value[i][j]=with_out;
                SINR_value_cpy[i][j]=with_out;
            }
        }
    }
}

void calculateSINR_graph_based_ul(BaseStation &base)
{
    int i,j;
    double sinr_c,sinr_d;
    double with_out,with;
    base.initializeSharedStatus();
    for(i=0; i<CELLULAR; i++)
    {
        for(j=0; j<CELLULAR; j++)
        {
            if(i<D2D)
            {
                double x,y,z,s_0;
                sinr_c = base.SINR_ul_eNB(j,i);
                sinr_d = base.SINR_ul_d(j,i);
                s_0 = base.SINR_ul_eNB(j,-1);
                x = base.getSumRate(sinr_c);
                y = base.getSumRate(sinr_d);
                z = base.getSumRate(s_0);
                SINR_value[i][j]=x+y-z;
                SINR_value_cpy[i][j]=x+y-z;
            }
            else
            {
                with_out = 0;
                SINR_value[i][j]=with_out;
                SINR_value_cpy[i][j]=with_out;
            }
        }
    }
}

void calculateSINR_graph_based_inter(BaseStation &base)
{
    int i,j;
    double sinr_c,sinr_d;
    double with_out,with;
    base.initializeSharedStatus();
    for(i=0; i<CELLULAR; i++)
    {
        for(j=0; j<CELLULAR; j++)
        {
            if(i<D2D)
            {
                double x;
                x=base.getInterference_ul(j,i);
                SINR_value[i][j]=x;
                SINR_value_cpy[i][j]=x;
            }
            else
            {
                with_out = MAX_INTER;
                SINR_value[i][j]=with_out;
                SINR_value_cpy[i][j]=with_out;
            }
        }
    }
}

void calculateSINR_coral_inter(BaseStation &base)
{
    int i,j;
    double sinr_c,sinr_d;
    double with_out,with;
    base.initializeSharedStatus();
    for(i=0; i<CELLULAR; i++)
    {
        for(j=0; j<CELLULAR; j++)
        {
            if(i<D2D)
            {
                with=0;
                double x,y,z,s_0;
                sinr_c = base.SINR_ul_eNB(j,i);
                sinr_d = base.SINR_ul_d(j,i);
                s_0 = base.SINR_ul_eNB(j,-1);
                x = base.getSumRate(sinr_c);
                y = base.getSumRate(sinr_d);
                z = base.getSumRate(s_0);
                with=x+y;
                with_out=z;
                if(with-with_out>=0)
                {
                    x=base.getInterference_ul(j,i);
                    SINR_value[i][j]=x;
                    SINR_value_cpy[i][j]=x;
                }
                else
                {
                    SINR_value[i][j]=MAX_INTER;
                    SINR_value_cpy[i][j]=MAX_INTER;
                }
            }
            else
            {
                with_out = MAX_INTER;
                SINR_value[i][j]=with_out;
                SINR_value_cpy[i][j]=with_out;
            }
        }
    }
}

void calculateSINR_coral_ul(BaseStation &base)
{
    int i,j;
    double sinr_c,sinr_d;
    double with_out,with;
    base.initializeSharedStatus();
    for(i=0; i<CELLULAR; i++)
    {
        for(j=0; j<CELLULAR; j++)
        {
            if(i<D2D)
            {
                with=0;
                double x,y,z,s_0;
                sinr_c = base.SINR_ul_eNB(j,i);
                sinr_d = base.SINR_ul_d(j,i);
                s_0 = base.SINR_ul_eNB(j,-1);
                x = base.getSumRate(sinr_c);
                y = base.getSumRate(sinr_d);
                z = base.getSumRate(s_0);
                with=x+y;
                with_out=z;
                //cout<<with-with_out<<endl;
                if(with-with_out>=0)
                {
                    SINR_value[i][j]=with-with_out;
                    SINR_value_cpy[i][j]=with-with_out;
                }
                else
                {
                    //cout<<with-with_out<<endl;
                    SINR_value[i][j]=0;
                    SINR_value_cpy[i][j]=0;
                }
            }
            else
            {
                double s_0 = base.SINR_ul_eNB(j,-1);
                double z = 0;
                SINR_value[i][j]=z;
                SINR_value_cpy[i][j]=z;
            }
        }
    }
}
void calculateSINR_coral(BaseStation &base)
{
    int i,j;
    double sinr_c,sinr_d;
    double with_out,with;
    base.initializeSharedStatus();
    for(i=0; i<CELLULAR; i++)
    {
        for(j=0; j<CELLULAR; j++)
        {
            if(i<D2D)
            {
                double checkCoral=base.SINR_dt_c_coral(j,i);
                with = 0;
                sinr_c = base.SINR_dl_c(j,i);
                sinr_d = base.SINR_dl_d(i);
                double x,y;
                if(base.constraintCheck(sinr_c,j,sinr_d,i))
                {
                    with += base.getSumRate(sinr_c);
                    with += base.getSumRate(sinr_d);
                }
                double check = base.getSumRate(sinr_d);
                checkCoral = base.getSumRate(checkCoral);
                sinr_c = base.SINR_dl_c(j,-1);
                with_out = base.getSumRate(sinr_c);
                if(check-checkCoral>=0)
                {
                    SINR_value[i][j]=with;
                    SINR_value_cpy[i][j]=with;
                }
                else
                {
                    SINR_value[i][j]=with_out;
                    SINR_value_cpy[i][j]=with_out;
                }
            }
            else
            {
                sinr_c = base.SINR_dl_c(j,-1);
                with_out = base.getSumRate(sinr_c);
                SINR_value[i][j]=with_out;
                SINR_value_cpy[i][j]=with_out;
            }
        }
    }
}

void calculateSINR_coral_multipleAssignment(BaseStation &base)
{
    int i,j,newAssignment,l=0;
    double sinr_c_shared,sinr_c_unshared,sinr_d;
    double _unsharedSumrate,_sharedSumrate;

    for(i=0; i<D2D; i++) ///D2D
    {
        for(j=0; j<CELLULAR; j++)
        {
            sinr_c_shared = base.SINR_dl_c(j,i);
            sinr_c_unshared= base.SINR_dl_c(j,-1);
            sinr_d = base.SINR_dl_d(i);

            _sharedSumrate = base.getSumRate(sinr_c_shared)+base.getSumRate(sinr_d);
            _unsharedSumrate = base.getSumRate(sinr_c_unshared);
            if(base.constraintCheck(sinr_c_shared,j,sinr_d,i) && (_sharedSumrate > _unsharedSumrate))
            {
                SINR_value[i][j]=_sharedSumrate;
                SINR_value_cpy[i][j]=_sharedSumrate;
            }
            else
            {
                SINR_value[i][j]=0;
                SINR_value_cpy[i][j]=0;
            }
        }
    }
    for(i=D2D; i<CELLULAR; i++)
    {
        for(j=0; j<CELLULAR; j++)
        {
            sinr_c_unshared= base.SINR_dl_c(j,-1);
            _unsharedSumrate = base.getSumRate(sinr_c_unshared);
            SINR_value_cpy[i][j]=SINR_value[i][j]=_unsharedSumrate;
        }

    }

}

double print_Weight_MATRIX(BaseStation base)
{
    cout<<endl;
    for(int i=0; i<CELLULAR; i++)
    {
        for(int j=0; j<CELLULAR; j++)
        {
            cout<<SINR_value_cpy[i][j]<<"\t";
        }
        cout<<endl;
    }
    for(int i=0; i<CELLULAR; i++)
    {
        double sinr_c = base.SINR_dl_c(i,-1);
        cout<<base.getSumRate(sinr_c)<<"\t";

    }
    cout<<endl;
    cout<<endl;
    cout<<endl;
}

double runHungarianAlgorithm(BaseStation &base,int algorithm)
{
    for(int i=0; i<2*CELLULAR; i++)
    {
        SINR_value[i]= new double[CELLULAR];
        SINR_value_cpy[i]= new double[CELLULAR];
    }
    int nrows;
    int ncols;
    nrows = CELLULAR;
    ncols = CELLULAR;
    if(algorithm==NORMAL_ALGORITHM)
    {
        calculateSINR_normal(base);
    }
    else if(algorithm==GRAPH_BASED_ALGORITHM)
    {
        calculateSINR_graph_based(base);
    }
    else if(algorithm==CORAL_ALGORITHM)
    {
        calculateSINR_coral(base);
    }
    else if(algorithm==UPLINK_ALGORITHM)
    {
        calculateSINR_graph_based_ul(base);
    }
    else if(algorithm==UPLINK_ALGORITHM_CORAL)
    {
        calculateSINR_coral_ul(base);
    }
    //print_Weight_MATRIX(base);
    Matrix<double> matrix(nrows, ncols);
    double value;
    double in[CELLULAR][CELLULAR];
    double max_value=MIN_INT;
    for ( int row = 0 ; row < nrows ; row++ )
    {
        for ( int col = 0 ; col < ncols ; col++ )
        {
            value = SINR_value[row][col];
            in[row][col]=value;
            if(value>max_value)max_value=value;
            matrix(row,col) = value;
        }
    }
    for ( int row = 0 ; row < nrows ; row++ )
    {
        for ( int col = 0 ; col < ncols ; col++ )
        {
            matrix(row,col) = max_value-matrix(row,col);
        }
    }

    Munkres<double> m;
    m.solve(matrix);

    base.initializeSharedStatus();

    for ( int row = 0 ; row < D2D ; row++ )
    {
        for ( int col = 0 ; col < ncols ; col++ )
        {
            if(algorithm==CORAL_ALGORITHM)
            {
                double sinr_c = base.SINR_dl_c(col,-1);
                double without = base.getSumRate(sinr_c) + 0;
                if(matrix(row,col)==0&&SINR_value_cpy[row][col]!=without)/// saki....here we need to add the condition for SINR check
                {
                    base.sharedStatus[row][col]=1;
                }
                else
                {
                    base.sharedStatus[row][col]=0;
                }
            }
            else if(algorithm==NORMAL_ALGORITHM)
            {
                if(matrix(row,col)==0)/// saki...here we need to add the condition for SINR check
                {
                    base.sharedStatus[row][col]=1;
                }
                else
                {
                    base.sharedStatus[row][col]=0;
                }

            }
        }
    }
    if(algorithm==UPLINK_ALGORITHM || algorithm==UPLINK_ALGORITHM_CORAL)
    {
        return base.calculateSystemSumRate_ul();
    }
    else
    {
        return base.calculateSystemSumRate();
    }

}
double runHungarianAlgorithm_inter(BaseStation &base,int algorithm)
{

    for(int i=0; i<2*CELLULAR; i++)
    {
        SINR_value[i]= new double[CELLULAR];
        SINR_value_cpy[i]= new double[CELLULAR];
    }
    base.initializeSharedStatus();
    int nrows;
    int ncols;
    nrows = CELLULAR;
    ncols = CELLULAR;
    if(algorithm==INTERFERENCE_NORMAL)
    {
        calculateSINR_graph_based_inter(base);
    }
    else if(algorithm==INTERFERENCE_CORAL)
    {
        calculateSINR_coral_inter(base);
    }
    //print_Weight_MATRIX();
    Matrix<double> matrix(nrows, ncols);
    double value;
    double in[CELLULAR][CELLULAR];
    double max_value=MIN_INT;
    for ( int row = 0 ; row < nrows ; row++ )
    {
        for ( int col = 0 ; col < ncols ; col++ )
        {
            value = SINR_value[row][col];
            in[row][col]=value;
            if(value>max_value)max_value=value;
            matrix(row,col) = value;
        }
    }
    for ( int row = 0 ; row < nrows ; row++ )
    {
        for ( int col = 0 ; col < ncols ; col++ )
        {
            //matrix(row,col) = max_value-matrix(row,col);
        }
    }

    Munkres<double> m;
    //print_Weight_MATRIX();
    m.solve(matrix);

    base.initializeSharedStatus();

    for ( int row = 0 ; row < D2D ; row++ )
    {
        for ( int col = 0 ; col < ncols ; col++ )
        {
            //cout<<matrix(row,col)<<" ";
            if(matrix(row,col)==0&&SINR_value_cpy[row][col]!=MAX_INTER)
            {
                base.sharedStatus[row][col]=1;
            }
            else
            {
                base.sharedStatus[row][col]=0;
            }

        }
        //cout<<endl;
    }
    return base.calculateSystemSumRate_ul();

}

double runHungarianAlgorithmMultipleTime(BaseStation &base)
{
    for(int i=0; i<2*CELLULAR; i++)
    {
        SINR_value[i]= new double[CELLULAR];
        SINR_value_cpy[i]= new double[CELLULAR];
    }

    int nrows;
    int ncols;
    int newAssignment;
    char b;
    nrows = CELLULAR;
    ncols = CELLULAR;

    base.initializeSharedStatus();
    calculateSINR_coral_multipleAssignment(base);
    while(1)
    {
        newAssignment = 0; /// check is there new assignment in this iteration.

        Matrix<double> matrix(nrows, ncols);
        double value;
        double in[CELLULAR][CELLULAR];
        double max_value=MIN_INT;
        for ( int row = 0 ; row < nrows ; row++ )
        {
            for ( int col = 0 ; col < ncols ; col++ )
            {
                value = SINR_value[row][col];
                in[row][col]=value;
                if(value>max_value)max_value=value;
                matrix(row,col) = value;
            }
        }
        for ( int row = 0 ; row < nrows ; row++ )
        {
            for ( int col = 0 ; col < ncols ; col++ )
            {
                matrix(row,col) = max_value-matrix(row,col);
            }
        }

        Munkres<double> m;
        m.solve(matrix);

        for ( int row = 0 ; row < D2D ; row++ )
        {
            for ( int col = 0 ; col < ncols ; col++ )
            {
                if(matrix(row,col)==0&&SINR_value_cpy[row][col]!=0)
                {
                    base.sharedStatus[row][col]=1;
                    newAssignment++;
                    for(int x=0; x<CELLULAR; x++)
                    {
                        SINR_value_cpy[x][col] = 0;
                    }
                }
                else
                {
                    //base.sharedStatus[row][col]=0;
                }

            }
        }
        if(newAssignment == 0)
        {
            break;
        }
        else
        {
            for(int i=0; i<CELLULAR; i++)
            {
                for(int j=0; j<CELLULAR; j++)
                {
                    SINR_value[i][j] = SINR_value_cpy [i][j];
                }
            }
        }

    }
    return base.calculateSystemSumRate();
}


double runHungarianPlusGreedyAlgorithm(BaseStation &base)
{
    for(int i=0; i<2*CELLULAR; i++)
    {
        SINR_value[i]= new double[CELLULAR];
        SINR_value_cpy[i]= new double[CELLULAR];
    }
    int nrows;
    int ncols;
    int r,c;
    int newAssignment;
    double maximumValue;
    char b;
    nrows = CELLULAR;
    ncols = CELLULAR;

    base.initializeSharedStatus();
    calculateSINR_coral_multipleAssignment(base);

    newAssignment = 0; /// check is there new assignment in this iteration.

    Matrix<double> matrix(nrows, ncols);
    double value;
    double in[CELLULAR][CELLULAR];
    double max_value=MIN_INT;
    for ( int row = 0 ; row < nrows ; row++ )
    {
        for ( int col = 0 ; col < ncols ; col++ )
        {
            value = SINR_value[row][col];
            in[row][col]=value;
            if(value>max_value)max_value=value;
            matrix(row,col) = value;
        }
    }
    for ( int row = 0 ; row < nrows ; row++ )
    {
        for ( int col = 0 ; col < ncols ; col++ )
        {
            matrix(row,col) = max_value-matrix(row,col);
        }
    }

    Munkres<double> m;
    m.solve(matrix);

    for ( int row = 0 ; row < D2D ; row++ )
    {
        for ( int col = 0 ; col < ncols ; col++ )
        {
            if(matrix(row,col)==0&&SINR_value_cpy[row][col]!=0)
            {
                base.sharedStatus[row][col]=1;
                for(int x=0; x<CELLULAR; x++)
                {
                    SINR_value_cpy[x][col] = 0;
                }
            }
            else
            {
                base.sharedStatus[row][col]=0;
            }

        }
    }
    for(int i=0; i<CELLULAR; i++)
    {
        for(int j=0; j<CELLULAR; j++)
        {
            SINR_value[i][j] = SINR_value_cpy [i][j];
        }
    }
    while(1)
    {
        maximumValue = -1;
        for(int i=0; i<D2D; i++)
        {
            for(int j=0; j<CELLULAR; j++)
            {
                if(SINR_value[i][j] > maximumValue)
                {
                    maximumValue = SINR_value[i][j];
                    r=i;
                    c=j;
                }
            }
        }
        if(maximumValue > 0)
        {
            base.sharedStatus[r][c] = 1;
            /// make all the weight of that cellular 0 as it should not share with any D2D
            for(int k=0; k<D2D; k++)
            {
                SINR_value[k][c]=0;
            }
        }
        else
        {
            break;
        }

    }
    return base.calculateSystemSumRate();
}


double runGREEDYmultipleAssignmentAlgorithm(BaseStation &base)
{
    for(int i=0; i<2*CELLULAR; i++)
    {
        SINR_value[i]= new double[CELLULAR];
        SINR_value_cpy[i]= new double[CELLULAR];
    }
    char a;
    base.initializeSharedStatus();
    calculateSINR_coral_multipleAssignment(base);
    double maximumValue = -1;
    int r,c,bb;
    while(1)
    {

        maximumValue = -1;
        for(int i=0; i<D2D; i++)
        {
            for(int j=0; j<CELLULAR; j++)
            {
                if(SINR_value[i][j] > maximumValue)
                {
                    maximumValue = SINR_value[i][j];
                    r=i;
                    c=j;
                }
            }
        }

        if(maximumValue > 0)
        {
            base.sharedStatus[r][c] = 1;
            /// make all the weight of that cellular 0 as it should not share with any D2D
            for(int k=0; k<D2D; k++)
            {
                SINR_value[k][c]=0;
            }
        }
        else
        {
            break;
        }

    }

    return base.calculateSystemSumRate();
}


#endif
