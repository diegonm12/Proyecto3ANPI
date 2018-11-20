/*
 * liebmann.cpp
 *
 *  Created on: 28 may. 2017
 *      Author: ariel
 */
#include<iostream>
#include<iomanip>
#include<cmath>
#include <armadillo>
using namespace arma;
using namespace std;
int liebmannSolver2(int tam,mat matriz)
{
	double error=0.0001,y; //error= máximo error permitido
	cout.precision(4);
    cout.setf(ios::fixed);
    int n,flag=0,count=0;
    //n=((tam-1)*(tam-1));
    n=3;
    mat a(n,n+1);           // matriz de ecuacion
    a=matriz;
    mat x(1,n);                //matriz con la ecuaciones
    x.zeros();


    for (int i=0;i<n;i++)                    //Pivotisation(partial) to make the equations diagonally dominant
        for (int k=i+1;k<n;k++)
            if (abs(a.at(i,i))<abs(a.at(k,i)))
                for (int j=0;j<=n;j++)
                {
                    int temp=a.at(i,j);
                    a.at(i,j)=a.at(k,j);
                    a.at(k,j)=temp;
                }
    cout<<"Iter"<<setw(10);
    for(int i=0;i<n;i++)
        cout<<"x"<<i<<setw(18);
    cout<<"\n----------------------------------------------------------------------";
    do                            //Perform iterations to calculate x1,x2,...xn
    {
        cout<<"\n"<<count+1<<"."<<setw(16);
        for (int i=0;i<n;i++)                //Loop that calculates x1,x2,...xn
        {
            y=x.at(0,i);
            x.at(0,i)=a.at(i,n);
            for (int j=0;j<n;j++)
            {
                if (j!=i)
                	x.at(0,i)=x.at(0,i)-a.at(i,j)*x.at(j);
            }
            x.at(0,i)=x.at(0,i)/a.at(i,i);
            if (abs(x.at(0,i)-y)<=error)            //Compare the ne value with the last value
                flag++;
            cout<<x.at(0,i)<<setw(18);
        }
        cout<<"\n";
        count++;
    }while(flag<n);                        //If the values of all the variables don't differ from their previious values with error more than eps then flag must be n and hence stop the loop

    cout<<"\n The solution is as follows:\n";
    for (int i=0;i<n;i++)
        cout<<"x"<<i<<" = "<<x.at(0,i)<<endl;        //Print the contents of x[]
    return 0;
}



