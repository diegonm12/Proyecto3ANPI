/*
 * main.cpp
 *
 *  Created on: 28 may. 2017
 *      Author: ariel
 */
#include <stdio.h>
#include "iostream"
#include<iomanip>
#include<cmath>
#include <armadillo>
using namespace std;
using namespace arma;


mat temperaturaPunto(int I,int J,int fila,int tam,colvec* colum,double arriba,double abajo,double izquierda,
		double derecha){
	int i=I;
	int j=J;
	mat mat(tam,tam);
	mat.zeros();
	mat.at(i,j)=-4;

	if(i-1==0){

		colum->at(fila,0)=colum->at(fila,0)+izquierda;
		;
	}
	else{
		mat.at(i-1,j)=1;
	}
	if (i+1==(tam)){
		colum->at(fila,0)=colum->at(fila,0)+derecha;
	}
	else{
		mat.at(i+1,j)=1;
	}
	if(j-1==0){
			colum->at(fila,0)=colum->at(fila,0)+abajo;
		}

	else{
			mat.at(i,j-1)=1;
		}
	if(j+1==tam){
				colum->at(fila,0)=colum->at(fila,0)+arriba;
			}
	else{
				mat.at(i,j+1)=1;
			}
return mat;
}

void matrizToFilaToMatriz(mat parcial,mat* final,int tam,int fila){

	int z=0;
	for(int i=1;i<tam;i++){
		for(int j=1;j<tam;j++){
		final->at(fila,z)=parcial.at(i,j);
		z++;
		}
	}


}


mat  liebmannSolver(int tam,mat matriz)
{

	double error=0.0001,y; //error= máximo error permitido
    int numEcua,flag=0,count=0;
    numEcua=((tam-1)*(tam-1));
    mat marizFinal(numEcua,numEcua+1);           // matriz de ecuacion
    marizFinal=matriz;
    mat x(1,numEcua);                //matriz con la ecuaciones
    x.zeros();


    for (int i=0;i<numEcua;i++)                    //matriz diagonalmente dominante
        for (int k=i+1;k<numEcua;k++)
            if (abs(marizFinal.at(i,i))<abs(marizFinal.at(k,i)))
                for (int j=0;j<=numEcua;j++)
                {
                    int temp=marizFinal.at(i,j);
                    marizFinal.at(i,j)=marizFinal.at(k,j);
                    marizFinal.at(k,j)=temp;
                }

    do                            //cálculo de temperaturas
    {

        for (int i=0;i<numEcua;i++)
        {
            y=x.at(0,i);
            x.at(0,i)=marizFinal.at(i,numEcua);
            for (int j=0;j<numEcua;j++)
            {
                if (j!=i)
                	x.at(0,i)=x.at(0,i)-marizFinal.at(i,j)*x.at(j);
            }
            x.at(0,i)=x.at(0,i)/marizFinal.at(i,i);
            if (abs(x.at(0,i)-y)<=error)            //Comparación con el valor anterior
                flag++;

        }

        count++;
    }while(flag<numEcua);
    /*cout<<"\n vector solución:\n";
    for (int i=0;i<numEcua;i++)
        cout<<"x"<<i<<" = "<<x.at(0,i)<<endl;*/
    return x;
}

mat temperatura(int arriba,int abajo, int izquierda,int derecha,int tam){
	mat matriz((tam-1)*(tam-1),(tam-1)*(tam-1));
	matriz.zeros();
	colvec column((tam-1)*(tam-1),1);
	column.zeros();

	int z=0;
	for(int i=1;i<=tam-1;i++){
			for(int j=1;j<=tam-1;j++){
				matrizToFilaToMatriz(temperaturaPunto(j,i,z,tam,&column,arriba,abajo,izquierda,derecha),&matriz,tam,z);
			z++;

			}
		}


     matriz.insert_cols(((tam-1)*(tam-1)),-1*column);




    return matriz;
}

mat rowToMatrix(mat resp,int tam){
	mat matriz(tam-1,tam-1);
	int z=0;
	for(int i=tam-2;i>=0;i--){
			for(int j=0;j<tam-1;j++){
			matriz.at(i,j)=resp.at(0,z);
			z++;
			}
		}
		return matriz;
}
int main(int argc, char** argv) {
//temperaturaPunto(1,1,&A);
//matrizToFilaToMatriz(temperaturaPunto(1,1,0,4,&column,100,0,75,50),&matriz,4,0);
int tam=4;
cout<<"matriz de temperaturas: \n"<<rowToMatrix(liebmannSolver(tam,temperatura(100,0,75,50,tam)),tam)<<endl;
}




