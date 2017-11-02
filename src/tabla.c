//#include "tablas.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double* LUT_exp(int N){
  double step = 3.0/(N-1);
  int m= floor(2.5/step);
  double r=0;
  double* res=malloc(N*sizeof(double));
  double E = exp(-25.0/4)/11;
  res[0] = 1;
  for(int i=1;i<m;i++){
    r = r+step;
    res[i] = exp(-r*r);
  }
  for(int i=m;i<N;i++){
    r = r+step;
    res[i] = -E*(r-3.0)*(6*(r-3.0)+25);  // En r=2.5 vale Fo*(0.5)² = -Fo/4
  }
  return res;
}

int guardar_tablas(double *LUT, int N){
  FILE* fp = fopen("tabla.txt", "w");
  fprintf(fp, "%d\n", N);   // Guardo la longitud de la tabla en el primer renglon
  for(int i=0;i<N;i++){
    fprintf(fp, "%f ", LUT[i]);  // Guardo la tabla con potenciales separadas por espacios
  }
  fprintf(fp, "\n"); // Salto de linea para que sepa que termino la tabla
  fclose(fp);
  return 0;
}

int leer_tablas(double **LUT){
  FILE* fp = fopen("tabla.txt","r");
  int N,k=0;
  k=fscanf(fp,"%d\n",&N); // Leo la longitud de mi tabla
  double* tabP = malloc(N*sizeof(double));
  double aux;
  for(int i=0;i<N;i++){
    k=k+fscanf(fp,"%lg",&aux);
    tabP[i]=aux;
  }
  fclose(fp);
  *LUT = tabP;
  if (k==N+1){
    return N;
  }else{ // Si no lei 2*N elementos, algo salió mal y devuelvo 0
    printf("Error al leer tabla: %d en lugar de %d\n",k,2*N+1);
    return 0;
  }
}

int main(int argc, char **argv){
  if (argc==2){
    int N;
    sscanf(argv[1],"%d", &N);
    printf("1\n");
    double* LUT1 = LUT_exp(N);
    guardar_tablas(LUT1,N);
    free(LUT1);
  }else{
    printf("Error: Introduzca cantidad de puntos\n");
  }
  return 0;
}
