#include <iostream> 
#include <cmath>
#include <limits>
#include <vector>
using namespace std;
class gradConj{

public:
    gradConj(){};
    
    double* NoCSR(const double**, const double*, int);
    double* CSR(const double*, const int*, const int*, int, const double*);

};
int main(){
    //Aqui solo estaba probando con ciertas matrices para comprobar el funcionamiento
    gradConj G;
    int m;
    do{
    cout<<"Ingrese las dimensiones de la matriz simetrica definida positiva A, para aproximar Ax=b; m: "; cin>>m;
    m = abs(m);
    if(m == 0){
        m = -1;
    }
    }while(m == -1);
    double**A = new double*[m];
    double *b = new double[m];
    cout<<"Creando A y b..."<<endl;
    for(int i=0; i<m; i++){
        A[i] = new double[m];
       for(int j=0; j<m; j++){
        cout<<"A["<<i<<"]["<<j<<"]: "; cin>>A[i][j];
       }
    }
    cout<<endl<<endl;
    for(int i=0; i<m; i++){
        cout<<"b["<<i<<"]: "; cin>>b[i];
    }
    cout<<"A y b creadas..."<<endl;
    cout<<"Llamando a NoCSR..."<<endl;
    double*x = G.NoCSR(const_cast<const double**>(A), b, m);
    cout<<"Llamada completa, imprimiendo aproximacion: x = (";
    for(int i=0; i<m-1; i++){ 
        cout<<x[i]<<", ";
    }
    cout<<x[m-1]<<")\nTermino la impresion..."<<endl;
    
    cout<<"Cambiando A a formato CSR..."<<endl;
    vector<double> values;
    vector<int> columnas;
    int* rows = new int[m+1];
    rows[0] = 0;
    for(int i=0; i<m; i++){
        int aux =0;
        for(int j=0; j<m; j++){
            if(A[i][j]!=0){
                values.push_back(A[i][j]);
                columnas.push_back(j);
                aux++;
            }
        }
        rows[i+1] = rows[i] + aux;
    }

    rows[m] = rows[m-1] + 1;
    double *valores = new double[values.size()];
    int*columnaz = new int[values.size()];
    for(int i=0; i<values.size(); i++){
        valores[i] = values[i];
        columnaz[i] = columnas[i]; 

    }
    cout<<"Se ha creado A en formato CSR..."<<endl;
    cout<<"\nFilas\tColumnas\tvalores"<<endl;
    for(int i=0; i<m; i++){
        int inicio = rows[i];
        int fin = rows[i+1];
        for(int j=inicio; j<fin; j++){
            cout<<i<<"\t\t"<<columnaz[j]<<"\t\t"<<valores[j]<<endl;;
        }
    }
        
    cout<<"Llamando a G.CSR..."<<endl;
    double*y = G.CSR(valores, rows, columnaz, m, b);
    cout<<"Llamada completa, imprimiendo aproximacion: y = (";
    for(int i=0; i<m; i++){ 
        cout<<y[i]<<", ";
    }
    for(int i=0; i<m; i++){
        delete A[i];
    }
    delete A;
    delete b;
    delete x;
    delete columnaz;
    delete rows;
    delete valores;
    delete y;
    return 0;
}


double *gradConj::NoCSR(const double **A, const double *b, int m)
{

    double*x = new double[m];
    double*r = new double[m];
    double*p = new double[m];
    double*Ap = new double[m];
    double B = 0.0, B1 = 0.0, nr = 0.0, nr_1=0.0, a, aux, Em;
    int limite = 0;
    Em = numeric_limits<double>::epsilon(); //Epsilon de maquina referente a los numeros tipos double
    for(int i=0; i<m; i++){
        x[i] = 0.0; //Inicio igualando x y p = 0, y r = b 
        p[i] = 0.0;
        r[i] = b[i];
        nr+= pow(r[i], 2);
    }
    nr = sqrt(nr); //La norma de r
    while(true){
        if(nr<Em){
            cout<<"El error se torno menor a epsilon de maquina en la iteracion: "<<limite<<endl;
            break;
        }
        else if(limite>10000){
            cout<<"Se superaron los 10,000 intentos, algo pasa aqui..."<<endl;
            break;
        }
        if(limite>0){ //No estamos en la primera iteracion
            B = pow(nr/nr_1, 2); //Bi = r*r / r_1*r_1 = (norm(r)/norm(r_1))^2
        }
        for(int i=0; i<m; i++){
            p[i] = r[i] + B*p[i];   
        }

        a = 0.0; 
            
        for(int i=0; i<m; i++){ //No se pueden juntar estos for porque es necesario conocer todo p
            aux = 0.0;



            Ap[i] = 0.0; 
            for(int j=0; j<m; j++){
                aux+= p[i]*A[j][i];
                Ap[i] += A[i][j]*p[i]; //El vector Ap se ocupa mas abajo
            }
            a+= p[i]*aux;
        }
        a = pow(nr,2) / a; //a = r*r / p*Ap
        nr_1 = nr; //Guardo la norma actual de r porque en la siguiente iteracion la necesito como norm(r_1)
        nr = 0.0; 
        for(int i=0; i<m; i++){ //No se pueden juntar estos for porque ocupo conocer "a"
            x[i] = x[i] + a*p[i]; //Encuentro la aproximación de x
            r[i] = r[i] - a*Ap[i]; //Calculo el error de la solucion
            nr += pow(r[i], 2);
        }
        nr = sqrt(nr); //Calculo la norma del nuevo vector r para verificar si la aproximacion es buena en la siguiente iteracion
        limite++;
    }
    delete r;
    delete p;
    delete Ap;
    return x;
}

double *gradConj::CSR(const double* values, const int* rows_ptr, const int* col_ind, int m, const double*b)
{

    double*x = new double[m];
    double*r = new double[m];
    double*p = new double[m];
    double*Ap = new double[m];
    double B = 0.0, B1 = 0.0, nr = 0.0, nr_1=0.0, a, aux, Em;
    int limite = 0;
    Em = numeric_limits<double>::epsilon(); //Epsilon de maquina referente a los numeros tipos double
    for(int i=0; i<m; i++){
        x[i] = 0.0; //Inicio igualando x y p = 0, y r = b 
        p[i] = 0.0;
        r[i] = b[i];
        nr+= pow(r[i], 2);
    }
    nr = sqrt(nr); //La norma de r
    while(true){
        if(nr<Em){
            cout<<"El error se torno menor a epsilon de maquina en la iteracion: "<<limite<<endl;
            break;
        }
        else if(limite>10000){
            cout<<"Se superaron los 10,000 intentos, algo pasa aqui..."<<endl;
            break;
        }
        if(limite>0){ //No estamos en la primera iteracion
            B = pow(nr/nr_1, 2); //Bi = r*r / r_1*r_1 = (norm(r)/norm(r_1))^2
        }
        for(int i=0; i<m; i++){
            p[i] = r[i] + B*p[i];   
        }

        a = 0.0; 
        for(int i=0; i<m; i++){ //No se pueden juntar estos for porque es necesario conocer todo p
            //Calcular alpha y Ap
            Ap[i] = 0.0;
            int inicio_fila = rows_ptr[i];
            int fin_fila = rows_ptr[i+1]; 
            for(int j=inicio_fila; j<fin_fila; j++){ //Calcular Ap[i], A en formato CSR
                int indice = col_ind[j];
                double valor = values[j];
                Ap[i]+= valor*p[j];
            }
            a+= p[i]*Ap[i]; //p* Ap
        }
        a = pow(nr,2) / a; //a = r*r / p*Ap
        nr_1 = nr; //Guardo la norma actual de r porque en la siguiente iteracion la necesito como norm(r_1)
        nr = 0.0; 
        for(int i=0; i<m; i++){ //No se pueden juntar estos for porque ocupo conocer "a"
            x[i] = x[i] + a*p[i]; //Encuentro la aproximación de x
            r[i] = r[i] - a*Ap[i]; //Calculo el error de la solucion
            nr += pow(r[i], 2);
        }
        nr = sqrt(nr); //Calculo la norma del nuevo vector r para verificar si la aproximacion es buena en la siguiente iteracion
        limite++;
    }
    delete r;
    delete p;
    delete Ap;
    return x;
}
