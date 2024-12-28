#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;
const double eps       = 0.00001;  // epsilon
const double a_         = 0.0;      //
const double b_         = 1.0;
const int N_            = 2000;
const double lambda_a  = 0.0;
const double lambda_b  = 1.0;
const double mu_a  = 1.0;
const double mu_b  = 0.0;
const double psi_a  = 1.0;
const double psi_b  = 1.38294;

double MAX_ITERATIONS  = 1000;

double x_i(int N, int i){
    return a_ + ((b_ - a_)/(N))*i;
}

double p(int N, int i){
    return x_i(N, i) + 1.0;
}

double q(int i){
    return -1.0;
}

double g(int N, int i){
    double x = x_i(N, i);
    return (x*x + 2.0*x + 2.0)/(x + 1.0);
}

struct three_diagonals_matrix {
private:
    double a[10001];
    double b[10001];
    double c[10001];
    double f[10001];
public:
    three_diagonals_matrix(){
       for(int i = 0; i < 10001; i++){
           a[i] = 0;
           b[i] = 0;
           c[i] = 0;
           f[i] = 0;
       }
    }
    void set_values(double ai, double bi, double ci, double fi, int i){
        if (i < 0){
            cout << " i < 0" << endl;
        }
        a[i] = ai;
        b[i] = bi;
        c[i] = ci;
        f[i] = fi;
    }
    double ai(int i){
        return a[i];
    }
    double bi(int i){
        return b[i];
    }
    double ci(int i){
        return c[i];
    }
    double fi(int i){
        return f[i];
    }

};

void solve_three_diagonals_matrix(three_diagonals_matrix* M, double* u, int N)
{
    double temp;
    double alpha[N + 2];//первое поле использовать не буду, ленивый ыыыыы, лень двигать индексы, говнокод, да, но я очень спешу
    double beta[N + 2];
    //double u[N + 1];
    alpha[1] = M->bi(0)/M->ci(0);
    beta[1] = M->fi(0)/M->ci(0);

    for (int i = 1; i < N - 1; i++)
    {
        alpha[i + 1] = M->bi(i)/(M->ci(i) - alpha[i]*M->ai(i));
        beta[i + 1] =(M->fi(i) + beta[i]*M->ai(i))/(M->ci(i) - alpha[i]*M->ai(i));
    }
    beta[N] =(M->fi(N-1) + beta[N-1]*M->ai(N-1))/(M->ci(N-1) - alpha[N-1]*M->ai(N-1));
    temp = M->fi(N-1) ;
    temp = beta[N-1]*M->ai(N-1);
    temp = M->ci(N-1);
    temp = alpha[N-1]*M->ai(N-1);

    u[N-1] = beta[N];
    temp = beta[N];

    for (int i = N - 1; i > 0; i--)
    {
        //temp = *(u + i);
        u[i - 1] = alpha[i]*u[i] + beta[i];
        temp = u[i - 1];
    }
}

void dummy_node_method(double* u, int N) {
    double a,b,c,f;
    double h = (b_ - a_)/N;
    three_diagonals_matrix M;
    //double u[N_ + 1];
    //заполнение трехдиагональной матрицы со страницы 7
    a = 0;
    b = mu_a/(h*h);
    c = mu_a/(h*h) - lambda_a/h + (lambda_a*p(N, 0)) - mu_a*q(0)/2;
    f = -(1/h + p(N,0)/2)*psi_a - mu_a*g(N,0)/2;
    M.set_values(a, b, c,f,0);
    for(int i = 1; i < N - 1; i++){
        a = 1.0/(h*h) - p(N,i)/(2 * h);
        b = 1.0/(h*h) + p(N,i)/(2 * h);
        c = 2.0/(h*h) - q(i);
        f = -g(N,i);
        M.set_values(a, b, c,f,i);
    }
    // N -1 строка
    a = 1.0/(h*h) - p(N,N-1)/(2 * h);
    b = 0;
    c = 2.0/(h*h) - q(N - 1);
    f = (1.0/(h*h) + p(N,N - 1)/(2 * h))*(psi_b/lambda_b) - g(N,N - 1);
    M.set_values(a, b, c,f,N - 1);
    solve_three_diagonals_matrix(&M,u, N );

}

void asymmetric_difference_derivative_method(double* u, int N){
    double a,b,c,f, temp;
    double h = (b_ - a_)/N;
    three_diagonals_matrix M;

    //заполнение трехдиагональной матрицы со страницы 7
    a = 0;
    b = mu_a/(h*h) + (mu_a*p(N,1))/h + mu_a*q(1)/2;
    c = mu_a/(h*h) + ((mu_a*p(N,1)) - lambda_a)/h - lambda_a*p(N,1)/2;
    f = -(1.0/h + p(N,1)/2.0)*psi_a - mu_a*g(N,1)/2.0;
    M.set_values(a, b, c,f,0);
    for(int i = 1; i < N - 1; i++){
        a = 1.0/(h*h) - p(N,i)/(2.0 * h);
        b = 1.0 / (h * h) + p(N,i) / (2.0 * h);
        c = 2.0/(h*h) - q(i);
        f = -g(N,i);
        M.set_values(a, b, c,f,i);
    }
    // N -1 строка
    a = 1.0/(h*h) - p(N,N-1)/(2.0 * h);
    b = 0;
    c = 2.0/(h*h) - q(N - 1);
    f = (1.0/(h*h) + p(N,N - 1)/(2.0 * h))*(psi_b/lambda_b) - g(N,N - 1);
    M.set_values(a, b, c,f,N - 1);

    solve_three_diagonals_matrix(&M, u, N);

}

double infinity_norm(double* u, double* v, int N) {
    double max = 0.0;
    for(int i = 0; i < N; i++) {
        if (fabs(u[i] - v[i]) > max) {
            max = fabs(u[i] - v[i]);
        }
    }
    return max;
}

double one_norm(double* u, double* v, int N) {
    double norm = 0.0;
    double h = (b_ - a_)/N;
    for(int i = 0; i < N; i++) {
        norm += fabs(u[i] - v[i])*h;
    }
    return norm;
}

double two_norm(double* u, double* v, int N) {
    double norm = 0.0;
    double h = (b_ - a_)/N;
    for(int i = 0; i < N; i++) {
        norm += ((u[i] - v[i]))*((u[i] - v[i]))*h;
    }
    return sqrt(norm);
}

//вычисление отношения норм
double n1_n2(int mode, int n1) {
    double u[n1], v[n1], u2[2*n1], v2[2*n1] ;
    double norm1 = 0.0, norm2 = 0.0;
    if(mode == 0) {
        asymmetric_difference_derivative_method(u,n1);
        dummy_node_method(v, n1);
        norm1 = infinity_norm(u,v,n1);
        asymmetric_difference_derivative_method(u2,2*n1);
        dummy_node_method(v2,2*n1);
        norm2 = infinity_norm(u2,v2,2*n1);
        if (norm2 == 0) {
            return 100;
        }
        return norm1/norm2;
    }
    if(mode == 1) {
        asymmetric_difference_derivative_method(u,n1);
        dummy_node_method(v,n1);
        norm1 = one_norm(u,v,n1);
        asymmetric_difference_derivative_method(u2,2*n1);
        dummy_node_method(v2,2*n1);
        norm2 = one_norm(u2,v2,2*n1);
        if (norm2 == 0) {
            return 101;
        }
        return 2*norm1/norm2;
    }
    if(mode == 2) {
        asymmetric_difference_derivative_method(u,n1);
        dummy_node_method(v,n1);
        norm1 = two_norm(u,v,n1);
        asymmetric_difference_derivative_method(u2,2*n1);
        dummy_node_method(v2,2*n1);
        norm2 = two_norm(u2,v2,2*n1);
        if (norm2 == 0) {
            return 102;
        }
        return 2*norm1/norm2;
    }
    return 0;
}

void table_of_n1_n2() {
    int arr_n1[8] = {25,50,100,200,500,1000,2000,5000};
    ofstream res_out;
    res_out.open("result_table.txt");
    res_out << "eps   = " << eps << endl;
    res_out << "+================================================================================================================+" << endl;
    res_out << "|  n1/n2         | 25/50     | 50/100    | 100/200   | 200/400   | 500/10^3  | 5000/10^3 | 2000/4000 | 5000/10^4 |"<< endl;
    res_out << "+================+===========+===========+===========+===========+===========+===========+===========+===========+" << endl;
    res_out << "|" << "infinity        " << "|" << setw(11) << n1_n2(0,25) << "|"<< setw(11) << n1_n2(0,50)<< "|" << setw(11) << n1_n2(0,100)<< "|" ;
    res_out << setw(11) << n1_n2(0,200)<<  "|" << setw(11) << n1_n2(0,500)<<"|" << setw(11) << n1_n2(0,1000)<< "|" << setw(11) << n1_n2(0,2000)<<"|" << setw(11) << n1_n2(0,5000)<<"|"<< endl;
    res_out << "+================+===========+===========+===========+===========+===========+===========+===========+===========+"  << endl;
    res_out << "|" << "||u||_1         " << "|" << setw(11) << n1_n2(1,25) << "|"<< setw(11) << n1_n2(1,50)<< "|" << setw(11) << n1_n2(1,100)<< "|" ;
    res_out << setw(11) << n1_n2(1,200)<<  "|" << setw(11) << n1_n2(1,500)<<"|" << setw(11) << n1_n2(1,1000)<< "|" << setw(11) << n1_n2(1,2000)<<"|" << setw(11) << n1_n2(1,5000)<<"|"<< endl;
    res_out << "+================+===========+===========+===========+===========+===========+===========+===========+===========+"  << endl;
    res_out << "|" << "||u||_2         " << "|" << setw(11) << n1_n2(2,25) << "|"<< setw(11) << n1_n2(2,50)<< "|" << setw(11) << n1_n2(2,100)<< "|" ;
    res_out << setw(11) << n1_n2(2,200)<<  "|" << setw(11) << n1_n2(2,500)<<"|" << setw(11) << n1_n2(2,1000)<< "|" << setw(11) << n1_n2(2,2000)<<"|" << setw(11) << n1_n2(2,5000)<<"|"<< endl;
    res_out << "+================+===========+===========+===========+===========+===========+===========+===========+===========+"  << endl;
    res_out.close();


}


int main() {
    double u[N_ + 1];
    ofstream res_out;
    res_out.open("result_u.txt");
    asymmetric_difference_derivative_method(u, N_ + 1);
    for(int i = 0; i < N_;i++){
        cout << u[i] << endl;
    }

    for (int i = 0; i < N_; i++) {
        res_out << x_i(N_,i)<< " ";
    }
    res_out << endl;
    for (int i = 0; i < N_; i++) {
        res_out << u[i]<< " ";
    }
    res_out.close();

    res_out.open("result_v.txt");
    dummy_node_method(u, N_);
    for(int i = 0; i < N_;i++){
        cout << u[i] << endl;
    }

    for (int i = 0; i < N_; i++) {
        res_out << x_i(N_,i)<< " ";
    }
    res_out << endl;
    for (int i = 0; i < N_; i++) {
        res_out << u[i]<< " ";
    }
    res_out.close();
    table_of_n1_n2();
    return 0;
    //zeidel_method(res_out);
}

// TIP See CLion help at <a
// href="https://www.jetbrains.com/help/clion/">jetbrains.com/help/clion/</a>.
//  Also, you can try interactive lessons for CLion by selecting
//  'Help | Learn IDE Features' from the main menu.