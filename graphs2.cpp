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

double x_i(int i){
    return a_ + ((b_ - a_)/(N_))*i;
}

double p(int i){
    return x_i(i) + 1.0;
}

double q(int i){
    return -1.0;
}

double g(int i){
    double x = x_i(i);
    return (x*x + 2.0*x + 2.0)/(x + 1.0);
}

struct three_diagonals_matrix {
private:
    double a[N_+1];
    double b[N_+1];
    double c[N_+1];
    double f[N_+1];
public:
    three_diagonals_matrix(){
       for(int i = 0; i < N_ + 2; i++){
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
        if (i > N_ + 1){
            cout << " i > N + 1" << endl;
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

void solve_three_diagonals_matrix(three_diagonals_matrix* M, ofstream & file)
{
    double temp;
    double alpha[N_ + 2];//первое поле использовать не буду, ленивый ыыыыы, лень двигать индексы, говнокод, да, но я очень спешу
    double beta[N_ + 2];
    double u[N_ + 1];
    alpha[1] = M->bi(0)/M->ci(0);
    beta[1] = M->fi(0)/M->ci(0);

    for (int i = 1; i < N_ - 1; i++)
    {
        alpha[i + 1] = M->bi(i)/(M->ci(i) - alpha[i]*M->ai(i));
        beta[i + 1] =(M->fi(i) + beta[i]*M->ai(i))/(M->ci(i) - alpha[i]*M->ai(i));
    }
    beta[N_] =(M->fi(N_-1) + beta[N_-1]*M->ai(N_-1))/(M->ci(N_-1) - alpha[N_-1]*M->ai(N_-1));

    u[N_-1] = beta[N_];
    temp = beta[N_];

    for (int i = N_ - 1; i > 0; i--)
    {
        //temp = *(u + i);
        u[i - 1] = alpha[i]*u[i] + beta[i];
        temp = u[i - 1];
    }

    for(int i = 0; i < N_;i++){
        cout << u[i] << endl;
    }

    for (int i = 0; i < N_; i++) {
        file << x_i(i)<< " ";
    }
    file << endl;
    for (int i = 0; i < N_; i++) {
        file << u[i]<< " ";
    }


}

void dummy_node_method() {
    double a,b,c,f;
    double h = (b_ - a_)/N_;
    three_diagonals_matrix M;

    //заполнение трехдиагональной матрицы со страницы 7
    a = 0;
    b = mu_a/(h*h);
    c = mu_a/(h*h) - lambda_a/h + (lambda_a*p(0)) - mu_a*q(0)/2;
    f = -(1/h + p(0)/2)*psi_a - mu_a*g(0)/2;
    M.set_values(a, b, c,f,0);
    for(int i = 1; i < N_ - 1; i++){
        a = 1.0/(h*h) - p(i)/(2 * h);
        b = 1.0/(h*h) + p(i)/(2 * h);
        c = 2.0/(h*h) - q(i);
        f = -g(i);
        M.set_values(a, b, c,f,i);
    }
    // N -1 строка
    a = 1.0/(h*h) - p(N_-1)/(2 * h);
    b = 0;
    c = 2.0/(h*h) - q(N_ - 1);
    f = (1.0/(h*h) + p(N_ - 1)/(2 * h))*(psi_b/lambda_b) - g(N_ - 1);
    M.set_values(a, b, c,f,N_ - 1);


    ofstream res_out;
    res_out.open("result_v.txt");
    solve_three_diagonals_matrix(&M, res_out);
    res_out.close();
}

void asymmetric_difference_derivative_method(){
    double a,b,c,f, temp;
    double h = (b_ - a_)/N_;
    three_diagonals_matrix M;

    //заполнение трехдиагональной матрицы со страницы 7
    a = 0;
    b = mu_a/(h*h) + (mu_a*p(1))/h + mu_a*q(1)/2;
    c = mu_a/(h*h) + ((mu_a*p(1)) - lambda_a)/h - lambda_a*p(1)/2;
    f = -(1.0/h + p(1)/2.0)*psi_a - mu_a*g(1)/2.0;
    M.set_values(a, b, c,f,0);
    for(int i = 1; i < N_ - 1; i++){
        a = 1.0/(h*h) - p(i)/(2.0 * h);
        b = 1.0 / (h * h) + p(i) / (2.0 * h);
        c = 2.0/(h*h) - q(i);
        f = -g(i);
        M.set_values(a, b, c,f,i);
    }
    // N -1 строка
    a = 1.0/(h*h) - p(N_-1)/(2.0 * h);
    b = 0;
    c = 2.0/(h*h) - q(N_ - 1);
    f = (1.0/(h*h) + p(N_ - 1)/(2.0 * h))*(psi_b/lambda_b) - g(N_ - 1);
    M.set_values(a, b, c,f,N_ - 1);

    ofstream res_out;
    res_out.open("result_u.txt");
    solve_three_diagonals_matrix(&M, res_out);
    res_out.close();

}

int main() {
    ofstream res_out;
    res_out.open("result.txt");
    res_out << "eps   = " << eps << endl;
    res_out << "+================================================================================================================+" << endl;
    res_out << "|  n1/n2         | 25/50 | 50/100 | 100/200 | 200/400 | 500/10^3 | 10^3/2 / 10^3 | 2*10^3 / 4*10^3 | 5*10^3/10^4 |"<< endl;
    res_out << "+================+=======+========+=========+=========+==========+===============+=================+=============+" << endl;
    //zeidel_method(res_out);
    asymmetric_difference_derivative_method();
    dummy_node_method();
    res_out.close();
    return 0;
}

// TIP See CLion help at <a
// href="https://www.jetbrains.com/help/clion/">jetbrains.com/help/clion/</a>.
//  Also, you can try interactive lessons for CLion by selecting
//  'Help | Learn IDE Features' from the main menu.