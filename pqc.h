#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <algorithm>
#include <bitset>
#include <cmath>
#include "assert.h"


using namespace std;


#define q  3329
#define n 256

//kyber512
#define version 512
#define k 2
#define eta1 3
#define eta2 2
#define du 10
#define dv 4
/*
//kyber768
#define version 768
#define k 3
#define eta1 2
#define eta2 2
#define du 10
#define dv 4
*/
/*
//kyber1024
#define version 1024
#define k 4
#define eta1 2
#define eta2 2
#define du 11
#define dv 5
*/

//#define numIneq_each 25 

mt19937 generator;
//default_random_engine generator;
bernoulli_distribution distribution(0.5); // this distributil will call random engine defined in previous line
uniform_int_distribution<int> uniform(0, q-1);

string curtime(){
    char name[256] = {0};
    time_t timep;
    time(&timep);
    strftime(name, sizeof(name),"_%Y%m%d%H%M%S",localtime(&timep) );
    return name;
}


int B_eta(const int eta){
    int res=0; 
    for (int i = 0; i < eta; i++){
        res += distribution(generator) - distribution(generator);
    }
    return res;
}


vector<int> sample_poly_B(int eta){
    vector<int> poly(n);
    for (int i = 0; i< n; i++){
        poly[i] = B_eta(eta);
    }
    return poly;
}


vector<int> sample_poly_U(){
    vector<int> poly(n);
    for (int i = 0; i< n; i++){
        poly[i] = uniform(generator);
    }
    return poly;
}

int mod(int a, int b){
    return ((a % b) + b) % b;
}
vector<int> mod(vector<int> a, int b){
    vector<int> c(a.size());
    for (int i = 0; i < a.size(); i++){
        c[i] = mod(a[i], b);
    }
    return c;
}
int Compress_q(int x, float d){ // x in Z_q
    assert(x < q);       
    assert(x >= 0); 
    return mod((int) round(pow(2,d)/q * x ), (int) pow(2, d));
}

int Decompress_q(int x, float d){
    assert(x < pow(2, d));       
    assert(x >= 0); 

    return (int) round(q/pow(2, d) * x);
}

vector<int> Compress_q(vector<int> a, float d){
    vector<int> res(a.size());
    for (int i = 0; i < a.size(); i++){
        res[i] = Compress_q(a[i], d);
    }
    return res;
}

vector<int> Decompress_q(vector<int> a, float d){
    vector<int> res(a.size());
    for (int i = 0; i < a.size(); i++){
        res[i] = Decompress_q(a[i], d);
    }
    return res;
}

int centered_mod(int x, int a){
    int m;
    m = mod(x, a);
    return (m > a/2 )? (m-a) : m;
}

int biased_mod(int x, int a){
    //range in [-1/4 * a, 3/4 * a)
    int m; 
    m = mod(x, a);
    return (m >= 0.75 * a) ? (m - a): m;
}
vector<int> biased_mod(vector<int> x, int a){
    vector<int> res(x.size());
    for (int i = 0; i < x.size(); i++){
        res[i] = biased_mod(x[i], a);
    }
    return res;
}
vector<vector<int>> biased_mod(vector<vector<int>> x, int a){
    vector<vector<int>> res(x.size());
    for (int i = 0; i < x.size(); i++){
        res[i] = biased_mod(x[i], a);
    }
    return res;
}
void print_matrix(int * M, int r, int c){
    for (int i = 0; i<r; i++){
        for (int j = 0; j < c; j++){
            cout<< * (M + i*c +j) << " ";
        }
        cout << endl;
    }
}
void print_vector(vector<int> a){
    for (int i = 0; i < a.size(); i++){
        cout<<a[i]<<" ";
    }
    cout<<endl;
}

vector<vector<int>> coeff2matrix(vector<int> a){
    vector<vector<int>> coeff_matrix(n, vector<int>(n));
    
    for (int i = 0 ; i < n; i++){
        coeff_matrix[i][0] = a[i];
    }
    int j=1;
    while (j < n){
        coeff_matrix[0][j] = mod(-1 * coeff_matrix[n-1][j-1], q);
        for (int i = 1; i < n; i++){
            coeff_matrix[i][j] = coeff_matrix[i-1][j-1];
        }       
        j++;
    }

    return coeff_matrix;
}
vector<int> operator*(vector<int> a, vector<int> b){
    
    int coeff_matrix[n][n];
    
    for (int i = 0 ; i < n; i++){
        coeff_matrix[i][0] = a[i];
    }
    int j=1;
    while (j < n){
        coeff_matrix[0][j] = mod(-1 * coeff_matrix[n-1][j-1], q);
        for (int i = 1; i < n; i++){
            coeff_matrix[i][j] = coeff_matrix[i-1][j-1];
        }       
        j++;
    }
    vector<int> c(n);
    for (int i= 0; i < n; i++){
        c[i] = 0;
        for (j = 0; j < n; j++){
            c[i] = mod(c[i] + coeff_matrix[i][j] * b[j], q);
        }
    }
    return c;
}

vector<int> plain_mult(vector<int> a, vector<int> b){
    
    int coeff_matrix[n][n];
    
    for (int i = 0 ; i < n; i++){
        coeff_matrix[i][0] = a[i];
    }
    int j=1;
    while (j < n){
        coeff_matrix[0][j] = -1 * coeff_matrix[n-1][j-1];
        for (int i = 1; i < n; i++){
            coeff_matrix[i][j] = coeff_matrix[i-1][j-1];
        }       
        j++;
    }
    vector<int> c(n);
    for (int i= 0; i < n; i++){
        c[i] = 0;
        for (j = 0; j < n; j++){
            c[i] = c[i] + coeff_matrix[i][j] * b[j];
        }
    }
    return c;
}

vector<int> operator*(int coe, vector<int> b){
    vector<int> c(b.size());
    for (int i = 0; i< b.size(); i++){
        
        c[i] = coe * b[i];
       
    }
    return c;
}

vector<int> operator+(vector<int> a, vector<int> b){
    vector<int> c(a.size());
    for (int i = 0; i< a.size(); i++){
        
        c[i] = mod(a[i] + b[i], q);
       
    }
    return c;
}
vector<int> plain_add(vector<int> a, vector<int> b){
    vector<int> c(a.size());
    for (int i = 0; i< a.size(); i++){
        
        c[i] = a[i] + b[i];
       
    }
    return c;
}
vector<int> operator-(vector<int> a, vector<int> b){
    vector<int> c(a.size());
    for (int i = 0; i< a.size(); i++){
        c[i] = mod(a[i] - b[i], q);
    }
    return c;
}

vector<vector<int>> operator-(vector<vector<int>> a, vector<vector<int>> b){
    vector<vector<int>> c(a.size());
    for (int i = 0; i< a.size(); i++){
        c[i] = mod(a[i] - b[i], q);
    }
    return c;
}
vector<int> plain_minus(vector<int> a, vector<int> b){
    vector<int> c(a.size());
    for (int i = 0; i< a.size(); i++){
        
        c[i] = a[i] - b[i];
       
    }
    return c;
}
vector<int> dotproduct(vector<vector<int>> a, vector<vector<int>> b){
    vector<int> c(n);
    c = a[0] * b[0];
   // vector<int> temp(n);
    for (int i = 1; i < k; i++){
        
       // temp = a[i] * b[i];
        c = c +  a[i] * b[i];
    }
    return c;
}
vector<int> plain_dotproduct(vector<vector<int>> a, vector<vector<int>> b){
    vector<int> c(n);
    c = plain_mult(a[0] , b[0]);
   // vector<int> temp(n);
    for (int i = 1; i < k; i++){
        
       // temp = a[i] * b[i];
        c = plain_add(c,  plain_mult(a[i], b[i]));
    }
    return c;
}
void test_plain_offset(){
    //sample e from B_eta1;
    vector<vector<int>> e(k, vector<int>(n));
    for (int i = 0; i < k; i++){
        e[i] = sample_poly_B(eta1);
    }
    
    //sample r from B_eta1
    vector<vector<int>> r(k, vector<int>(n));
    for (int i = 0; i < k; i++){
        r[i] = sample_poly_B(eta1);
    }
    //sample e2 from B_eta2
    vector<int> e2 = sample_poly_B(eta2);
    //sample  s from B_eta1
    vector<vector<int>> s(k, vector<int>(n));
    for (int i = 0; i < k; i++){
        s[i] = sample_poly_B(eta1);
    }
    //sample e1 from B_eta2
    vector<vector<int>> e1(k, vector<int>(n));
    for (int i = 0; i < k; i++){
        e1[i] = sample_poly_B(eta2);
    }
    
    //sample m from {0,1}
    vector<int> m(n);
    for (int i = 0; i < n; i++){
        m[i] = distribution(generator);
    }

    vector<int> w(n);
    
    w = dotproduct(e, r) + e2 - dotproduct(s,e1);
    
    vector<int> m_prime(n);

    m_prime = Compress_q(w + Decompress_q(m, 1), 1);
    
    for (int i = 0; i < n; i++){
              
        if ((w[i] >= q/pow(2, 2)) & (w[i] <=  3*q/pow(2, 2))){
            cout << "test w " << w[i] << " "<< Decompress_q(m, 1)[i] <<endl;
        }

        if(m_prime[i] != m[i]){
            cout << i <<endl;
        }
    }

    
}
int distance_DecompCompr_q(int x, float d){
    //calculate res = DecompCompr_q(x) - x as an integer. maybe negetive or positive
    //Note that x and DecompCompr_q(x) are in Z_q according to our calculation rule
    //z in (-q, q) and z satisfy that |z centered_mod q| <= q/(2^{d+1})
    return Decompress_q(Compress_q(x, d), d) - x;
}

vector<int> distance_DecompCompr_q(vector<int> x, float d){
    vector<int> res(x.size());
    for (int i = 0; i < x.size(); i++){
        res[i] = distance_DecompCompr_q(x[i], d);
    }
    return res;
}
int quotient(int a, int b){
    // a = xb + r, 0<= r < b.
    // return x
    return (a - mod(a, b))/b;
}
int biased_quotient(int a, int b){
    // a = xb + r, -b/4<= r < 3b/4.
    // return x
    return (a - biased_mod(a, b))/b;
}
vector<int> quotient(vector<int> a, int b){
    vector<int> res(a.size());
    for (int i = 0; i < a.size(); i++){
        res[i] = quotient(a[i], b);
    }
    return res;
}
void test_DecCompr_offset(){
    //generate A from R_q
    
    vector<vector<vector<int>>> A(k, vector<vector<int>> (k, vector<int> (n)));
   
    for (int i = 0; i < k; i++){
        for (int j = 0; j < k ; j++){
            A[i][j] = sample_poly_U();
        }
    }
    
    //sample  s from B_eta1
    vector<vector<int>> s(k, vector<int>(n));
    for (int i = 0; i < k; i++){
        s[i] = sample_poly_B(eta1);
    }
    //sample e from B_eta1;
    vector<vector<int>> e(k, vector<int>(n));
    for (int i = 0; i < k; i++){
        e[i] = sample_poly_B(eta1);
    }
    //calculate t =  As + e
    vector<vector<int>> t(k, vector<int>(n));
    for (int i = 0 ; i < k; i++){
        t[i] = dotproduct(A[i], s) + e[i];
    }

    //sample r from B_eta1
     vector<vector<int>> r(k, vector<int>(n));
    for (int i = 0; i < k; i++){
        r[i] = sample_poly_B(eta1);
    }
    //sample e1 from B_eta2
    vector<vector<int>> e1(k, vector<int>(n));
    for (int i = 0; i < k; i++){
        e1[i] = sample_poly_B(eta2);
    }

    //calculate u = A'r + e1. Note: Tansform of A
    vector<vector<int>> u(k, vector<int>(n));
    for (int i = 0 ; i < k; i++){
        vector<vector<int>> A_current_col(k, vector<int> (n));
        for (int j = 0; j < k; j++){
            A_current_col[j] = A[j][i];
        }
        u[i] = dotproduct(A_current_col, r) + e1[i];
        
    }

    //sample e2 from B_eta2
    vector<int> e2(n);
    e2 = sample_poly_B(eta2);

    //sample m from {0,1}
    vector<int> m(n);
    for (int i = 0; i < n; i++){
        m[i] = distribution(generator);
    }
    //calculate v = t'r + e2 +Decompress_q(m, 1)
    vector<int> v(n);
    v = dotproduct(t, r) + e2 + Decompress_q(m, 1);

    //calculate u' = Compress_q(u)
    vector<vector<int>> u_prime(k, vector<int>(n));
    for (int i = 0; i < k; i++){
        u_prime[i] = Compress_q(u[i], du);
    }
    //calculate v' = Compress_q(v)
    vector<int> v_prime(n);
    v_prime = Compress_q(v, dv);



    //calculate u'' = Decompress_q(u') = Decompress_q(Compress_q(u))
    vector<vector<int>> u_primeprime(k, vector<int>(n));
    for (int i = 0; i < k; i++){
        u_primeprime[i] = Decompress_q(u_prime[i], du);
    }
    //calculate v'' =  Decompress_q(v') = Decompress_q(Compress_q(v))
    vector<int> v_primeprime(n);
    v_primeprime = Decompress_q(v_prime, dv);

    
    //calculate m'
    vector<int> m_prime(n);
    m_prime = Compress_q(v_primeprime - dotproduct(s, u_primeprime), 1);

    //test if m'=m
    for (int i = 0; i < n; i++){
        if (m_prime[i] != m[i]){
            cout<< i <<endl;
        }
    }
}



vector<vector<int>> gen_coeffMatrix(vector<int> a){
    assert (a.size() == n);

    vector<vector<int>> coeff_matrix(n, vector<int> (n));
    
    for (int i = 0 ; i < n; i++){
        coeff_matrix[i][0] = a[i];
    }
    int j=1;
    while (j < n){
        coeff_matrix[0][j] = -1 * coeff_matrix[n-1][j-1], q;//Note this is not mod value
        for (int i = 1; i < n; i++){
            coeff_matrix[i][j] = coeff_matrix[i-1][j-1];
        }       
        j++;
    }

    return coeff_matrix;
}
vector<int> coeffe_coeffs(vector<vector<int>> r, vector<vector<int>> d, int index){
    assert(r.size() == k);
    assert(d.size() == k);
    assert(r[0].size() == n);
    assert(d[0].size() == n);

    vector<int> res(k*n + k*n);
    vector<int> temp(n);
    for (int i = 0; i < k; i++){
       
        temp = gen_coeffMatrix(r[i])[index];
        for (int j = 0; j < n; j++){
            res[i*n +j] = temp[j];
        }     
    }
    for (int i = 0; i < k; i++){
        temp = gen_coeffMatrix(d[i])[index];
        for (int j = 0; j < n; j++){
            res[k*n + i*n +j] = temp[j];
        }     
    }
    return res;
}

void Decompress_q_truthtable(int d){
    for (int i = 0; i < pow(2, d); i++){
        if (d==4)
            cout << bitset<4>(i) << " " << Decompress_q(i, d) << endl;
        else if (d==10)
            cout << bitset<10>(i) << " " << Decompress_q(i, d) << endl;
    }
}
