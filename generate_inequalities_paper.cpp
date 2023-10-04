#include <string_view>

#include "pqc.h"

using namespace std::literals;

void generate_inequations(bool filtered, long long int sd, int numIneq_each){
    ////Alice generate and fix public and private keys
    
    generator.seed(sd);

    string fileprefix;
    
 

    if (filtered) {
        fileprefix = "inequalities/filtered/filter_Kyber"+ to_string(version) + "Simplified_seed" + to_string(sd);      
    } else {
        fileprefix = "inequalities/unfiltered/unfilter_Kyber"+ to_string(version) + "Simplified_seed" + to_string(sd);
    }

    ofstream myout_pkA(fileprefix + "_pkA.txt"); //file for public key A
    ofstream myout_pkt(fileprefix + "_pkt.txt"); //file for public key t
    ofstream myout_Coeff(fileprefix + ".txt"); //file for coefficient matrix
    ofstream myout_e_s(fileprefix + "_es.txt"); // file for solution
    ofstream myout_lvu(fileprefix +"_lvu.txt"); //file for lower bound, true value, upper bound
    
    
    
    

    //generate A from R_q   
    vector<vector<vector<int>>> A(k, vector<vector<int>> (k, vector<int> (n)));
    //store lwe A
    vector<vector<int>> lwe_A(n*k, vector<int>(n*k));
    vector<vector<int>> temp(n, vector<int>(n));
    for (int i = 0; i < k; i++){
        for (int j = 0; j < k ; j++){
            A[i][j] = sample_poly_U();
            temp = coeff2matrix(A[i][j]);
            for (int ii = 0; ii < n; ii++){
                for (int jj = 0; jj < n; jj++){
                    lwe_A[i * n + ii][j * n + jj] = temp[ii][jj];
                }                
            }
        }
    }
    for (int i = 0; i < k*n; i++){
        for (int j = 0; j < k*n; j++){
            myout_pkA << lwe_A[i][j] << " ";
        }
        myout_pkA << endl;
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

    for (int i = 0 ; i < k ; i++){
        for (int j = 0; j < n; j++){
            myout_e_s << e[i][j] <<endl;
        }
        
    }

    for (int i = 0 ; i < k ; i++){
        for (int j = 0; j < n; j++){
            myout_e_s << s[i][j] <<endl;
        }
        
    }
    

    //calculate t =  As + e and store lwe_t
    vector<vector<int>> t(k, vector<int>(n));
    for (int i = 0 ; i < k; i++){
        t[i] = dotproduct(A[i], s) + e[i];
    }
    for (int i = 0; i < k; i++){
        for (int j = 0; j < n; j++){
            myout_pkt << t[i][j] << endl;
        }
    }

    
    
    
    int numIneq = 0;
    int flag;
    int index = 0; //v' index
    int collected = 0; 
   
    vector<vector<int>> addN(n, vector<int>(numIneq_each)); 
    
    
    int lowerbound, upperbound; 
    int delta_N; 
    int delta_Nadd1; 
    
    vector<int> m(n);    
    vector<int> m_prime(n);
    vector<int> a(n);
    vector<int> diff_v; //v'' - v 
    vector<vector<int>> diff_u(k, vector<int>(n)); //u'' - u 
    vector<vector<int> > Coeff_eq; //coefficient matrix
	int d, d_un, d_es;
    //int num_eq = 0;
    while ((index < n) & (collected < numIneq_each) ){
        // Bob sample a fixed m and all random values
        //sample m from {0,1}
        vector<int> m(n);
        for (int i = 0; i < n; i++){
            m[i] = distribution(generator); 
        }
        
        //sample e2 from B_eta2
        vector<int> e2(n);
        e2 = sample_poly_B(eta2);
        
        
        //sample r from B_eta1
        vector<vector<int>> r(k, vector<int>(n));
        for (int i = 0; i < k; i++){
            r[i] = sample_poly_B(eta1);
        }

        //calculate v = t'r + e2 +Decompress_q(m, 1)
        vector<int> v(n);
        v = dotproduct(t, r) + e2 + Decompress_q(m, 1);

        //calculate v' = Compress_q(v)
        vector<int> v_prime(n);
        v_prime = Compress_q(v, dv);

        //V''
        vector<int> v_primeprime(n);
        v_primeprime = Decompress_q(v_prime, dv);

        // diff_v
        diff_v = biased_mod(plain_minus(v_primeprime, v), q); 
		//print_vector(diff_v);
        //filture ciphertexts

        
        if (filtered) { //for filtered ciphertexts

            if (abs(e2[index] + diff_v[index]) > 10 ){ 
                continue; 
            }
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
        

        //calculate u' = Compress_q(u)
        vector<vector<int>> u_prime(k, vector<int>(n));
        for (int i = 0; i < k; i++){
            u_prime[i] = Compress_q(u[i], du);
        }
        

        
        //calculate u'' = Decompress_q(u') = Decompress_q(Compress_q(u))
        vector<vector<int>> u_primeprime(k, vector<int>(n)); //u''
        for (int i = 0; i < k; i++){
            u_primeprime[i] = Decompress_q(u_prime[i], du);
        }
        

        // diff_u
        for (int i = 0; i < k; i++){
            diff_u[i] = biased_mod(plain_minus(u_primeprime[i], u[i]), q) ;
        }        
        
	   
       	vector<vector<int>> neg_e1_plus_diffu(k, vector<int>(n));
       	for (int i = 0; i < k; i++){
        	neg_e1_plus_diffu[i] = plain_minus(-1 * e1[i], diff_u[i]); //plain value
            
       	}
		
		// noise d
        d = (plain_add(plain_dotproduct(r, e), plain_dotproduct(neg_e1_plus_diffu, s)), plain_add(e2, diff_v))[index];
	
		if ((d < q/4.0) && (d > -q/4.0)){

		} else {
			printf("d = %d \n", d);
            for (int i = 0; i < k; i++){
                print_vector(e1[i]);
            }

			//printf("d_un = %d \n", d_un);
			//printf("d_es = %d \n", d_es);
			printf("s * diff_u = %d\n", diff_v[index] - plain_dotproduct(diff_u, s)[index]);
			return;
		}
        
       

        vector<int> coeff_plain(k*n + k*n);

		 
        
        coeff_plain = coeffe_coeffs(r, neg_e1_plus_diffu, index);
        
        for (int i = 0; i < n*k *2; i++){
            myout_Coeff <<  coeff_plain[i] <<" "; //
        }
        myout_Coeff << endl;

        //Bob alter v' step by step
        vector<int> v_prime_error(n);
        v_prime_error = v_prime;
        v_prime_error[index]++;
        int v_primeprime_error;

        
        //add N
        v_prime_error[index] = mod(v_prime[index]+ 2, pow(2, dv));
        v_primeprime_error = Decompress_q(v_prime_error[index], dv);
        m_prime[index] = Compress_q(mod(v_primeprime_error - dotproduct(s, u_primeprime)[index], q), 1);
        
        while (m[index] == m_prime[index]) {           
            v_prime_error[index]++; 
            
            v_primeprime_error = Decompress_q(mod(v_prime_error[index], pow(2, dv)), dv);
            //calculate m'
            m_prime[index] = Compress_q(mod(v_primeprime_error - dotproduct(s, u_primeprime)[index], q), 1);
        };
        
        addN[index][collected] = v_prime_error[index] - v_prime[index] - 1;
        
      
        delta_N = Decompress_q(mod(v_prime[index] + addN[index][collected], pow(2, dv)), dv) - Decompress_q(v_prime[index], dv);
        delta_Nadd1 = Decompress_q(mod(v_prime[index] + addN[index][collected]+1, pow(2, dv)), dv) - Decompress_q(v_prime[index], dv);
       
        if (m[index] == 0){   // 0 -> 1        
 
            lowerbound = biased_mod(ceil(q/4.0) - delta_Nadd1, q);
            upperbound = biased_mod(floor( q/4.0) - delta_N, q);

        }else{// 1 -> 0
            //cout << addN[index][collected] << endl; 
            lowerbound = biased_mod(ceil(3*q/4.0) - round(q/2.0) -delta_Nadd1, q);
            upperbound =  biased_mod(floor(3*q/4.0) - round(q/2.0) - delta_N, q);   
                          
        }      
            
        
           
    

        //verify inequalities
        int value_plain;
        value_plain = 0;
        for (int i = 0; i < k; i++){
            for (int j = 0; j < n; j++){
                value_plain += coeff_plain[i*n + j]*e[i][j];
                value_plain += coeff_plain[k*n + i*n + j]* s[i][j];
            }
            
        }
        if ((value_plain <= upperbound - e2[index] - diff_v[index]) & (value_plain >= lowerbound - e2[index] - diff_v[index])){ 
            myout_lvu << lowerbound - e2[index] - diff_v[index] <<" "<< value_plain <<" "<< upperbound - e2[index] - diff_v[index]<<endl;//
            
        }
        else{
            cout<<"error value_plain = "<< value_plain <<", bound = [" << lowerbound - e2[index] - diff_v[index] <<", " <<upperbound - e2[index] - diff_v[index]<<"], "<<endl;  
            cout<<"d = "<< d << " [" << lowerbound << ", "<< upperbound << "]"<<endl;
           
        }
        index++;
        if (index == n){                
            collected++;
            cout<<"Collected "<<collected<<" inequalities for all indexes"<<endl;
            index = 0;
        }

            
        
    }

    
    myout_pkA.close();
    myout_pkt.close();
    myout_Coeff.close();
    myout_e_s.close();
    myout_lvu.close();
    
}

int main(int argc,char * argv[]){

    //generator.seed(10*time(NULL));
    bool filtered;
    int numIneq_each;
    
    std::string s(argv[1]);
    if ( s == "filtered") {
        filtered = true;
        printf("Generating filterd inequalities...");
        numIneq_each = 50; 
    }else if (s == "unfiltered"){
        filtered = false;
        printf("Generating unfilterd inequalities...");
        numIneq_each = 50; 
    }else {
        printf("No indication of filtering. Give 'filtered' or 'unfiltered' \n");
        return 0;
    }


    long long int seed;
    seed = atoll(argv[2]);
    printf("Using seed %lld \n", seed);
    
    
    
    //test_diffIneq_validity();
   
    generate_inequations(filtered, seed, numIneq_each);
   
    return 0;
}
