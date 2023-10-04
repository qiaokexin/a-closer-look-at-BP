Codes for generating inequalities on secret variables with intervel q/16 for Kyber512.
To compile the codes, run
'''
g++ generate_inequalities_paper.cpp -o gen 
'''
Then to generate inequalities for unfiltered ciphertexts with random seed 12345, run
'''
./gen unfiltered 12345
'''
To generate inequalities for filtered ciphertexts with random seed 12345, run
'''
./gen filtered 12345
'''
