#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <queue>
#include <limits.h>
#include "umfpack.h"


using namespace std;

#define number_of_pins 510
#define number_of_nets 510
#define max_connection 100
#define max_fanout 20


struct Fixed_pin{
    int pin_number;
    int x_loc;
    int y_loc;
};

struct Pin{
    int pin_number;
    int connection_number;
    //net connection
    int* connection; 
    bool fixed;
};

struct Net{
    int connected_block_len;
    int* connected_block;
    double* weight;
};


//Global Variables
Fixed_pin* fp;
Pin* p;
Net* n;
int num_of_blocks = 0;
int num_of_fixed = 0;


void init(){
    p = new Pin[number_of_pins];
    fp = new Fixed_pin[number_of_pins];
    n = new Net[number_of_nets];
    for(int i = 0; i < number_of_pins; i++){
        n[i].connected_block = new int[max_fanout];
        n[i].weight = new double[max_fanout];
        n[i].connected_block_len = 0;
    }
}

//create clique model
void make_clique(){
    for(int i = 0; i < number_of_nets; i++){
        if(n[i].connected_block_len == 0){
            continue;
        }
        else{
            for(int j = 0; j < n[i].connected_block_len; j++){
                n[i].weight[j] = (2.0/ n[i].connected_block_len);
            }
        }   
    }
}

double find_weight(int net_number, int target_block_number){
    for(int i = 0; i < n[net_number].connected_block_len; i++){
        if(n[net_number].connected_block[i] == target_block_number)
            return n[net_number].weight[i];
    }
    return -1;
}

int find_pin(int pin_number){
    for(int i = 0; i < num_of_blocks; i++)
        if(p[i].pin_number == pin_number)
            return i;
    return -1;
}


int find_fixed_pin(int pin_number){
    for(int i = 0; i < num_of_blocks; i++)
        if(fp[i].pin_number == pin_number)
            return i;
    return -1;
}


//Calculate Weight Matrix
double** cal_weight(){

    double** pre_weight_mat = new double*[num_of_blocks];
    double** weight_mat = new double*[num_of_blocks - num_of_fixed];
    
    for(int i = 0 ;  i < num_of_blocks - num_of_fixed; i++)
        weight_mat[i] = new double[num_of_blocks - num_of_fixed];
    
    for(int i = 0 ; i < num_of_blocks; i++){
        pre_weight_mat[i] = new double[num_of_blocks];
    }
    //init to zero
    for(int i = 0; i < num_of_blocks; i++){
        for(int j = 0; j < num_of_blocks; j++){
            pre_weight_mat[i][j] = 0;
        }
    } 
    for(int i = 0; i < num_of_blocks; i++){
        for(int j = 0; j < num_of_blocks; j++){
            if(i == j){//diagonal elements
                pre_weight_mat[i][j] = 0;
            }
            else{
                for(int k = 0; k < p[i].connection_number; k++){
                    double weight_temp = find_weight(p[i].connection[k],p[j].pin_number);
                    if(weight_temp != -1)
                        pre_weight_mat[i][j] += weight_temp;
                }
            }
        }
    }
    /************************************************************Debug*****************************************************************/
    // for(int i = 0; i < num_of_blocks; i++){
    //     for(int j = 0; j < num_of_blocks; j++){
    //         cout << pre_weight_mat[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    
    int ii = 0;
    int jj = 0;
    for(int i = 0; i < num_of_blocks; i++){
        for(int j = 0; j < num_of_blocks; j++){
            if(p[i].fixed || p[j].fixed)
                continue;
            else{
                weight_mat[ii][jj] = pre_weight_mat[i][j];
                if(jj < num_of_blocks - num_of_fixed - 1){
                    jj++;
                }
                else{
                    jj = 0;
                    ii++;
                }
            }
        }
    }
    int p_index = 0;
    for(int i = 0; i < num_of_blocks - num_of_fixed; i++){
        for(int j = 0; j < num_of_blocks- num_of_fixed; j++){
            if(i != j){
                if(weight_mat[i][j] != 0)
                    weight_mat[i][j] = -1 * weight_mat[i][j];
            }
            else{
                while (p[p_index].fixed)
                {
                    p_index++;
                }
                
                double w = 0;
                for(int k = 0; k < num_of_blocks; k++){
                    w += pre_weight_mat[p[p_index].pin_number-1][k];
                }
                p_index++;
                weight_mat[i][j] = w;
            }
        }
    }
    /************************************************************Debug*****************************************************************/
    // for(int i = 0; i < num_of_blocks-num_of_fixed; i++){
    //     for(int j = 0; j < num_of_blocks-num_of_fixed; j++){
    //         cout << weight_mat[i][j] << "  ";
    //     }
    //     cout << endl;
    // }
    return weight_mat;
}

double* cal_fixed(int dim){//dim 0 : x_vector dim 1: y_vecotr
    double* dim_vec = new double[num_of_blocks - num_of_fixed];
    int index_dim_vec = 0;
    for(int i = 0; i < num_of_blocks - num_of_fixed; i++)
        dim_vec[i] = 0;
    
    for(int i = 0; i < num_of_blocks; i++){
        if(p[i].fixed)
            continue;
        for(int j = 0; j < p[i].connection_number; j++){
            int net_id = p[i].connection[j];
            if(n[net_id].connected_block_len == 0)
                continue;
            else{
                for(int k = 0; k < n[net_id].connected_block_len; k++){
                    int p_id = find_pin(n[net_id].connected_block[k]);
                    if(p[p_id].fixed){//p[i] is connected to a fixed block
                        int fp_in = find_fixed_pin(p[p_id].pin_number);
                        double pre_weight = (dim < 1) ? (n[net_id].weight[k] * fp[fp_in].x_loc) : (n[net_id].weight[k] * fp[fp_in].y_loc);
                        dim_vec[index_dim_vec] = dim_vec[index_dim_vec] + pre_weight;
                    }
                }
            }
        }
        index_dim_vec++;
    }
    
    /************************************************************Debug*****************************************************************/
    // for(int i = 0; i < num_of_blocks - num_of_fixed; i++){
    //     cout << " " << dim_vec[i] << " ";
    // }
    // cout << endl;
    return dim_vec;
}

double* solve(int dim){
    double** weight_mat = cal_weight();
    int non_zero_elem = 0;
    for(int i = 0 ; i < num_of_blocks - num_of_fixed; i++){
        for(int j = 0; j < num_of_blocks - num_of_fixed; j++){
            if(weight_mat[i][j] != 0)
                non_zero_elem++;
        }
    }
    //preprocess matrix for solver
    int* Ap = new int[num_of_blocks - num_of_fixed + 1];
    Ap[0] = 0;
    int Ap_index = 1;
    int* Ai = new int[non_zero_elem];
    int Ai_Ax_index = 0;
    double* Ax = new double[non_zero_elem];
    for(int j = 0; j < num_of_blocks - num_of_fixed; j++){
        Ap[Ap_index] = Ap[Ap_index-1];
        for(int i = 0; i < num_of_blocks - num_of_fixed; i++){
            if(weight_mat[i][j] != 0){
                Ai[Ai_Ax_index] = i;
                Ax[Ai_Ax_index] = weight_mat[i][j];
                Ai_Ax_index++;
                Ap[Ap_index]++;
            }
        }
        Ap_index++;
    }
    int n = num_of_blocks - num_of_fixed;
    double* in_vector = cal_fixed(dim);
    double* x = new double [num_of_blocks - num_of_fixed];
    //solve for x dimension
    double *null = (double *) NULL ;
    void *Symbolic, *Numeric ;
    (void) umfpack_di_symbolic (n, n, Ap, Ai, Ax, &Symbolic, null, null) ;
    (void) umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, null, null) ;
    umfpack_di_free_symbolic (&Symbolic) ;
    (void) umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, x, in_vector, Numeric, null, null) ;
    umfpack_di_free_numeric (&Numeric) ;
    /************************************************************Debug*****************************************************************/
    // for(int i = 0 ; i < n ; i++){
    //     cout << x[i] << " ";
    // }
    // cout << endl;
    return x;
}


int main(){
    string fileName;
    cout << "Enter the name of test file (cct1.txt)" << endl;
    cin >> fileName; 
    ifstream testFile(fileName);
    
    //init global variables
    init();

    //load input into data structures
    int input;
    int counter = 0;  
    bool end_of_netlist = false;
    while(testFile >> input){
        if(input == -1 && end_of_netlist)
            break;
        else if (input == -1 && !end_of_netlist){
            end_of_netlist = true;
            counter = 0;
            num_of_blocks++;
            continue;
        }
        else{
            end_of_netlist = false;
            if(counter == 0){//new row
                p[num_of_blocks].pin_number= input;
                p[num_of_blocks].connection_number = 0;
                p[num_of_blocks].connection = new int[max_connection];
                p[num_of_blocks].fixed = false;
                counter++;
            }
            else{
                p[num_of_blocks].connection[counter-1]= input;
                p[num_of_blocks].connection_number++;
                n[input].connected_block[n[input].connected_block_len] = p[num_of_blocks].pin_number;
                n[input].connected_block_len++;
                counter++;
            }
        }
    }
    counter = 0;
    int index;
    while(testFile >> input){
        if(input == -1)
            break;
        else{
            switch (counter)
            {
            case 0:
                fp[num_of_fixed].pin_number = input; 
                index = find_pin(input); 
                p[index].fixed = true;
                break;
            case 1:
                fp[num_of_fixed].x_loc = input;
                break;
            case 2:
                fp[num_of_fixed].y_loc = input;
                break;
            default:
                break;
            }
            if(counter == 2){
                counter = 0;
                num_of_fixed++;
            }
            else{
                counter++;
            }
        }
    }
    /************************************************************Debug*****************************************************************/
    //print input netlist
    // for(int i = 0; i < num_of_blocks ; i++){
    //     cout << "Pin Number: " << p[i].pin_number << endl;
    //     cout << "Connection:" << endl;
    //     for(int j = 0; j < p[i].connection_number; j++){           
    //         cout << p[i].connection[j] << " ";
    //     }
    //     cout << endl;
    // }
    //print netlists 
    // for(int i = 0; i < number_of_nets; i++){
    //     if(n[i].connected_block_len != 0){
    //         cout << i << " : " << endl;
    //         for(int j = 0; j < n[i].connected_block_len; j++){
    //             cout << n[i].connected_block[j] << " " << n[i].weight[j] << " ";
    //         }
    //         cout << endl;
    //     }
    // }
    //print fixed items
    // for(int i = 0; i < num_of_fixed; i++){
    //     cout << fp[i].pin_number << " " << fp[i].x_loc << " " << fp[i].y_loc << endl;
    // }
    
    make_clique();
    double* x_vec = solve(0);
    double* y_vec = solve(1);
    /************************************************************Debug*****************************************************************/
    for(int i = 0 ; i < num_of_blocks - num_of_fixed; i++){
        cout << x_vec[i] << " " << y_vec[i] << endl; 
    }
    
    return 0;
}