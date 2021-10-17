#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <queue>
#include <limits.h>
#include <math.h>
#include <bits/stdc++.h>
#include "umfpack.h"

using namespace std;

#define number_of_pins 510
#define number_of_nets 510
#define max_connection 100
#define max_fanout 20

struct Bin{
    int x;
    int y;
    int block_count;
    int* blocks;
};

struct Fixed_pin{
    int pin_number;
    double x_loc;
    double y_loc;
};

struct Pin{
    int pin_number;
    int connection_number;
    //net connection
    int* connection; 
    bool fixed;
    double x_loc;
    double y_loc;
    double moved_x_loc;
    double moved_y_loc;
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
Bin* b;
Bin* overfilled_b;
int num_of_blocks = 0;
int num_of_fixed = 0;
int num_of_bins = 0;
int overfilled_bins_num =0;
int max_x = -1;
int max_y = -1;


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
    //delete pre_mat
    for(int i = 0; i < num_of_blocks; i++)
        delete[] pre_weight_mat[i];
    delete[] pre_weight_mat;
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
    delete[] Ap;
    delete[] Ax;
    delete[] Ai;
    return x;
}

int find_bin(int x, int y){
    for(int i = 0; i < num_of_bins; i++){
        if(b[i].x == x && b[i].y == y)
            return i;
    }
    return -1;
}

void create_bin(){
    b = new Bin [max_x*max_y];
    
    for(int i = 0 ;  i < max_x; i++){
        for(int j = 0; j < max_y; j++){
            b[num_of_bins].x = i;
            b[num_of_bins].y = j;
            b[num_of_bins].block_count = 0;
            b[num_of_bins].blocks = new int[number_of_pins];
            num_of_bins++;
        }
    }
    for(int i = 0; i < num_of_blocks; i++){
        int pid = p[i].pin_number;
        int p_index = i;
        int b_x = -1;
        int b_y = -1;
        int b_index;
        if(p[i].fixed){
            p_index = find_fixed_pin(pid);
            b_x = floor(fp[p_index].x_loc - 0.001 > 0 ? fp[p_index].x_loc - 0.001 : 0);
            b_y = floor(fp[p_index].y_loc - 0.001 > 0 ? fp[p_index].y_loc - 0.001 : 0);
        }
        else{
            b_x = floor(p[i].x_loc> 0 ? p[i].x_loc - 0.001 : 0);
            b_y = floor(p[i].y_loc> 0 ? p[i].y_loc - 0.001 : 0);
        }
        b_index = find_bin(b_x,b_y);
        b[b_index].blocks[b[b_index].block_count] = p[i].pin_number;
        b[b_index].block_count++;
    }
    /************************************************************Debug*****************************************************************/
    // for(int i = 0; i < num_of_bins; i++){
    //     if(b[i].block_count == 0)
    //         continue;
    //     cout << b[i].x << " " << b[i].y << ":" << endl;
    //     for(int j = 0; j < b[i].block_count; j++){
    //         int index = find_pin(b[i].blocks[j]);
    //         cout << b[i].blocks[j] << " : " << p[index].x_loc << " " << p[index].y_loc << endl; 
    //     }
    //     cout << "=================================================" << endl;
    // }
}

bool compare_two_bins(Bin b1, Bin b2){
    return b1.block_count < b2.block_count;
}

Bin* sort_overfilled_bins(){
    overfilled_bins_num = 0;
    for(int i = 0; i < num_of_bins; i++)
        if(b[i].block_count > 1)
            overfilled_bins_num++;
    Bin* overfilled_b = new Bin[overfilled_bins_num];
    overfilled_bins_num = 0;
    for(int i = 0; i < num_of_bins; i++)
        if(b[i].block_count > 1){
            overfilled_b[overfilled_bins_num] = b[i];
            overfilled_bins_num++;
        }
    sort(overfilled_b,overfilled_b+overfilled_bins_num,compare_two_bins);
    return overfilled_b;
}


void spread(){
    int itr = 0;
    int max_allowed_movement = pow(itr+1,2);
    overfilled_b = sort_overfilled_bins();
    for(int i = 0; i < overfilled_bins_num; i++){
        //TODO: find path and update
        
    }
}

//TODO: CHECK COMPUTE COST
int compute_cost(Bin tail_bin, Bin bk, int max_allowed_movement){
    double cost = DBL_MAX; 
    double distance = 0;
    int move = -1; //0 : up, 1 : down, 2 : right, 3: left
    if(tail_bin.x == bk.x){
        if(bk.y > tail_bin.y){ //right
            move = 2;
        }
        else{//left
            move = 3;
        }
    }
    else if(tail_bin.y == bk.y){
        if(bk.x > tail_bin.x) {//up
            move = 0;
        }
        else{//down
            move = 1;
        }
    }
    for(int i = 0; i < tail_bin.block_count; i++){//blocks placed in tail bin
        int pid = find_pin(tail_bin.blocks[i]);
        //compute x and y future positions to calculate quardaric distance
        double move_x = p[pid].moved_x_loc;
        double move_y = p[pid].moved_y_loc;
        switch (move)
        {
        case 0: move_x++;
            break;
        case 1: move_x--;
            break;
        case 2: move_y++;
            break;
        case 3: move_y--;
            break;
        default: 
            break;
        }
        distance = pow((p[pid].x_loc - move_x),2) + pow((p[pid].y_loc - move_y),2);
        if(distance < max_allowed_movement){
            if(distance < cost)
                cost = distance;
        }
    }
    return cost;
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
                if(input > max_x)
                    max_x = input;
                fp[num_of_fixed].x_loc = input;
                break;
            case 2:
                if(input > max_y)
                    max_y = input;
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
    
    make_clique();
    double* x_vec = solve(0);
    double* y_vec = solve(1);
    /************************************************************Debug*****************************************************************/
    // for(int i = 0 ; i < num_of_blocks - num_of_fixed; i++){
    //     cout << x_vec[i] << " " << y_vec[i] << endl; 
    // }
    int cnt = 0;
    for(int i = 0; i < num_of_blocks; i++){
        if(p[i].fixed)
            continue;
        else{
            p[i].x_loc = x_vec[cnt];
            p[i].y_loc = y_vec[cnt];
            p[i].moved_x_loc = x_vec[cnt];
            p[i].moved_y_loc = y_vec[cnt];
            cnt++;
        }
    }
    
    //free allocated spaces
    delete[] x_vec;
    delete[] y_vec;

    create_bin();
   


    return 0;
}