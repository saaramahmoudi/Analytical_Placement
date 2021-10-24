#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <queue>
#include <limits.h>
#include <math.h>
#include <bits/stdc++.h>
#include "umfpack.h"
#include "graphics.h"

using namespace std;

#define number_of_pins 510
#define number_of_nets 510
#define max_connection 100
#define max_fanout 20

void drawscreen (void);
void act_on_button_press (float x, float y);

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
    double x_loc;
    double y_loc;
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


double cal_HPWL(){
    double HPWL = 0; 
    double p_x_max = DBL_MIN;
    double p_y_max = DBL_MIN; 
    double p_x_min = DBL_MAX;
    double p_y_min = DBL_MAX;
    for(int i = 0; i < number_of_nets; i++){
        if(n[i].connected_block_len == 0)
            continue;
        for(int j = 0; j < n[i].connected_block_len ; j++){
            int pid_temp = find_pin(n[i].connected_block[j]);
            int pid_fixed_temp = find_fixed_pin(n[i].connected_block[j]);
            double x_loc = (!p[pid_temp].fixed) ? (p[pid_temp].x_loc) : (fp[pid_fixed_temp].x_loc); 
            double y_loc = (!p[pid_temp].fixed) ? (p[pid_temp].y_loc) : (fp[pid_fixed_temp].y_loc); 
            // cout << "net number: " << i << " connected to pin : " << p[pid_temp].pin_number << endl;
            // cout << "pin loc : " << x_loc << " " << y_loc << endl;
            // cout << "========================================================" << endl; 
            if(p_x_max < x_loc){
                p_x_max = x_loc;
            }
            if(p_y_max < y_loc){
                p_y_max = y_loc;
            }
            if(p_x_min > x_loc){
                p_x_min = x_loc;
            }
            if(p_y_min > y_loc){
                p_y_min = y_loc;
            }
        }
        HPWL += (p_y_max - p_y_min) + (p_x_max-p_x_min);
        p_x_max = DBL_MIN;
        p_y_max = DBL_MIN; 
        p_x_min = DBL_MAX;
        p_y_min = DBL_MAX;
    }
    return HPWL;
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
    for(int i = 0; i < num_of_fixed; i++){
        cout << fp[i].pin_number << " " << fp[i].x_loc << " " << fp[i].y_loc << endl;
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
            cnt++;
        }
    }
    delete[] x_vec;
    delete[] y_vec;
    double HPWL = cal_HPWL();
    cout << "HPWL: " << HPWL << endl;

    init_graphics("Assignment 2 - Part1", WHITE);
    init_world (0.,0.,5000.,5000.);
    update_message("Fatemehsadat(Sara) Mahmoudi - Placement");
    event_loop(act_on_button_press, NULL, NULL, drawscreen); 
    return 0;
}


void drawscreen (void) {
    set_draw_mode (DRAW_NORMAL);
    clearscreen(); 
    setfontsize (10);
    setlinestyle (SOLID);
    setlinewidth (2);
    setcolor (BLACK);

    //DRAW GRID
    setcolor(CYAN);
    float left_corner_x = 300;
    float left_corner_y = 2000;
    float len = 300;
    for(int i = 0; i < max_x; i++){
        for(int j = 0; j < max_y; j++){        
            drawrect (left_corner_x + j*len, left_corner_y + i*len, left_corner_x + j*len + len, left_corner_y + i*len - len); 
        }
    }

    //DRAW FIXED ELEMENTS
    setcolor(MAGENTA);
    for(int i = 0; i < num_of_fixed; i++){
        float x_margin = fp[i].x_loc;
        float y_margin = fp[i].y_loc;
        fillarc(left_corner_x + x_margin*len , left_corner_y + (max_y-1-y_margin)*len ,len/15,0,360);
    }
    
    //DRAW OTHER ELEMENTS
    setcolor(BLUE);
    for(int i = 0; i < num_of_blocks; i++){
        if(p[i].fixed)
            continue;    
        float x_bin = floor(p[i].x_loc);
        float y_bin = floor(p[i].y_loc);
        float x_margin = p[i].x_loc - x_bin;
        float y_margin = p[i].y_loc - y_bin;
        fillarc(left_corner_x + x_bin*len + len*x_margin , left_corner_y + (max_y-2-y_bin)*len + y_margin*len ,len/20,0,360);
    }

    //DRAW NETS
    setcolor(BLACK);
}


void act_on_button_press (float x, float y) {
    printf("User clicked a button at coordinates (%f, %f)\n", x, y);
}