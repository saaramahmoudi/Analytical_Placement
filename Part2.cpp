#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <queue>
#include <limits.h>
#include <math.h>
#include <bits/stdc++.h>
#include<string>  
#include "umfpack.h"
#include "graphics.h"

using namespace std;

#define number_of_pins 510
#define number_of_nets 510
#define max_connection 100
#define max_fanout 20

void drawscreen (void);
void act_on_button_press (float x, float y);
void act_on_net_button_function (void (*drawscreen_ptr) (void));
void act_on_one_net_button_function (void (*drawscreen_ptr) (void));
void act_on_spread_button_function (void (*drawscreen_ptr) (void));
void act_on_solve_button_function (void (*drawscreen_ptr) (void));
void act_on_anchor_button_function (void (*drawscreen_ptr) (void));

struct Bin{
    int x;
    int y;
    int block_count;
    int* blocks;
    bool visited;
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
Pin* anchor_p;
Net* n_draw;
bool anchor_applied = false;
double anchor_weight = 100;
int n_draw_len = 0;
int num_of_blocks = 0;
int num_of_fixed = 0;
int num_of_bins = 0;
int overfilled_bins_num =0;
int max_x = -1;
int max_y = -1;
int nets_clicked = 0;
int one_net_clicked = -1;
int spread_clicked = 0;
int solve_clicked = 0;
int anchor_clicked = 0;

void init(){
    p = new Pin[number_of_pins];
    anchor_p = new Pin[number_of_pins];
    fp = new Fixed_pin[number_of_pins];
    n = new Net[number_of_nets];
    for(int i = 0; i < number_of_pins; i++){
        n[i].connected_block = new int[max_fanout];
        n[i].weight = new double[max_fanout];
        n[i].connected_block_len = 0;
    }
}

void make_clique(){
    //create clique model
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


double** cal_weight(){
    //Calculate Weight Matrix
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
    if(anchor_applied){
        int index = 0;
        for(int i = 0; i < num_of_blocks; i++){
            if(p[i].fixed)
                continue;
            else{
                if(dim < 1){
                    dim_vec[index] += anchor_weight * p[i].moved_x_loc;
                }
                else{
                    dim_vec[index] += anchor_weight * p[i].moved_y_loc;
                }
                index++;
            }
        }
    }
    
    /************************************************************Debug*****************************************************************/
    // for(int i = 0; i < num_of_blocks - num_of_fixed; i++){
    //     cout << " " << dim_vec[i] << " ";
    // }
    // cout << endl;
    return dim_vec;
}

double* solve(int dim, double** weight_mat){
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
    delete[] in_vector;
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
            b[num_of_bins].visited = false;
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
    /************************************************************Debug*****************************************************************/
    // cout << "/////////////////BINS/////////////////" << endl;
    // for(int i = 0; i < max_x*max_y; i++){
    //     if(b[i].block_count == 0){
    //         continue;
    //     }
    //     cout << b[i].x << " " << b[i].y << endl;
    //     cout << "CELL IN THIS BIN" << endl;
    //     for(int j = 0; j < b[i].block_count; j++){
    //         int pid_temp = find_pin(b[i].blocks[j]);
    //         int pid_fixed_temp = find_fixed_pin(b[i].blocks[j]);
    //         if(!p[pid_temp].fixed)
    //             cout << p[pid_temp].pin_number << " : " << p[pid_temp].moved_x_loc <<" " << p[pid_temp].moved_y_loc << endl;
    //         else
    //             cout << fp[pid_fixed_temp].pin_number << " : " << fp[pid_fixed_temp].x_loc <<" " << fp[pid_fixed_temp].y_loc << endl;
    //     }
    //     cout << "===============================" << endl;
    // }
    return overfilled_b;
}

pair<double,int> compute_cost(Bin tail_bin, Bin bk, int max_allowed_movement){
    double cost = DBL_MAX; 
    int pin_num_with_least_cost = -1;
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
        if(p[pid].fixed)
            continue;
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
        if(distance <= max_allowed_movement){
            if(distance < cost){
                cost = distance;
                pin_num_with_least_cost = p[pid].pin_number;
            }
        }
    }
    pair<double,int> cost_pair;
    cost_pair.first = cost;
    cost_pair.second = pin_num_with_least_cost;
    return cost_pair;
}

pair<double*,queue<Bin>*> find_path(int supply_bk, Bin bk, int max_allowed_movement){
    int demand = 0;
    int max_path = 500;
    queue<Bin>* path = new queue<Bin>[max_path];
    queue<Bin>* final_path = new queue<Bin>[supply_bk];
    double* final_path_cost = new double[supply_bk];
    double* path_cost = new double[max_path];
    for(int i = 0; i < max_path; i++)
        path_cost[i] = 0;
    int counter = 0;
    int counter_neighbor = 0;
    for(int i = 0; i < num_of_bins; i++)
        b[i].visited = false;
    int b_index = find_bin(bk.x,bk.y);
    b[b_index].visited = true;
    queue<queue<Bin>> q;
    //add bk into an empty path
    path[counter].push(bk);
    path_cost[counter] = 0;
    //add path into a FIFO queue
    q.push(path[counter]);
    while(demand < supply_bk && !q.empty()){  
        queue<Bin> current_path = q.front();
        q.pop();
        Bin tail_bin = current_path.back(); //tail bin
        /************************************************************Debug*****************************************************************/
        // cout << "+++++++++++++++++++++++++++++++source+++++++++++++++++++++++++++" << endl;
        // cout << bk.x << " " << bk.y << endl;
        // cout << "TAIL OF CURRENT PATH: " << tail_bin.x << " " << tail_bin.y << endl;
        //tail bin neighbors
        int cnt = 0;
        int b_neighbor = -1;
        while(cnt < 4 && demand < supply_bk){
            switch (cnt)
            {
            case 0: b_neighbor = find_bin(tail_bin.x+1,tail_bin.y);
                break;
            case 1: b_neighbor = find_bin(tail_bin.x-1,tail_bin.y);
                break;
            case 2: b_neighbor = find_bin(tail_bin.x,tail_bin.y+1);
                break;
            case 3: b_neighbor = find_bin(tail_bin.x,tail_bin.y-1);
                break;
            default:
                break;
            }
            cnt++;         
            //if a neighbor does not exist continue
            if(b_neighbor == -1)
                continue;
            //if it is already visited continue
            if(b[b_neighbor].visited)
                continue;
            pair<double,int> cost_res = compute_cost(tail_bin,b[b_neighbor],max_allowed_movement);
            double cost = cost_res.first;
            if(cost < DBL_MAX){//path cost is less than infinity
                queue<Bin> pcopy = current_path;
                pcopy.push(b[b_neighbor]);
                double cost_of_pcopy = path_cost[counter_neighbor] + cost;
                /************************************************************Debug*****************************************************************/
                // cout << "src: " << tail_bin.x << " " << tail_bin.y << endl;
                // cout << max_allowed_movement << endl;
                // cout << path_cost[counter_neighbor] << endl;
                // cout << cost << endl;
                // cout << b[b_neighbor].x << " " << b[b_neighbor].y << endl;
                // cout << "end of here" << endl;
                b[b_neighbor].visited = true;
                if(b[b_neighbor].block_count == 0){// b neighbor is empty
                    while(!pcopy.empty()){
                        Bin b_copy = pcopy.front();
                        pcopy.pop();
                        final_path[demand].push(b_copy);    
                    }
                    final_path_cost[demand] = cost_of_pcopy;
                    /************************************************************Debug*****************************************************************/
                    //     queue<Bin> temp = final_path[demand];
                    //     cout << "================================================================================" << endl;
                    //     cout << "source bin " << bk.x << " " << bk.y << endl;
                    //     cout << "neighbor: " << b[b_neighbor].x << " " << b[b_neighbor].y << endl;
                    //     cout << cost_of_pcopy << endl;
                    //     cout << "size of path: " << temp.size() << endl;
                    //     while(!temp.empty()){
                    //         Bin b_print = temp.front();
                    //         temp.pop();
                    //         cout << "bin location in path " << b_print.x << " " << b_print.y << endl;
                    // }
                    demand++;
                }
                else{
                    q.push(pcopy);
                    counter++;
                    path_cost[counter] = cost;  
                }
            }
        }
        counter_neighbor++;
    }
    pair<double*,queue<Bin>*> res;
    res.second = final_path;
    res.first = final_path_cost;
    delete[] path_cost;
    delete[] path;

    /************************************************************Debug*****************************************************************/
    // print all the possible path
    // cout << "!!!path for : " << bk.x << " " << bk.y << " !!!" << endl;
    // cout << "number of path found: " << demand << endl;
    // cout << "max movement: " << max_allowed_movement << endl;
    // double* temp2 = res.first;
    // queue<Bin>* temp1 = res.second;
    // for(int i = 0; i < demand; i++){
    //     queue<Bin> temp = temp1[i];
    //     cout << "=============================================" << endl;
    //     while(!temp.empty()){
    //         Bin b_print = temp.front();
    //         temp.pop();
    //         cout << "bin location in path " << b_print.x << " " << b_print.y << endl;
    //     }
    //     cout << "cost " << temp2[i] << endl;
    // }  
    // cout << "=============================================" << endl;  
    return res;
}

void sort_queue(double* &q, queue<Bin>* &q_parallel,int n){
    for (int i = 0; i < n-1; i++){    
        for (int j = 0; j < n-i-1; j++){
            if(q_parallel[j].empty() || q_parallel[j+1].empty())
                continue;
            if (q[j] > q[j+1]){
                double temp = q[j];
                q[j] = q[j+1];
                q[j+1] = temp;
                queue <Bin> temp1 = q_parallel[j];
                q_parallel[j] = q_parallel[j+1];
                q_parallel[j+1] = temp1;
            }
        }
    }
}

void reverse_queue(queue<Bin>& Queue)
{
    stack<Bin> Stack;
    while (!Queue.empty()) {
        Stack.push(Queue.front());
        Queue.pop();
    }
    while (!Stack.empty()) {
        Queue.push(Stack.top());
        Stack.pop();
    }
}

int find_cell_in_bin(int pin_number,int bin_number){
    for(int i = 0; i < b[bin_number].block_count; i++){
        if(b[bin_number].blocks[i] == pin_number){
            return i;
        }
    }
    return -1;
}

void remove_cell_from_bin(int pin_number,int bin_number){
    int index = find_cell_in_bin(pin_number,bin_number);
    for(int i = index; i < b[bin_number].block_count; i++){
        b[bin_number].blocks[i] = b[bin_number].blocks[i+1];
    }
    b[bin_number].block_count--;
}

void spread(){
    //must go inside a loop untill no overfilled bin is available
    int itr = 0;
    overfilled_b = sort_overfilled_bins();
    while(overfilled_bins_num != 0){//there is no empty bin left!
        int max_allowed_movement = pow(itr+1,2);
        for(int i = 0; i < overfilled_bins_num; i++){
            int overfilled_b_id = find_bin(overfilled_b[i].x,overfilled_b[i].y);
            pair<double*,queue<Bin>*> possible_path_pair = find_path(b[overfilled_b_id].block_count-1,overfilled_b[i],max_allowed_movement);
            queue<Bin>* possible_path = possible_path_pair.second;
            double* possible_path_cost = possible_path_pair.first;
            sort_queue(possible_path_cost,possible_path,b[overfilled_b_id].block_count-1);
            int index = 0;  
             /************************************************************Debug*****************************************************************/
            // cout << "========================BIN ID========================" << endl;
            // cout << "bin loc: " << b[overfilled_b_id].x << " " << b[overfilled_b_id].y << endl; 
            // cout << "Max Movement: " << max_allowed_movement << endl;
            // cout << "Demand " << b[overfilled_b_id].block_count-1 << endl;
            
            for(int j = 0; j < overfilled_b[i].block_count-1;j++){//loop over possible path
                /************************************************************Debug*****************************************************************/
                //print possible path!
                // int inp = 0;
                // cout << "Bin loc: " << overfilled_b[i].x << " " << overfilled_b[i].y << endl;
                // cout << "Number of path " << b[overfilled_b_id].block_count-1 << endl;
                // for(int i = 0; i < b[overfilled_b_id].block_count-1; i++){
                //     queue<Bin> possible_path_copy = possible_path[inp];
                //     cout << "=========================================================" << endl;
                //     while (!possible_path_copy.empty()){
                //         Bin b_print = possible_path_copy.front();
                //         cout << "bin location in path " << b_print.x << " " << b_print.y << endl;
                //         possible_path_copy.pop();
                //     }
                // }
                // cout << "========================================================" << endl;
                // inp++;
                if(b[overfilled_b_id].block_count > 1){//bin is still overfilled
                    // cout << b[overfilled_b_id].block_count << endl;      
                    queue<Bin> possible_path_copy = possible_path[index];
                    if(possible_path_copy.empty()){//if no path has been found break the loop
                        break;
                    }
                    bool valid_path = true; //cells can be moved without violating max_allowed_movement along the path
                    stack<Bin> bin_stack;
                    Bin b_source = possible_path_copy.front();
                    possible_path_copy.pop();
                    bin_stack.push(b_source);
                    // cout << "path size: " << possible_path_copy.size() << endl;
                    while(!possible_path_copy.empty()){
                        Bin b_sink = possible_path_copy.front();
                        possible_path_copy.pop();
                        Bin b_src = bin_stack.top();
                        bin_stack.pop();
                        pair<double,int> cost_revised = compute_cost(b_src,b_sink,max_allowed_movement);
                        double cost = cost_revised.first;
                        /************************************************************Debug*****************************************************************/
                        // cout << "Source: " << b_src.x << " " << b_src.y << endl;
                        // cout << "Target: " << b_sink.x << " " << b_sink.y << endl;
                        // cout << "Cost: " << cost << endl;
                        if(cost == DBL_MAX || cost_revised.second <= 0){
                            /************************************************************Debug*****************************************************************/
                            // cout << "INVALID PATH CATCH" << endl;
                            // cout << index << endl;
                            // cout << "=========" << endl;
                            valid_path = false;
                            break;
                        }
                        bin_stack.push(b_sink);
                    }
                    if(valid_path){//start moving cell along path[index]
                        possible_path_copy = possible_path[index];
                        reverse_queue(possible_path_copy);
                        Bin b_sink = possible_path_copy.front();
                        possible_path_copy.pop();
                        while(!possible_path_copy.empty()){
                            Bin b_src = possible_path_copy.front();
                            possible_path_copy.pop();
                            pair<double,int> cost_revised = compute_cost(b_src,b_sink,max_allowed_movement);
                            if(cost_revised.second <= 0){
                                cout << "ERROR: No cell found to move!!" << endl;
                            }
                            else{
                                int pid = find_pin(cost_revised.second);
                                int b_src_id = find_bin(b_src.x,b_src.y);
                                int b_sink_id = find_bin(b_sink.x,b_sink.y);
                                // cout << "=====================BMOVE==============" << endl;
                                // cout << "moving cell " << p[pid].pin_number <<" " << cost_revised.second << endl;
                                // cout << "source: " << b_src.x << " " << b_src.y << endl; 
                                // cout << "source block count: " << b[b_src_id].block_count << endl;
                                // for(int i = 0; i < b[b_src_id].block_count; i++){
                                //     cout << b[b_src_id].blocks[i] << " ";
                                // }
                                // cout << endl;
                                // cout << "sink: " << b_sink.x << " " << b_sink.y << endl;
                                // cout << "sink block count: " << b[b_sink_id].block_count << endl;
                                // for(int j =0 ;j < b[b_sink_id].block_count;j++){
                                //     cout << b[b_sink_id].blocks[j] << " ";
                                // }
                                // cout << endl;
                                remove_cell_from_bin(p[pid].pin_number,b_src_id);
                                b[b_sink_id].blocks[b[b_sink_id].block_count] = p[pid].pin_number;
                                b[b_sink_id].block_count++;
                                /************************************************************Debug*****************************************************************/     
                                // cout << "=====================MOVE==============" << endl;
                                // cout << "moving cell " << p[pid].pin_number << endl;
                                // cout << "source: " << b_src.x << " " << b_src.y << endl; 
                                // cout << "source block count: " << b[b_src_id].block_count << endl;
                                // for(int i = 0; i < b[b_src_id].block_count; i++){
                                //     cout << b[b_src_id].blocks[i] << " ";
                                // }
                                // cout << endl;
                                // cout << "sink: " << b_sink.x << " " << b_sink.y << endl;
                                // cout << "sink block count: " << b[b_sink_id].block_count << endl;
                                // for(int j =0 ;j < b[b_sink_id].block_count;j++){
                                //     cout << b[b_sink_id].blocks[j] << " ";
                                // }
                                // cout << endl;
                                if(b_src.x == b_sink.x){
                                    if(b_sink.y > b_src.y){ //right
                                        p[pid].moved_y_loc++;
                                    }
                                    else{//left
                                        p[pid].moved_y_loc--;
                                    }
                                }
                                else if(b_src.y == b_sink.y){
                                    if(b_sink.x > b_src.x) {//up
                                        p[pid].moved_x_loc++;
                                    }
                                    else{//down
                                        p[pid].moved_x_loc--;
                                    }
                                }
                            }
                            b_sink = b_src;
                        }
                    }
                    if(index > overfilled_b[i].block_count-1){//possible path ended
                        break;
                    }
                    index++;
                }
            }
        }
        itr++;
        overfilled_b = sort_overfilled_bins();
        // cout << overfilled_bins_num<< endl;
    }
    // cout << "done" << endl;
}

double** add_anchor(double** weight_mat,double anchor_weight){
    double** anchor_weight_mat = new double*[num_of_blocks - num_of_fixed];
    for(int i = 0; i < num_of_blocks-num_of_fixed; i++)
        anchor_weight_mat[i] = new double[num_of_blocks - num_of_fixed]; 
    for(int i = 0; i < num_of_blocks - num_of_fixed; i++){
        for(int j = 0; j < num_of_blocks-num_of_fixed; j++){
            if(i != j){
                anchor_weight_mat[i][j] = weight_mat[i][j];
            }
            else{
                anchor_weight_mat[i][j] = weight_mat[i][j] + anchor_weight;
            }
        }
    }
    return anchor_weight_mat;
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
            double x_loc = (!p[pid_temp].fixed) ? (p[pid_temp].moved_x_loc) : (fp[pid_fixed_temp].x_loc); 
            double y_loc = (!p[pid_temp].fixed) ? (p[pid_temp].moved_y_loc) : (fp[pid_fixed_temp].y_loc); 
            if(anchor_applied && !p[pid_temp].fixed){
                x_loc = anchor_p[pid_temp].x_loc;
                y_loc = anchor_p[pid_temp].y_loc;
            }
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

double cal_displacement(){
    double displacement_value = 0;
    for(int i = 0; i < num_of_blocks; i++){
        if(p[i].fixed)
            continue;
        displacement_value += (abs(p[i].moved_x_loc - p[i].x_loc)) + abs((p[i].moved_y_loc - p[i].y_loc));
    }
    return displacement_value;
}

void create_nets_for_drawing(){
    for(int i = 0; i < number_of_nets; i++){
        if(n[i].connected_block_len != 0)
            n_draw_len++;
    }
    n_draw = new Net[n_draw_len];
    int cnt =0;
    for(int i = 0; i < number_of_nets; i++){
        if(n[i].connected_block_len != 0){
            n_draw[cnt] = n[i];
            cnt++;
        }
    }
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
                anchor_p[num_of_blocks].pin_number = input;
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
    double** weight_mat = cal_weight();
    double* x_vec = solve(0,weight_mat);
    double* y_vec = solve(1,weight_mat);
    /************************************************************Debug*****************************************************************/
    // for(int i = 0 ; i < num_of_blocks - num_of_fixed; i++){
    //     cout << x_vec[i] << " " << y_vec[i] << endl; 
    // }
    int cnt = 0;
    for(int i = 0; i < num_of_blocks; i++){
        if(p[i].fixed){
            continue;
        }
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
    double HPWL_before_spread = cal_HPWL();
    cout << "HPWL: " << HPWL_before_spread << endl;
    create_bin();
    spread();
    double HPWL_after_spread = cal_HPWL();
    cout << "HPWL AFTER SPREAD: " << HPWL_after_spread << endl;
    double displacement_value = cal_displacement();
    cout << "Displacement Value: " << displacement_value << endl; 
    anchor_applied = true;
    double** anchor_weight_mat = add_anchor(weight_mat,anchor_weight);
    double* x_vec1 = solve(0,anchor_weight_mat);
    double* y_vec1 = solve(1,anchor_weight_mat);  
    cnt = 0;
    for(int i = 0; i < num_of_blocks; i++){
        if(p[i].fixed){
            continue;
        }
        else{
            anchor_p[i].x_loc = x_vec1[cnt];
            anchor_p[i].y_loc = y_vec1[cnt];
            cnt++;
        }
    }
    double HPWL_after_anchor_applied = cal_HPWL();
    cout << "HPWL After Anchor Applied " << HPWL_after_anchor_applied  << endl;
    delete[] x_vec1;
    delete[] y_vec1;
    create_nets_for_drawing();
    init_graphics("Assignment 2 - Part2", WHITE);
    init_world (0.,0.,5000.,5000.);
    update_message("Fatemehsadat(Sara) Mahmoudi - Placement");
    create_button ("Proceed", "Solve", act_on_solve_button_function);
    create_button ("Solve", "Spread", act_on_spread_button_function);
    create_button ("Spread", "Anchor", act_on_anchor_button_function);
    create_button ("Anchor", "NETS", act_on_net_button_function);
    create_button ("NETS", "1 NET", act_on_one_net_button_function);
    event_loop(act_on_button_press, NULL, NULL, drawscreen); 
    return 0;
}

void connect_two_pin(int pid1, int pid2){
    float left_corner_x = 300;
    float left_corner_y = 2000;
    float len = 300;
    float x_bin1,y_bin1,x_margin1,y_margin1;
    float x_bin2,y_bin2,x_margin2,y_margin2; 
    float p1_x, p1_y, p2_x, p2_y;
    if(solve_clicked!= 0){
        if(!p[pid1].fixed){
            p1_x = p[pid1].x_loc;
            p1_y = p[pid1].y_loc;
            x_bin1 = floor(p[pid1].x_loc);
            y_bin1 = floor(p[pid1].y_loc);
            x_margin1 = p[pid1].x_loc - x_bin1;
            y_margin1 = p[pid1].y_loc - y_bin1;
        }
        else{
            int fp1 = find_fixed_pin(p[pid1].pin_number);
            p1_x = fp[fp1].x_loc;
            p1_y = fp[fp1].y_loc - 1;
            x_bin1 = fp[fp1].x_loc;
            y_bin1 = fp[fp1].y_loc -1 ;
            x_margin1 = 0;
            y_margin1 = 0; 
            }
            if(!p[pid2].fixed){
                p2_x = p[pid2].x_loc;
                p2_y = p[pid2].y_loc;
                x_bin2 = floor(p[pid2].x_loc);
                y_bin2 = floor(p[pid2].y_loc);
                x_margin2 = p[pid2].x_loc - x_bin2;
                y_margin2 = p[pid2].y_loc - y_bin2;
            }
            else{
                int fp2 = find_fixed_pin(p[pid2].pin_number);
                p2_x = fp[fp2].x_loc;
                p2_y = fp[fp2].y_loc - 1;
                x_bin2 = fp[fp2].x_loc;
                y_bin2 = fp[fp2].y_loc - 1;
                x_margin2 = 0;
                y_margin2 = 0; 
            }
    }
    else if (spread_clicked != 0){
        if(!p[pid1].fixed){
            p1_x = p[pid1].moved_x_loc;
            p1_y = p[pid1].moved_y_loc;
            x_bin1 = floor(p[pid1].moved_x_loc);
            y_bin1 = floor(p[pid1].moved_y_loc);
            x_margin1 = p[pid1].moved_x_loc - x_bin1;
            y_margin1 = p[pid1].moved_y_loc - y_bin1;
        }
        else{
            int fp1 = find_fixed_pin(p[pid1].pin_number);
            p1_x = fp[fp1].x_loc;
            p1_y = fp[fp1].y_loc - 1;
            x_bin1 = fp[fp1].x_loc;
            y_bin1 = fp[fp1].y_loc -1 ;
            x_margin1 = 0;
            y_margin1 = 0; 
            }
            if(!p[pid2].fixed){
                p2_x = p[pid2].moved_x_loc;
                p2_y = p[pid2].moved_y_loc;
                x_bin2 = floor(p[pid2].moved_x_loc);
                y_bin2 = floor(p[pid2].moved_y_loc);
                x_margin2 = p[pid2].moved_x_loc - x_bin2;
                y_margin2 = p[pid2].moved_y_loc - y_bin2;
            }
            else{
                int fp2 = find_fixed_pin(p[pid2].pin_number);
                p2_x = fp[fp2].x_loc;
                p2_y = fp[fp2].y_loc - 1;
                x_bin2 = fp[fp2].x_loc;
                y_bin2 = fp[fp2].y_loc - 1;
                x_margin2 = 0;
                y_margin2 = 0; 
            }
    }
    else if (anchor_clicked != 0){
        if(!p[pid1].fixed){
            p1_x = anchor_p[pid1].x_loc;
            p1_y = anchor_p[pid1].y_loc;
            x_bin1 = floor(anchor_p[pid1].x_loc);
            y_bin1 = floor(anchor_p[pid1].y_loc);
            x_margin1 = anchor_p[pid1].x_loc - x_bin1;
            y_margin1 = anchor_p[pid1].y_loc - y_bin1;
        }
        else{
            int fp1 = find_fixed_pin(p[pid1].pin_number);
            p1_x = fp[fp1].x_loc;
            p1_y = fp[fp1].y_loc - 1;
            x_bin1 = fp[fp1].x_loc;
            y_bin1 = fp[fp1].y_loc -1 ;
            x_margin1 = 0;
            y_margin1 = 0; 
            }
            if(!p[pid2].fixed){
                p2_x = anchor_p[pid2].x_loc;
                p2_y = anchor_p[pid2].y_loc;
                x_bin2 = floor(anchor_p[pid2].x_loc);
                y_bin2 = floor(anchor_p[pid2].y_loc);
                x_margin2 = anchor_p[pid2].x_loc - x_bin2;
                y_margin2 = anchor_p[pid2].y_loc - y_bin2;
            }
            else{
                int fp2 = find_fixed_pin(p[pid2].pin_number);
                p2_x = fp[fp2].x_loc;
                p2_y = fp[fp2].y_loc - 1;
                x_bin2 = fp[fp2].x_loc;
                y_bin2 = fp[fp2].y_loc - 1;
                x_margin2 = 0;
                y_margin2 = 0; 
            }
    }
    //locate their positions
    if(p1_x < p2_x && abs(p1_y - p2_y) > 0.5 * len){ //right
        drawline(left_corner_x + x_bin1*len + len*x_margin1 + len/20, left_corner_y + (max_y-2-y_bin1)*len + y_margin1*len,left_corner_x + x_bin2*len + len*x_margin2 -len/20 , left_corner_y + (max_y-2-y_bin2)*len +y_margin2*len);
    }
    else if(p1_x > p2_x && abs(p1_y - p2_y) > 0.5 * len){//left
        drawline(left_corner_x + x_bin1*len + len*x_margin1 - len/20, left_corner_y + (max_y-2-y_bin1)*len + y_margin1*len,left_corner_x + x_bin2*len + len*x_margin2 +len/20 , left_corner_y + (max_y-2-y_bin2)*len +y_margin2*len);
    }
    else if(p1_y < p2_y){//up
        drawline(left_corner_x + x_bin1*len + len*x_margin1 , left_corner_y + (max_y-2-y_bin1)*len + y_margin1*len - len/20,left_corner_x + x_bin2*len + len*x_margin2  , left_corner_y + (max_y-2-y_bin2)*len + + y_margin2*len +len/20);
    }
    else{
        drawline(left_corner_x + x_bin1*len + len*x_margin1 , left_corner_y + (max_y-2-y_bin1)*len + y_margin1*len + len/20,left_corner_x + x_bin2*len + len*x_margin2  , left_corner_y + (max_y-2-y_bin2)*len + + y_margin2*len -len/20);
    }
}

void draw_nets(){
    for(int i = 0; i < number_of_nets; i++){
        if(n[i].connected_block_len == 0){
            continue;
        }
        if(n[i].connected_block_len == 2){
            setcolor(LIGHTGREY);
            int pid1 = find_pin(n[i].connected_block[0]);
            int pid2 = find_pin(n[i].connected_block[1]); 
            connect_two_pin(pid1,pid2);
        }
        else{//clique model
            setcolor(YELLOW);
            for(int j = 0; j < n[i].connected_block_len; j++){
                for(int k = j+1; k < n[i].connected_block_len; k++){
                    int pid1 = find_pin(n[i].connected_block[j]);
                    int pid2 = find_pin(n[i].connected_block[k]);
                    connect_two_pin(pid1,pid2);
                }
            }
        }
    }
}
void draw_one_net(){
    if(n_draw[one_net_clicked].connected_block_len == 2){
        setcolor(LIGHTGREY);
        int pid1 = find_pin(n_draw[one_net_clicked].connected_block[0]);
        int pid2 = find_pin(n_draw[one_net_clicked].connected_block[1]); 
        connect_two_pin(pid1,pid2);
    }
    else{//clique model
        setcolor(YELLOW);
        for(int j = 0; j < n_draw[one_net_clicked].connected_block_len; j++){
            for(int k = j+1; k < n_draw[one_net_clicked].connected_block_len; k++){
                int pid1 = find_pin(n_draw[one_net_clicked].connected_block[j]);
                int pid2 = find_pin(n_draw[one_net_clicked].connected_block[k]);
                connect_two_pin(pid1,pid2);
            }
        }
    }
}

void spread_draw(){
    float left_corner_x = 300;
    float left_corner_y = 2000;
    float len = 300;
    //DRAW OTHER ELEMENTS
    setcolor(BLACK);
    for(int i = 0; i < num_of_blocks; i++){
        if(p[i].fixed)
            continue;    
        float x_bin = floor(p[i].moved_x_loc);
        float y_bin = floor(p[i].moved_y_loc);
        float x_margin = p[i].moved_x_loc - x_bin;
        float y_margin = p[i].moved_y_loc - y_bin;
        fillarc(left_corner_x + x_bin*len + len*x_margin , left_corner_y + (max_y-2-y_bin)*len + y_margin*len ,len/20,0,360);
        char* s = &(to_string(p[i].pin_number)[0]);
        drawtext (left_corner_x + x_bin*len + len*x_margin, left_corner_y + (max_y-2-y_bin)*len + y_margin*len + len/12,s, 300);
    }
}

void solve_draw(){
    float left_corner_x = 300;
    float left_corner_y = 2000;
    float len = 300;
    //DRAW OTHER ELEMENTS
    setcolor(BLACK);
    for(int i = 0; i < num_of_blocks; i++){
        if(p[i].fixed)
            continue;    
        float x_bin = floor(p[i].x_loc);
        float y_bin = floor(p[i].y_loc);
        float x_margin = p[i].x_loc - x_bin;
        float y_margin = p[i].y_loc - y_bin;
        fillarc(left_corner_x + x_bin*len + len*x_margin , left_corner_y + (max_y-2-y_bin)*len + y_margin*len ,len/20,0,360);
        char* s = &(to_string(p[i].pin_number)[0]);
        drawtext (left_corner_x + x_bin*len + len*x_margin, left_corner_y + (max_y-2-y_bin)*len + y_margin*len + len/12,s, 300);
    }
}

void anchor_draw(){
    float left_corner_x = 300;
    float left_corner_y = 2000;
    float len = 300;
    //DRAW OTHER ELEMENTS
    setcolor(BLACK);
    for(int i = 0; i < num_of_blocks; i++){
        if(p[i].fixed)
            continue;    
        float x_bin = floor(anchor_p[i].x_loc);
        float y_bin = floor(anchor_p[i].y_loc);
        float x_margin = anchor_p[i].x_loc - x_bin;
        float y_margin = anchor_p[i].y_loc - y_bin;
        fillarc(left_corner_x + x_bin*len + len*x_margin , left_corner_y + (max_y-2-y_bin)*len + y_margin*len ,len/20,0,360);
        char* s = &(to_string(anchor_p[i].pin_number)[0]);
        drawtext (left_corner_x + x_bin*len + len*x_margin, left_corner_y + (max_y-2-y_bin)*len + y_margin*len + len/12,s, 300);
    }
}

void drawscreen (void) {
    set_draw_mode (DRAW_NORMAL);
    clearscreen(); 
    setfontsize (2);
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
        char* s = &(to_string(fp[i].pin_number)[0]);
        drawtext (left_corner_x + x_margin*len, left_corner_y + (max_y-1-y_margin)*len + len/12,s, 300);
    }

    //DRAW NETS
    if(solve_clicked != 0){
        solve_draw();
    }
    if(spread_clicked != 0){
        spread_draw();
    }
    if(anchor_clicked != 0){
        anchor_draw();
    }
    if(nets_clicked % 2 == 1 && (spread_clicked != 0 || solve_clicked!= 0 || anchor_clicked != 0)){
        draw_nets(); 
    }
    if(one_net_clicked != n_draw_len && one_net_clicked != -1 && (spread_clicked != 0 || solve_clicked!= 0 || anchor_clicked!= 0)){
        draw_one_net();
    }
}
void act_on_button_press (float x, float y) {
    printf("User clicked a button at coordinates (%f, %f)\n", x, y);
}
void act_on_spread_button_function (void (*drawscreen_ptr) (void)){
    if(solve_clicked != 0){
        solve_clicked = 0;
        one_net_clicked = -1;
        drawscreen();
    }
    if(anchor_clicked != 0){
        anchor_clicked = 0;
        one_net_clicked = -1;
        drawscreen();
    }
    if(spread_clicked == 0){
        spread_clicked++;
        spread_draw();
    }
    else{
        spread_clicked = 0;
        drawscreen();
    }
}
void act_on_anchor_button_function (void (*drawscreen_ptr) (void)){
    if(spread_clicked != 0){
        spread_clicked = 0;
        one_net_clicked = -1;
        drawscreen();
    }
    if(solve_clicked != 0){
        solve_clicked = 0;
        one_net_clicked = -1;
        drawscreen();
    }
    if(anchor_clicked == 0){
        anchor_clicked++;
        anchor_draw();
    }
    else{
        anchor_clicked = 0;
        drawscreen();
    }
}
void act_on_solve_button_function (void (*drawscreen_ptr) (void)){
    if(spread_clicked != 0){
        spread_clicked = 0;
        one_net_clicked = -1;
        drawscreen();
    }
    if(anchor_clicked != 0){
        anchor_clicked = 0;
        one_net_clicked = -1;
        drawscreen();
    }
    if(solve_clicked == 0){
        solve_clicked++;
        solve_draw();
    }
    else{
        solve_clicked = 0;
        drawscreen();
    }
}
void act_on_net_button_function (void (*drawscreen_ptr) (void)) {
    //DRAW NETS
    one_net_clicked = -1;
    drawscreen();
    nets_clicked++;
    if(nets_clicked % 2 == 1){
        draw_nets();
    }
    else{
        drawscreen();
    }
}
void act_on_one_net_button_function (void (*drawscreen_ptr) (void)){
    nets_clicked = 0;
    drawscreen();
    one_net_clicked++;
    if(one_net_clicked == n_draw_len){
        one_net_clicked = 0;
        drawscreen();
    }
    else{
        drawscreen();   
        draw_one_net();
    }        
}
