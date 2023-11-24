#include <barneshut.h>
using namespace std;

static int points_checked;

void tree_construct(unordered_map<uint64_t, struct TreeNode>& tree, vector<particle> particles, int n_val){
    //unordered_map<int, struct TreeNode> tree;
    particle p = particles[0];
    struct TreeNode root = {0, 0, MIN_X, MIN_Y, MAX_X, MAX_Y, true, p.px, p.py, p.mass, 0};
    tree.insert(make_pair(0, root));
    for (int i = 1; i < n_val; i++){
        insert_node(tree, particles, i, 0);
    }
    /*printf("tree constructed!\n");
    for(auto const &pair : tree){
        print_node(pair.second);
    }*/
}

void insert_node(unordered_map<uint64_t, struct TreeNode>& tree, vector<particle> particles, int index,  uint64_t key){
    
    if (particles[index].mass == -1){
        return;
    }
    if(tree.find(key) == tree.end()){
        printf("tree node not found\n");
        return;
    }
    TreeNode *node = &tree[key];
    if (node->isExternal){
        //printf("external node layer %d, index %d, p_index %d\n", node->layer, node->index_in_layer, index);
 
        node->isExternal = false;
        int old_index = node->index;
        node->cordX = 0;
        node->cordY = 0;
        node->mass = 0;
        node->index = 0;
        insert_node(tree, particles, old_index, key);
        insert_node(tree, particles, index, key);
        return;
    }
    if(!node->isExternal){
        //printf("non external node layer %d, index %d, p_index %d\n", node->layer, node->index_in_layer, index);
        
        node->cordX += particles[index].px * particles[index].mass;
        node->cordY += particles[index].py * particles[index].mass;
        node->mass += particles[index].mass;
        uint64_t child_key;
        int re;
        tie(child_key , re) = find_child(tree, particles[index], key);
        //printf("child key %lld ",(long long)child_key);
        if (re == -1){
            //printf("new child!\n");
            return;
        }else if(re < -1){
            printf("find child error! index = %d\n", index);
            
            return;
        }
        
        insert_node(tree, particles, index, child_key);
        return;
    }
}

tuple<uint64_t, int> find_child(unordered_map<uint64_t, struct TreeNode>& tree, particle p, uint64_t key){
    double minX=0, minY=0, maxX=0, maxY=0;
    //printf("find child\n");
    if (tree.find(key) == tree.end()){
        return make_tuple(0, -2);
    }
    TreeNode node = tree[key];
    uint16_t layer = node.layer + 1; 
    uint32_t index_in_layer = node.index_in_layer;
    for (int i = 0; i < 4; i++){
        uint32_t child_index = index_in_layer * 4 + i;
        uint64_t child_layer = (uint64_t)layer;
        uint64_t child_key =(child_layer<<32)| child_index;
        //printf("layer %d, u64 layer %lld, child_key %lld",layer, (long long)child_layer,(long long)child_key);
        switch (i){
            case 0:
                minX = (node.minX + node.maxX)/2;
                maxX = node.maxX;
                minY = (node.minY + node.maxY)/2;
                maxY = node.maxY;
                break;
            case 1:
                minX = node.minX;
                maxX = (node.minX + node.maxX)/2;
                minY = (node.minY + node.maxY)/2;
                maxY = node.maxY;
                break;
            case 2:
                minX = node.minX;
                maxX = (node.minX + node.maxX)/2;
                minY = node.minY;
                maxY = (node.minY + node.maxY)/2;
                break;
            case 3:
                minX = (node.minX + node.maxX)/2;
                maxX = node.maxX;
                minY = node.minY;
                maxY = (node.minY + node.maxY)/2;
                break;
        }
        if (p.px >= minX && p.px <= maxX &&
            p.py >= minY && p.py <= maxY){
            if (tree.find(child_key) == tree.end()){
                struct TreeNode child = {layer, child_index, minX, minY, maxX, maxY, true, p.px, p.py, p.mass, p.index};
                tree.emplace(child_key, child);
                //print_node(tree.at(child_key));
                return make_tuple(child_key, -1);
            }else{
                
                return make_tuple(child_key, 0);
            }
        }
        
    } 
    print_node(node);
    return make_tuple(0, -2);    
        
    
}

void print_node(struct TreeNode node){
    printf("layer %d; index_in_layer %d, minX %f; minY %f; maxX %f; maxY %f; isExtern %d; cordX %f; cordY %f; mass %f; index %d \n",
    node.layer, node.index_in_layer, node.minX, node.minY, node.maxX, node.maxY, node.isExternal,  node.cordX, node.cordY, node.mass, node.index);
}

void new_position(unordered_map<uint64_t, struct TreeNode>& tree, vector<particle>& particles, int n_val, double threshold, double dt){
    points_checked = 0;
    for (int i = 0; i < n_val; i++){
        

        if (particles[i].mass < 0)
            continue;
        //printf("old position index %d, px = %f py = %f %f %f\n", particles[i].index, particles[i].px, particles[i].py, particles[i].vx, particles[i].vy);
        update_position(tree, &particles[i], threshold, dt);
        //printf("new position index %d, px = %f py = %f %f %f\n", particles[i].index, particles[i].px, particles[i].py, particles[i].vx, particles[i].vy);
        #ifdef DEBUG
        double ax, ay;
            if (mass[i] > -0.1){
                double fx = 0, fy = 0;
                for (int j = 0; j < n_val; j++){
                    double lfx, lfy;
                    double d = 0.0;
                    if (j == i) {lfx = 0.0; lfy = 0.0;}
                    else { 
                        //double dx = position[j].first - position[i].first;
                        //double dy = position[j].second - position[i].second;
                        d = sqrt((position[j].first - position[i].first)*(position[j].first - position[i].first) + (position[j].second - position[i].second)*(position[j].second - position[i].second));                
                        d = (d < rlimit) ? rlimit : d;
                        
                        lfx = G * mass[i] * mass[j] * (position[j].first - position[i].first)/(d*d*d);
                        lfy = G * mass[i] * mass[j] * (position[j].second - position[i].second)/(d*d*d);
                        fx += lfx;
                        fy += lfy;
                    }
                    //printf("DEBUG %d to %d force %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n", j, i, lfx, lfy, d, position[j].first - position[i].first, position[j].second - position[i].second, position[j].first, position[i].first, position[j].second, position[i].second );
                }
                debug_force[i] = (make_pair(fx,fy));
                ax = fx / mass[i];
                ay = fy / mass[i];   
                velocity[i].first += ax * dt;
                velocity[i].second += ay * dt;
                new_pos[i].first += velocity[i].first * dt + 0.5 * ax * dt * dt;
                new_pos[i].second += velocity[i].second * dt + 0.5 * ay * dt * dt;
                particle checkpoint ={0,new_pos[i].first, new_pos[i].second,0,0,0};
                new_mass[i] = check_boundary(checkpoint) ? mass[i] : -1;
            }
            
        #endif
        
    }
    #ifdef DEBUG
        for (int i = 0; i < n_val; i++){
            position[i].first = new_pos[i].first;
            position[i].second = new_pos[i].second;
            new_mass[i] = mass[i];
        }
    #endif
    printf("Points checked %d\n", points_checked);
}

void update_position(unordered_map<uint64_t, struct TreeNode>& tree, particle* p, double threshold, double dt){
    double Fx, Fy, ax = 0.0, ay = 0.0;
    tie(Fx, Fy) = get_force(tree, *p, 0, threshold);
    ax = Fx / p->mass;
    ay = Fy / p->mass;
    p->vx += ax * dt;
    p->vy += ay * dt;
    p->px += p->vx * dt + 0.5 * ax * dt * dt;
    p->py += p->vy * dt + 0.5 * ay * dt * dt;
    p->mass = check_boundary(*p) ? p->mass : -1;
#ifdef DEBUG
    force[p->index] = make_pair(Fx, Fy);
#endif
    //printf("New point %f %f %f %f %f\n", p->mass, p->px, p->py, p->vx, p->vy);
}

tuple<double, double> get_force(unordered_map<uint64_t, struct TreeNode>& tree, particle p, uint64_t key, double threshold){
    double cordX = 0.0, cordY = 0.0, d, Fx = 0.0, Fy = 0.0;
    if (tree.find(key) == tree.end()){
        return make_pair(0, 0);
    }
    TreeNode node = tree[key];
    if (node.isExternal){
        if (p.index == node.index)
        {
            Fx = 0;
            Fy = 0;
        } else {
            cordX = node.cordX;
            cordY = node.cordY;
            d = sqrt((p.px-cordX)*(p.px-cordX)+(p.py-cordY)*(p.py-cordY));
            d = (d < rlimit) ? rlimit : d;
            Fx = G * p.mass * node.mass * (cordX - p.px) / (d*d*d);
            Fy = G * p.mass * node.mass * (cordY - p.py) / (d*d*d);
            points_checked++;
        }
        //printf("%d to %d force %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n", node.index, index, Fx, Fy, d, cordX - particle[0], cordY - particle[1], cordX, particle[0], cordY, particle[1]);
        
    }else{
        cordX = node.cordX / node.mass;
        cordY = node.cordY / node.mass;
        d = sqrt((p.px-cordX)*(p.px-cordX)+(p.py-cordY)*(p.py-cordY));
        if ((node.maxX-node.minX)/d < threshold){
            Fx = G * p.mass * node.mass * (cordX - p.px) / (d*d*d);
            Fy = G * p.mass * node.mass * (cordY - p.py) / (d*d*d);
            points_checked++;
        }else{
            uint16_t layer = node.layer + 1; 
            uint32_t index_in_layer = node.index_in_layer;
            for (int i = 0; i < 4; i++){
                uint32_t child_index = index_in_layer * 4 + i;
                uint64_t child_layer = (uint64_t)layer;
                uint64_t child_key =(child_layer<<32)| child_index;
                
                if (tree.find(child_key) == tree.end())
                    continue;
                double fx, fy;
                tie(fx, fy) = get_force(tree, p, child_key, threshold);
                Fx += fx;
                Fy += fy;
            }
            
        }
    }
    return make_tuple(Fx, Fy);
}

bool check_boundary(particle p){
    if (p.px >= MIN_X && p.px <= MAX_X && p.py >= MIN_Y && p.py <= MAX_Y)
        return true;
    return false;
}