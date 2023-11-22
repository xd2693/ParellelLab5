#include <barneshut.h>
using namespace std;
//using namespace std::tr1;

void tree_construct(unordered_map<int, struct TreeNode>& tree, double** particles, int n_val){
    //unordered_map<int, struct TreeNode> tree;
    struct TreeNode root = {0, 0, 0, MAX_X, MAX_Y, true, particles[0][0], particles[0][1], particles[0][2], 0};
    tree.insert(make_pair(0, root));
    for (int i = 1; i < n_val; i++){
        insert_node(tree, particles, i, 0);
    }
    printf("tree constructed!\n");
    for(auto const &pair : tree){
        print_node(pair.second);
    }
}

void insert_node(unordered_map<int, struct TreeNode>& tree, double** particles, int index, int key){
    
    if (particles[index][2] == -1){
        return;
    }

    if (tree[key].isExternal){
        //printf("external node %d, p_index %d\n", tree[key].key, index);
 
        tree[key].isExternal = false;
        int old_index = tree[key].index;
        tree[key].cordX = 0;
        tree[key].cordY = 0;
        tree[key].mass = 0;
        tree[key].index = 0;
        insert_node(tree, particles, old_index, key);
        insert_node(tree, particles, index, key);
        return;
    }
    if(!tree[key].isExternal){
        //printf("non external node %d, p_index %d\n", tree[key].key, index);
        //print_node(tree.at(3));
        tree[key].cordX += particles[index][0] * particles[index][2];
        tree[key].cordY += particles[index][1] * particles[index][2];
        tree[key].mass += particles[index][2];
        int child_key = find_child(tree, particles, index, key);
        //printf("child key %d ",child_key);
        if (child_key == -1){
            //printf("new child!\n");
            return;
        }else if(child_key < -1){
            printf("find child error!\n");
            return;
        }
        
        insert_node(tree, particles, index, child_key);
        return;
    }
}



int find_child(unordered_map<int, struct TreeNode>& tree, double** particles, int index, int key){
    double minX, minY, maxX, maxY;
    //printf("find child\n");
    for (int i = 1; i<=4; i++){
        int child_key = key * 4 + i;
        switch (i){
            case 1:
                minX = (tree[key].minX + tree[key].maxX)/2;
                maxX = tree[key].maxX;
                minY = (tree[key].minY + tree[key].maxY)/2;
                maxY = tree[key].maxY;
                break;
            case 2:
                minX = tree[key].minX;
                maxX = (tree[key].minX + tree[key].maxX)/2;
                minY = (tree[key].minY + tree[key].maxY)/2;
                maxY = tree[key].maxY;
                break;
            case 3:
                minX = tree[key].minX;
                maxX = (tree[key].minX + tree[key].maxX)/2;
                minY = tree[key].minY;
                maxY = (tree[key].minY + tree[key].maxY)/2;
                break;
            case 4:
                minX = (tree[key].minX + tree[key].maxX)/2;
                maxX = tree[key].maxX;
                minY = tree[key].minY;
                maxY = (tree[key].minY + tree[key].maxY)/2;
                break;
        }
        if (particles[index][0] >= minX && particles[index][0] <= maxX &&
            particles[index][1] >= minY && particles[index][1] <= maxY){
            if (tree.find(child_key) == tree.end()){
                struct TreeNode child = {child_key, minX, minY, maxX, maxY, true, particles[index][0], particles[index][1], particles[index][2], index};
                tree.insert({child_key, child});
                //print_node(tree.at(child_key));
                return -1;
            }else{
                return child_key;
            }
        }
        
    }    
    return -2;    
        
    
}
void print_node(struct TreeNode node){
    printf("key %d; minX %f; minY %f; maxX %f; maxY %f; isExtern %d; cordX %f; cordY %f; mass %f; index %d \n",
    node.key, node.minX, node.minY, node.maxX, node.maxY, node.isExternal,  node.cordX, node.cordY, node.mass, node.index);
}
