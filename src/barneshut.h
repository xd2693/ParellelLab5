#ifndef _BARNESHUT_H
#define _BARNESHUT_H

#include <iostream>
#include <unordered_map>
#include <stdio.h>
using namespace std;

#define MAX_X 4
#define MAX_Y 4

struct TreeNode {
    int key;
    double minX;
    double minY;
    double maxX;
    double maxY;
    bool isExternal;
    double cordX;
    double cordY;
    double mass;
    int index;
};

void tree_construct(unordered_map<int, struct TreeNode>& tree, double** particles, int n_val);

void insert_node(unordered_map<int, struct TreeNode>& tree, double** particles, int index, int key);

int find_child(unordered_map<int, struct TreeNode>& tree, double** particles, int index, int key);


void print_node(struct TreeNode node);

#endif