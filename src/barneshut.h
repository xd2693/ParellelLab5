#ifndef _BARNESHUT_H
#define _BARNESHUT_H

#include <iostream>
#include <unordered_map>
#include <stdio.h>
#include <cmath>
#include <vector>
using namespace std;

#define MIN_X       0
#define MIN_Y       0
#define MAX_X       4
#define MAX_Y       4
#define G           0.0001
#define rlimit      0.03
#define opt         false

#ifdef DEBUG
    extern vector<pair<double, double>> force;
    extern vector<pair<double, double>> debug_force;
    extern vector<pair<double, double>> position;
    extern vector<pair<double, double>> new_pos;
    extern vector<pair<double, double>> velocity;
    extern vector<double> mass;
    extern vector<double> new_mass;
#endif

struct particle{
    int index;
    double px;
    double py;
    double mass;
    double vx;
    double vy;
    int weight;
};

struct TreeNode {
    uint16_t layer; //From 0
    uint32_t index_in_layer;
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

void tree_construct(unordered_map<uint64_t, struct TreeNode>& tree, vector<particle>& particles, int n_val);

void insert_node(unordered_map<uint64_t, struct TreeNode>& tree, vector<particle>& particles, int index, uint64_t key);

tuple<uint64_t, int> find_child(unordered_map<uint64_t, struct TreeNode>& tree, particle p, uint64_t key);


void print_node(struct TreeNode node);
void print_particle(struct particle p);

void new_position(unordered_map<uint64_t, struct TreeNode>& tree, vector<particle>& particles, int n_val, double threshold, double dt);

void update_position(unordered_map<uint64_t, struct TreeNode>& tree, particle& p, double threshold, double dt);

tuple<double, double> get_force(unordered_map<uint64_t, struct TreeNode>& tree, particle& p, uint64_t key, double threshold);

bool check_boundary(particle p);

#endif