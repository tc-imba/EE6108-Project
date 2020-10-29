#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <queue>
#include <memory>
#include <fstream>
#include <sstream>
#include <limits>
#include <algorithm>
#include <getopt.h>

using namespace std;

class Graph {
private:
    struct Node {
        vector<pair<Node *, size_t> > neighbors;
        size_t id = 0;
        size_t sortedId = 0;
        size_t distance = 0;
        bool visited = false;

        Node() = default;
    };

    typedef pair<Node *, Node *> Edge;

    // Define a hash function for Edge
    struct EdgeHash {
        size_t operator()(const Edge &p) const {
            return hash<size_t>()(p.first->id) ^ hash<size_t>()(p.second->id);
        }
    };

    // Define a compare function for Node
    struct NodeComp {
        bool operator()(const pair<size_t, Node *> &a, const pair<size_t, Node *> &b) {
            if (a.first != b.first) return a.first > b.first;
            return a.second->id > b.second->id;
        }
    };

    typedef priority_queue<pair<size_t, Node *>, vector<pair<size_t, Node *>>, NodeComp> PriorityQueue;

    // we define the half of the maximum value of size_t be infinity
    // here half is used to prevent overflow
    constexpr static size_t INFINITY = numeric_limits<size_t>::max() / 4;

    unordered_map<size_t, unique_ptr<Node> > nodes;     // Save the unique pointer of all nodes
    vector<pair<size_t, Node *> > sortedNodes;          // Sort the nodes by id
    unordered_map<Edge, size_t, EdgeHash> edges;        // Make sure edges are unique
    vector<vector<size_t> > distances;                  // distance matrix for Floyd-Warshall algorithm
    size_t maxNodeNum = INFINITY;                       // Set maximum node number

    // Get (or create if not exist) a node by id
    Node *getNode(size_t nodeId, bool create = false) {
        auto it = nodes.find(nodeId);
        if (it == nodes.end()) {
            if (create) {
                if (nodes.size() >= maxNodeNum) {
                    return nullptr;                     // Limit maximum node number here
                }
                auto node = make_unique<Node>();        // Create a node
                node->id = nodeId;
                sortedNodes.emplace_back(nodeId, node.get());
                it = nodes.emplace_hint(it, nodeId, move(node));
            } else {
                throw runtime_error("Node " + to_string(nodeId) + " not found");
            }

        }
        return it->second.get();
    }

    // Create an edge between too nodes, with duplication check
    void createEdge(size_t nodeAId, size_t nodeBId, size_t distance) {
        auto nodeA = getNode(nodeAId, true);            // Get node A
        auto nodeB = getNode(nodeBId, true);            // Get node B
        if (!nodeA || !nodeB) return;                   // Limit maximum node number here
        auto edge = make_pair(nodeA, nodeB);            // Form the edge
        auto it = edges.find(edge);                     // Find the edge
        if (it != edges.end()) return;                  // Do nothing if edge already exists
        nodeA->neighbors.emplace_back(nodeB, distance); // Create node B in node A's neighbor list
        edges.emplace_hint(it, edge, distance);         // Make sure edges are unique
    }

    void printResult() {
        for (auto &p : sortedNodes) {                   // Iterate all nodes (sorted by id)
            auto node = p.second;                       // Get the node
            cout << node->id << " ";                    // Print node id
            if (node->distance < INFINITY) {
                cout << node->distance << endl;         // Print distance if it's not infinity
            } else {
                cout << "inf" << endl;                  // Print inf if distance is infinity
            }
        }
    }

public:
    void bellmanFord(size_t nodeId) {
        // Initialization
        for (auto &node : nodes) {                      // Iterate all nodes
            node.second->distance = INFINITY;           // Initialize all distance to be infinity
        }
        auto srcNode = getNode(nodeId);                 // Get the source node by its id
        srcNode->distance = 0;                          // Set the distance of first node from itself zero

        // Main loop
        for (int i = 0; i < nodes.size(); i++) {        // do |V|-1 iterations for calculation, one for error detection
            for (auto &edge : edges) {                  // Iterate all edges
                auto node = edge.first.first;           // The source node of the edge
                auto neighbor = edge.first.second;      // The destination node of the edge
                auto distance = edge.second;            // The weight of the edge
                auto newDistance = node->distance + distance;
                if (newDistance < neighbor->distance) {
                    if (i < nodes.size() - 1) {
                        neighbor->distance = newDistance;// Update if new distance is smaller
                    } else {
                        throw runtime_error("Error: negative cycle detected!");
                    }
                }
            }
        }

        printResult();
    }

    void floydWarshallInit() {
        distances.resize(nodes.size(), vector<size_t>(nodes.size(), INFINITY));

        // Initialize distances of edges
        for (auto &edge : edges) {
            auto u = edge.first.first->sortedId;
            auto v = edge.first.second->sortedId;
            distances[u][v] = edge.second;
        }
        // Initialize distances from node to itself
        for (auto &node : nodes) {
            auto v = node.second->sortedId;
            distances[v][v] = 0;
        }
        // Main loop
        for (size_t k = 0; k < nodes.size(); k++) {
            for (size_t i = 0; i < nodes.size(); i++) {
                for (size_t j = 0; j < nodes.size(); j++) {
                    auto newDistance = distances[i][k] + distances[k][j];
                    if (distances[i][j] > newDistance) {
                        distances[i][j] = newDistance;
                    }
                }
            }
        }
    }

    void floydWarshall(size_t nodeId) {
        // Initialization
        auto srcNode = getNode(nodeId);                 // Get the source node by its id
        auto u = srcNode->sortedId;                     // Get the sorted if of the source node
        if (distances.empty()) {
            floydWarshallInit();
        }

        // Main loop
        for (auto &node : nodes) {                      // Iterate all nodes
            auto v = node.second->sortedId;             // Get the sorted if of the destination node
            node.second->distance = distances[u][v];    // Set the distance by the matrix
        }

        printResult();
    }


    void dijkstra(size_t nodeId) {
        // Initialization
        for (auto &node : nodes) {                      // Iterate all nodes
            node.second->visited = false;               // Initialize all nodes not visited
            node.second->distance = INFINITY;           // Initialize all distance to be infinity
        }
        auto srcNode = getNode(nodeId);                 // Get the source node by its id
        srcNode->distance = 0;                          // Set the distance of first node from itself zero
        PriorityQueue pq;                               // Use a priority queue to save the distance of unvisited nodes
        pq.emplace(0, srcNode);                         // Insert the source node into the priority queue

        // Main loop
        while (!pq.empty()) {                           // While there are nodes not visited
            auto node = pq.top().second;                // Get the node with minimum distance
            pq.pop();                                   // Pop the node from the priority queue
            if (node->visited) continue;                // Skip already visited nodes
            node->visited = true;                       // Set the node visited
            for (auto &edge : node->neighbors) {        // Iterate all neighbors of the node
                auto neighbor = edge.first;             // One neighbor node
                if (neighbor->visited) continue;        // Skip already visited nodes
                auto newDistance = node->distance + edge.second;
                if (newDistance < neighbor->distance) {
                    neighbor->distance = newDistance;   // Update if new distance is smaller
                    pq.emplace(newDistance, neighbor);  // Insert the neighbor back to the priority queue
                }
            }
        }

        printResult();
    }

    explicit Graph(const string &filename, bool directed = false, size_t maxNodeNum = INFINITY)
            : maxNodeNum(maxNodeNum) {
        ifstream fin(filename);
        istringstream iss;
        string line;
        while (getline(fin, line)) {
            if (line.empty()) continue;
            iss.clear();
            iss.str(line);
            size_t nodeAId, nodeBId;
            iss >> nodeAId >> nodeBId;
            int distance;
            if (!(iss >> distance)) {
                distance = 1;
            }
            if (distance >= 0) {
                createEdge(nodeAId, nodeBId, distance);
                if (!directed) {
                    createEdge(nodeBId, nodeAId, distance);
                }
            }
        }
        sort(sortedNodes.begin(), sortedNodes.end());
        for (std::size_t i = 0; i < sortedNodes.size(); i++) {
            sortedNodes[i].second->sortedId = i;
        }
    }


};


int main() {
    Graph graph("facebook.weighted.txt");

//    graph.bellmanFord(0);
//    graph.dijkstra(0);

//    for (size_t i = 0; i < 256; i++) {
        graph.dijkstra(255);
//        graph.floydWarshall(i);
//    }

    return 0;
}
