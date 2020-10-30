#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <queue>
#include <set>
#include <memory>
#include <fstream>
#include <sstream>
#include <limits>
#include <algorithm>
#include <chrono>
#include <random>
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
    struct NodePairComp {
        bool operator()(const pair<size_t, Node *> &a, const pair<size_t, Node *> &b) {
            if (a.first != b.first) return a.first > b.first;
            return a.second->id > b.second->id;
        }
    };

    struct NodeComp {
        bool operator()(const Node *a, const Node *b) {
            if (a->distance != b->distance) return a->distance < b->distance;
            return a->id < b->id;
        }
    };

    typedef priority_queue<pair<size_t, Node *>, vector<pair<size_t, Node *>>, NodePairComp> PriorityQueue;

    // we define the half of the maximum value of size_t be infinity
    // here half is used to prevent overflow
    constexpr static size_t INF = numeric_limits<size_t>::max() / 4;

    istream &in;                                        // Input stream
    ostream &out;                                       // Output stream
    mt19937 g;                                          // Random number generator

    unordered_map<size_t, unique_ptr<Node> > nodes;     // Save the unique pointer of all nodes
    vector<pair<size_t, Node *> > sortedNodes;          // Sort the nodes by id
    unordered_map<Edge, size_t, EdgeHash> edges;        // Make sure edges are unique
    vector<vector<size_t> > distances;                  // distance matrix for Floyd-Warshall algorithm
    size_t maxNodeNum = INF;                            // Set maximum node number

    // Get (or create if not exist) a node by id
    Node *getNode(size_t nodeId, bool create = false) {
        auto it = nodes.find(nodeId);
        if (it == nodes.end()) {
            if (create) {
                if (maxNodeNum > 0 && nodes.size() >= maxNodeNum) {
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

public:
    // Print the shortest path to all nodes
    void printResult() {
        for (auto &p : sortedNodes) {                   // Iterate all nodes (sorted by id)
            auto node = p.second;                       // Get the node
            out << node->id << " ";                     // Print node id
            if (node->distance < INF) {
                out << node->distance << endl;          // Print distance if it's not infinity
            } else {
                out << "inf" << endl;                   // Print inf if distance is infinity
            }
        }
    }

    void bellmanFord(size_t nodeId) {
        // Initialization
        for (auto &node : nodes) {                      // Iterate all nodes
            node.second->distance = INF;                // Initialize all distance to be infinity
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
    }

    void floydWarshallInit() {
        distances.resize(nodes.size(), vector<size_t>(nodes.size(), INF));

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
    }

    void dijkstra(size_t nodeId) {
        // Initialization
        set<Node *, NodeComp> s;                        // Use a set to save the unvisited nodes
        for (auto &node : nodes) {                      // Iterate all nodes
            node.second->distance = INF;                // Initialize all distance to be infinity
            s.emplace(node.second.get());               // Insert all nodes into the set
        }
        auto srcNode = getNode(nodeId);                 // Get the source node by its id
        s.erase(srcNode);                               // Remove the source node from the set
        srcNode->distance = 0;                          // Set the distance of first node from itself zero
        for (auto &edge : srcNode->neighbors) {         // Iterate all neighbors of the source node
            auto neighbor = edge.first;                 // One neighbor node
            s.erase(neighbor);                          // Remove the neighbor from the set
            neighbor->distance = edge.second;           // Update distance with the edge
            s.emplace(neighbor);                        // Insert the neighbor back to the set
        }
        // Main loop
        while (!s.empty()) {                            // While there are nodes not visited
            auto node = *(s.begin());                   // Get the node with minimum distance
            s.erase(s.begin());                         // Remove the node from the set
            if (s.empty() || node->distance >= INF)     // Stop if the set is empty or all left distances are infinity
                break;
            for (auto &edge : node->neighbors) {        // Iterate all neighbors of the node
                auto neighbor = edge.first;             // One neighbor node
                auto it = s.find(neighbor);             // Find the neighbor in the set
                if (it == s.end()) continue;            // Skip already visited nodes
                auto newDistance = node->distance + edge.second;
                if (newDistance < neighbor->distance) {
                    s.erase(it);                        // Remove the neighbor from the set
                    neighbor->distance = newDistance;   // Update if new distance is smaller
                    s.emplace(neighbor);                // Insert the neighbor back to the set
                }
            }
        }
    }

    void dijkstraPQ(size_t nodeId) {
        // Initialization
        for (auto &node : nodes) {                      // Iterate all nodes
            node.second->visited = false;               // Initialize all nodes not visited
            node.second->distance = INF;                // Initialize all distance to be infinity
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
                    pq.emplace(newDistance, neighbor);  // Insert the neighbor to the priority queue
                }
            }
        }
    }

    // Get a random node for benchmark
    size_t getRandomNodeId() {
        static uniform_int_distribution<int> dis(0, nodes.size());
        size_t id = dis(g);
        return sortedNodes[id].first;
    }

    explicit Graph(istream &in, ostream &out, bool directed = false, size_t maxNodeNum = INF)
            : in(in), out(out), maxNodeNum(maxNodeNum) {
        istringstream iss;
        string line;
        // Read input
        while (getline(in, line)) {
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
        // Sort the nodes by id
        sort(sortedNodes.begin(), sortedNodes.end());
        for (std::size_t i = 0; i < sortedNodes.size(); i++) {
            sortedNodes[i].second->sortedId = i;
        }
    }


};

enum class Algorithm {
    DIJKSTRA,
    DIJKSTRA_PQ,
    BELLMAN_FORD,
    FLOYD_WARSHALL,
};


struct Options {
    string input;
    string output;
    size_t nodeNum = 0;
    size_t start = 0;
    Algorithm algorithm = Algorithm::DIJKSTRA_PQ;
    bool directed = false;
    bool benchmark = false;
};


// parse command line arguments
Options parseOptions(int argc, char **argv) {
    const static char *optstring = "i:o:n:s:a:db";
    const static option long_options[] = {
            {"input",     optional_argument, nullptr, 'i'},
            {"output",    optional_argument, nullptr, 'o'},
            {"node",      optional_argument, nullptr, 'n'},
            {"start",     optional_argument, nullptr, 's'},
            {"algorithm", optional_argument, nullptr, 'a'},
            {"directed",  no_argument,       nullptr, 'd'},
            {"benchmark", no_argument,       nullptr, 'b'},
            {nullptr, 0,                     nullptr, 0}
    };
    int opt, option_index = 0;
    Options options;
    while ((opt = getopt_long(argc, argv, optstring, long_options, &option_index)) != -1) {
        switch (opt) {
            case 'i':
                options.input = optarg;
                break;
            case 'o':
                options.output = optarg;
                break;
            case 'n':
                options.nodeNum = strtoul(optarg, nullptr, 10);
                break;
            case 's':
                options.start = strtoul(optarg, nullptr, 10);
                break;
            case 'a': {
                string algorithm = optarg;
                transform(algorithm.begin(), algorithm.end(), algorithm.begin(),
                          [](unsigned char c) { return std::tolower(c); });
                if (algorithm == "dijkstra") {
                    options.algorithm = Algorithm::DIJKSTRA;
                } else if (algorithm == "dijkstra-pq") {
                    options.algorithm = Algorithm::DIJKSTRA_PQ;
                } else if (algorithm == "bellman-ford") {
                    options.algorithm = Algorithm::BELLMAN_FORD;
                } else if (algorithm == "floyd") {
                    options.algorithm = Algorithm::FLOYD_WARSHALL;
                }
                break;
            }
            case 'd':
                options.directed = true;
                break;
            case 'b':
                options.benchmark = true;
                break;
            default:
                throw runtime_error("Unrecognized option");
        }
    }
    return options;
}


int main(int argc, char **argv) {
    ifstream fin;
    ofstream fout;

    auto options = parseOptions(argc, argv);

    if (!options.input.empty()) {
        fin.open(options.input);
        if (!fin.is_open()) {
            throw runtime_error("Input file " + options.input + "can't be opened!");
        }
    }
    if (!options.output.empty()) {
        fout.open(options.output);
        if (!fout.is_open()) {
            throw runtime_error("Output file " + options.output + "can't be opened!");
        }
    }

    istream &in = fin.is_open() ? fin : cin;
    ostream &out = fout.is_open() ? fout : cout;

    Graph graph(in, out, options.directed, options.nodeNum);

    if (!options.benchmark) {
        switch (options.algorithm) {
            case Algorithm::DIJKSTRA:
                out << "algorithm: dijkstra, start node: " << options.start << endl;
                graph.dijkstra(options.start);
                break;
            case Algorithm::DIJKSTRA_PQ:
                out << "algorithm: dijkstra (priority queue), start node: " << options.start << endl;
                graph.dijkstraPQ(options.start);
                break;
            case Algorithm::BELLMAN_FORD:
                out << "algorithm: bellman-ford, start node: " << options.start << endl;
                graph.bellmanFord(options.start);
                break;
            case Algorithm::FLOYD_WARSHALL:
                out << "algorithm: floyd-warshall, start node: " << options.start << endl;
                graph.floydWarshall(options.start);
                break;
        }
        graph.printResult();
    } else {
        auto start = chrono::system_clock::now();
        for (size_t i = 0; i < options.start; i++) {
            auto id = graph.getRandomNodeId();
            switch (options.algorithm) {
                case Algorithm::DIJKSTRA:
                    graph.dijkstra(id);
                    break;
                case Algorithm::DIJKSTRA_PQ:
                    graph.dijkstraPQ(id);
                    break;
                case Algorithm::BELLMAN_FORD:
                    graph.bellmanFord(id);
                    break;
                case Algorithm::FLOYD_WARSHALL:
                    graph.floydWarshall(id);
                    break;
            }
        }
        auto end = chrono::system_clock::now();
        auto time = chrono::duration_cast<chrono::milliseconds>(end - start).count();
        out << time << endl;
    }

    return 0;
}
