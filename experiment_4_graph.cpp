#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <set>
#include <queue>
#include <stack>
#include <limits>

using namespace std;

// 论文结构体
struct Paper {
    int id;
    string title;
    string source;
    string abstract;
    vector<string> keywords;
    string year;
    string cited;
    string doi;

    // 其他字段
    string name;
    string language;
    string type;
    string affiliations;
    string novelty;
};

// 关键词频率结构体
struct KeywordFreq {
    string keyword;
    int frequency;
};

// 快速排序，用降序
void quickSort(vector<KeywordFreq>& arr, int left, int right) {
    if (left >= right) return;

    int i = left, j = right;
    int pivot = arr[(left + right) / 2].frequency;

    while (i <= j) {
        while (arr[i].frequency > pivot) i++;
        while (arr[j].frequency < pivot) j--;

        if (i <= j) {
            swap(arr[i], arr[j]);
            i++;
            j--;
        }
    }

    if (left < j) quickSort(arr, left, j);
    if (i < right) quickSort(arr, i, right);
}

// 解析CSV字段（处理引号）
vector<string> parseCSVLine(const string& line) {
    vector<string> fields;
    string field;
    bool inQuotes = false;

    for (size_t i = 0; i < line.length(); i++) {
        char c = line[i];

        if (c == '"') {
            inQuotes = !inQuotes;
        } else if (c == ',' && !inQuotes) {
            fields.push_back(field);
            field.clear();
        } else {
            field += c;
        }
    }
    fields.push_back(field);

    return fields;
}

// 解析关键词字段 [keyword1, keyword2, ...]
vector<string> parseKeywords(const string& keywordStr) {
    vector<string> keywords;

    // 移除首尾的方括号和空格
    string content = keywordStr;
    size_t start = content.find('[');
    size_t end = content.rfind(']');

    if (start == string::npos || end == string::npos || start >= end) {
        return keywords;
    }

    content = content.substr(start + 1, end - start - 1);

    // 按逗号分割关键词
    stringstream ss(content);
    string keyword;

    while (getline(ss, keyword, ',')) {
        // 去除首尾空格
        size_t first = keyword.find_first_not_of(" \t\r\n");
        size_t last = keyword.find_last_not_of(" \t\r\n");

        if (first != string::npos && last != string::npos) {
            keyword = keyword.substr(first, last - first + 1);
            if (!keyword.empty()) {
                keywords.push_back(keyword);
            }
        }
    }

    return keywords;
}


// BFS
void BFS(const vector<vector<pair<int, int>>>& adjacencyList, int startNode, const vector<Paper>& graphPapers) {
    int n = adjacencyList.size();
    if (startNode < 0 || startNode >= n) {
        cout << "无效的起始节点！" << endl;
        return;
    }

    vector<bool> visited(n, false);
    vector<int> level(n, -1);
    vector<int> predecessor(n, -1);
    queue<int> q;

    visited[startNode] = true;
    level[startNode] = 0;
    q.push(startNode);

    cout << "\n========== BFS遍历结果 ==========" << endl;
    cout << "起始节点: " << startNode << " (论文ID: " << graphPapers[startNode].id << ")" << endl;
    cout << "\n层次遍历序列:" << endl;

    int visitCount = 0;
    while (!q.empty()) {
        int u = q.front();
        q.pop();

        cout << "节点 " << u << " [层数:" << level[u] << "] ";
        if (predecessor[u] != -1) {
            cout << "(前驱:" << predecessor[u] << ") ";
        }
        cout << "论文ID:" << graphPapers[u].id << endl;
        visitCount++;

        for (const auto& edge : adjacencyList[u]) {
            int v = edge.first;
            if (!visited[v]) {
                visited[v] = true;
                level[v] = level[u] + 1;
                predecessor[v] = u;
                q.push(v);
            }
        }
    }

    cout << "\n遍历统计: 从节点 " << startNode << " 可达 " << visitCount << " 个节点" << endl;

    // 输出层数统计
    unordered_map<int, int> levelCount;
    for (int i = 0; i < n; i++) {
        if (level[i] != -1) {
            levelCount[level[i]]++;
        }
    }
    cout << "\n各层节点数量:" << endl;
    for (int i = 0; i <= 10 && levelCount.count(i); i++) {
        cout << "第" << i << "层: " << levelCount[i] << "个节点" << endl;
    }
}

// DFS
void DFSUtil(const vector<vector<pair<int, int>>>& adjacencyList, int u, vector<bool>& visited,
             vector<int>& discoveryTime, vector<int>& finishTime, int& time, vector<int>& dfsSequence,
             const vector<Paper>& graphPapers) {
    visited[u] = true;
    discoveryTime[u] = ++time;
    dfsSequence.push_back(u);

    cout << "访问节点 " << u << " [发现时间:" << discoveryTime[u] << "] 论文ID:"
         << graphPapers[u].id << endl;

    for (const auto& edge : adjacencyList[u]) {
        int v = edge.first;
        if (!visited[v]) {
            DFSUtil(adjacencyList, v, visited, discoveryTime, finishTime, time, dfsSequence, graphPapers);
        }
    }

    finishTime[u] = ++time;
}

void DFS(const vector<vector<pair<int, int>>>& adjacencyList, int startNode, const vector<Paper>& graphPapers) {
    int n = adjacencyList.size();
    if (startNode < 0 || startNode >= n) {
        cout << "无效的起始节点！" << endl;
        return;
    }

    vector<bool> visited(n, false);
    vector<int> discoveryTime(n, -1);
    vector<int> finishTime(n, -1);
    vector<int> dfsSequence;
    int time = 0;

    cout << "\n========== DFS遍历结果 ==========" << endl;
    cout << "起始节点: " << startNode << " (论文ID: " << graphPapers[startNode].id << ")" << endl;
    cout << "\n深度优先访问序列:" << endl;

    DFSUtil(adjacencyList, startNode, visited, discoveryTime, finishTime, time, dfsSequence, graphPapers);

    // 统计连通分量
    int componentCount = 1;
    for (int i = 0; i < n; i++) {
        if (!visited[i] && !adjacencyList[i].empty()) {
            cout << "\n发现新的连通分量 #" << ++componentCount << ":" << endl;
            DFSUtil(adjacencyList, i, visited, discoveryTime, finishTime, time, dfsSequence, graphPapers);
        }
    }

    cout << "\n遍历统计: " << endl;
    cout << "- 访问节点总数: " << dfsSequence.size() << endl;
    cout << "- 连通分量数量: " << componentCount << endl;
}

// Prim算法
void Prim(const vector<vector<pair<int, int>>>& adjacencyList, int startNode, const vector<Paper>& graphPapers) {
    int n = adjacencyList.size();
    if (startNode < 0 || startNode >= n) {
        cout << "无效的起始节点！" << endl;
        return;
    }

    vector<bool> inMST(n, false);
    vector<int> key(n, numeric_limits<int>::max());
    vector<int> parent(n, -1);

    // 优先队列：pair<权重, 节点>，最小堆
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;

    key[startNode] = 0;
    pq.push({0, startNode});

    vector<pair<int, int>> mstEdges; // 存储MST的边 (u, v)
    int totalWeight = 0;

    cout << "\n========== Prim最小生成树 ==========" << endl;
    cout << "起始节点: " << startNode << " (论文ID: " << graphPapers[startNode].id << ")" << endl;
    cout << "\nMST构建过程:" << endl;

    while (!pq.empty()) {
        int u = pq.top().second;
        int weight = pq.top().first;
        pq.pop();

        if (inMST[u]) continue;  //跳过已经进入MST队列的结点

        inMST[u] = true;

        if (parent[u] != -1) {
            mstEdges.push_back({parent[u], u});
            totalWeight += weight;
            cout << "添加边: " << parent[u] << " -- " << u << " (权重:" << weight << ")" << endl;
        }

        for (const auto& edge : adjacencyList[u]) {
            int v = edge.first;
            int w = edge.second;

            if (!inMST[v] && w < key[v]) {
                key[v] = w;
                parent[v] = u;
                pq.push({w, v});
            }
        }
    }

    cout << "\nMST统计:" << endl;
    cout << "- MST边数: " << mstEdges.size() << endl;
    cout << "- MST总权重: " << totalWeight << endl;
    cout << "- 包含节点数: " << count(inMST.begin(), inMST.end(), true) << endl;

    // 保存MST到文件
    ofstream mstFile("prim_mst.txt");
    mstFile << "Prim算法最小生成树\n";
    mstFile << "起始节点: " << startNode << "\n";
    mstFile << "总权重: " << totalWeight << "\n\n";
    mstFile << "MST边列表:\n";
    for (const auto& edge : mstEdges) {
        mstFile << edge.first << " -- " << edge.second << "\n";
    }
    mstFile.close();
    cout << "\nMST已保存到 prim_mst.txt" << endl;
}

// Dijkstra算法
void Dijkstra(const vector<vector<pair<int, int>>>& adjacencyList, int startNode, const vector<Paper>& graphPapers) {
    int n = adjacencyList.size();
    if (startNode < 0 || startNode >= n) {
        cout << "无效的起始节点！" << endl;
        return;
    }

    vector<int> dist(n, numeric_limits<int>::max());
    vector<int> predecessor(n, -1);
    vector<bool> visited(n, false);

    // 优先队列：pair<距离, 节点>，最小堆
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;

    dist[startNode] = 0;
    pq.push({0, startNode});

    cout << "\n========== Dijkstra最短路径 ==========" << endl;
    cout << "起始节点: " << startNode << " (论文ID: " << graphPapers[startNode].id << ")" << endl;
    cout << "\n算法执行过程:" << endl;

    int step = 0;
    while (!pq.empty()) {
        int u = pq.top().second;
        int d = pq.top().first;
        pq.pop();

        if (visited[u]) continue;
        visited[u] = true;

        step++;
        cout << "步骤" << step << ": 处理节点 " << u << " (当前距离:" << d << ")" << endl;
        
        // 遍历所有从u出发的边
        for (const auto& edge : adjacencyList[u]) {
            int v = edge.first;
            int weight = edge.second;

            if (!visited[v] && dist[u] != numeric_limits<int>::max() && dist[u] + weight < dist[v]) {
                dist[v] = dist[u] + weight;
                pq.push({dist[v], v});
            }
        }
    }

    // 统计可达节点
    int reachableCount = 0;
    int totalDistance = 0;
    for (int i = 0; i < n; i++) {
        if (dist[i] != numeric_limits<int>::max()) {
            reachableCount++;
            totalDistance += dist[i];
        }
    }

    cout << "\n最短路径统计:" << endl;
    cout << "- 可达节点数: " << reachableCount << " / " << n << endl;
    cout << "- 平均距离: " << (reachableCount > 1 ? (double)totalDistance / (reachableCount - 1) : 0) << endl;

    
    // 保存结果到文件
    ofstream dijkstraFile("dijkstra_shortest_paths.txt");
    dijkstraFile << "Dijkstra单源最短路径\n";
    dijkstraFile << "起始节点: " << startNode << "\n\n";
    dijkstraFile << "所有可达节点的最短距离:\n";
    for (int i = 0; i < n; i++) {
        if (dist[i] != numeric_limits<int>::max()) {
            dijkstraFile << "节点 " << i << ": " << dist[i];
            if (predecessor[i] != -1) {
                dijkstraFile << " (前驱: " << predecessor[i] << ")";
            }
            dijkstraFile << "\n";
        }
    }
    dijkstraFile.close();
    cout << "\n最短路径结果已保存到 dijkstra_shortest_paths.txt" << endl;
}

int main() {
    ifstream file("E:/InnovationDataset/DeepInnovationAI/DeepDiveAI_10000.csv");
    if (!file.is_open()) return 1;

    vector<Paper> papers;
    unordered_map<string, int> keywordFreqMap;
    string line;
    getline(file, line);

    while (getline(file, line)) {
        if (line.empty()) continue;
        vector<string> fields = parseCSVLine(line);
        if (fields.size() < 14) continue;

        Paper paper;
        paper.id = stoi(fields[1]);
        paper.title = fields[2];
        paper.year = fields[6];
        paper.cited = fields[7];
        paper.keywords = parseKeywords(fields[8]);

        papers.push_back(paper);
        for (const auto& kw : paper.keywords) {
            keywordFreqMap[kw]++;
        }
    }
    file.close();

    vector<KeywordFreq> keywordFreqs;
    for (const auto& pair : keywordFreqMap) {
        keywordFreqs.push_back({pair.first, pair.second});
    }

    quickSort(keywordFreqs, 0, keywordFreqs.size() - 1);

    cout << "Top 3关键词:" << endl;
    vector<string> top3Keywords;
    for (int i = 0; i < 3; i++) {
        cout << keywordFreqs[i].keyword << " (" << keywordFreqs[i].frequency << ")" << endl;
        top3Keywords.push_back(keywordFreqs[i].keyword);
    }

    set<string> top3Set(top3Keywords.begin(), top3Keywords.end());
    vector<Paper> filteredPapers;

    for (const auto& paper : papers) {
        for (const auto& kw : paper.keywords) {
            if (top3Set.find(kw) != top3Set.end()) {
                filteredPapers.push_back(paper);
                break;
            }
        }
    }

    sort(filteredPapers.begin(), filteredPapers.end(),
         [](const Paper& a, const Paper& b) { return a.id < b.id; });

    int graphSize = min(1000, (int)filteredPapers.size());
    vector<Paper> graphPapers(filteredPapers.begin(), filteredPapers.begin() + graphSize);

    vector<vector<pair<int, int>>> adjacencyList(graphSize);

    for (int i = 0; i < graphSize; i++) {
        for (int j = i + 1; j < graphSize; j++) {
            set<string> kw1(graphPapers[i].keywords.begin(), graphPapers[i].keywords.end());
            int commonCount = 0;
            for (const auto& kw : graphPapers[j].keywords) {
                if (kw1.find(kw) != kw1.end()) commonCount++;
            }
            if (commonCount > 0) {
                adjacencyList[i].push_back({j, commonCount});
                adjacencyList[j].push_back({i, commonCount});
            }
        }
    }

    ofstream nodeFile("graph_nodes.txt");
    for (int i = 0; i < graphSize; i++) {
        nodeFile << "节点" << i << " ID:" << graphPapers[i].id << " 标题:" << graphPapers[i].title << "\n";
    }
    nodeFile.close();

    ofstream adjFile("adjacency_list.txt");
    int edgeCount = 0;
    for (int i = 0; i < graphSize; i++) {
        adjFile << i << ":";
        for (const auto& edge : adjacencyList[i]) {
            adjFile << " " << edge.first << "(" << edge.second << ")";
        }
        adjFile << "\n";
        edgeCount += adjacencyList[i].size();
    }
    adjFile.close();
    edgeCount /= 2;

    cout << "节点数:" << graphSize << " 边数:" << edgeCount << endl;

    int startNode;

    // 1. BFS
    cout << "\n[1/4] BFS" << endl;
    cout << "请输入起始节点 ID (0-" << (graphSize - 1) << "): ";
    cin >> startNode;
    if (startNode >= 0 && startNode < graphSize) {
        BFS(adjacencyList, startNode, graphPapers);
    } else {
        cout << "无效的节点 ID，跳过BFS" << endl;
    }

    // 2. DFS
    cout << "\n[2/4] DFS" << endl;
    cout << "请输入起始节点 ID (0-" << (graphSize - 1) << "): ";
    cin >> startNode;
    if (startNode >= 0 && startNode < graphSize) {
        DFS(adjacencyList, startNode, graphPapers);
    } else {
        cout << "无效的节点 ID，跳过DFS" << endl;
    }

    // 3. Prim
    cout << "\n[3/4] Prim" << endl;
    cout << "请输入起始节点 ID (0-" << (graphSize - 1) << "): ";
    cin >> startNode;
    if (startNode >= 0 && startNode < graphSize) {
        Prim(adjacencyList, startNode, graphPapers);
    } else {
        cout << "无效的节点 ID，跳过Prim" << endl;
    }

    // 4. Dijkstra
    cout << "\n[4/4] Dijkstra" << endl;
    cout << "请输入起始节点 ID (0-" << (graphSize - 1) << "): ";
    cin >> startNode;
    if (startNode >= 0 && startNode < graphSize) {
        Dijkstra(adjacencyList, startNode, graphPapers);
    } else {
        cout << "无效的节点 ID，跳过Dijkstra" << endl;
    }

    return 0;
}
