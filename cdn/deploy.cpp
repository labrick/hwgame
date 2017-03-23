#include "deploy.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <time.h>
#include <sys/timeb.h>

// #define _DEBUG
//
#undef _DEBUG

#ifdef _DEBUG
    #define PRINT printf
#else
    #define PRINT(fmt, ...)
#endif

#define USER_NODE_MAX_NUM 500
#define NETWORK_NODE_MAX_NUM 1000
#define MAX_LINK_NUM_PER_NODE 20

// 答案限制条件
#define NETWORK_PATH_MAX_NUM 50000
#define NETWORK_MAX_NUM_PER_PATH 1000

#define NOT_NETWORK_NODE_ID -1

#define MAXINT 65536

timeb startTime, curTime;

// 网络路径数量
int networkPathNum = 0;
char topo_file[NETWORK_PATH_MAX_NUM*NETWORK_MAX_NUM_PER_PATH*3]; 
char *topoFileCurPointer = topo_file;

using namespace std;

typedef struct UserNode{
    int curUserNodeID;
    int bandwidth;
    int conNetNodeID;
}UserNode;
UserNode userNode[USER_NODE_MAX_NUM];

// ===========

typedef struct Edge{
    int marked;
    int flow;
    int bandwidth;
    int costPerGB;
    int networkNodeID1;
    int networkNodeID2;
    struct Edge *edge1;
    struct Edge *edge2;
}Edge, *EdgePointer;
typedef struct NetworkNode{
    int curNetworkNodeID;
    EdgePointer nextEdge;
}NetworkNode, * NetworkNodePointer;
NetworkNode networkNode[NETWORK_NODE_MAX_NUM];

int networkNodeNum = 0;
int networkLinkNum = 0;
int userNodeNum = 0;
int costPerServer = 0;
int allCost = 0;

int tmp = 0;

#define ITERATION_NUM 2
#define ITERATION_TIME (75*1000)    // ms
#define CROSSOVER_PROBABILITY 0.3   // 交叉染色体个数
#define CROSSOVER_PROBABILITY_DOWN 0.2  // 每次交叉的基因个数 
#define VARIATION_PROBABILITY 0.3   // 变异染色体个数
#define VARIATION_PROBABILITY_DOWN 0.1  // 每次变异的基因个数
// 值越大，计算量越大，筛选的越少
// #define FILTER_COEFFICIENT 2.5
double filterCoefficient = 2.5;
int chromosomeAllNum;
int chromoKeepNum;
int geneNumPerChromo;
typedef struct Chromosome{
    bool haveCalc;      // 已经计算的就不要重复计算了
    bool *geneSeq;      // bool类型比char类型快了10ms左右
    int cost;
    double probability_up;
    double probability_down;
}Chromosome, *ChromosomePointer;
Chromosome chromosome[3*NETWORK_NODE_MAX_NUM];  // 交叉变异率和<2

void init() 
{
    // memset(userNode, 0, sizeof(struct UserNode)*USER_NODE_MAX_NUM);
    // memset(userNode, 0, sizeof(struct NetworkNode)*NETWORK_NODE_MAX_NUM);
    for (int i=0; i<NETWORK_NODE_MAX_NUM; i++){
        networkNode[i].curNetworkNodeID = i;
        networkNode[i].nextEdge = NULL;
    }
    for(int i=0; i<USER_NODE_MAX_NUM; i++){
        userNode[i].curUserNodeID = i;
    }
    strcpy(topo_file, "      \n\n");
    topoFileCurPointer = topo_file + 8;
}

// 用户需求量的大的排在前面以便优先处理
void bubbleSortUserNode()
{
    for(int i=0; i<userNodeNum; i++){
        for(int tmpi=i+1; tmpi<userNodeNum; tmpi++){
            if(userNode[i].bandwidth < userNode[tmpi].bandwidth){
                UserNode tmp = userNode[tmpi];
                userNode[tmpi] = userNode[i];
                userNode[i] = tmp;
            }
        }
    }
}

void readUserNodeInfo(char * topo[MAX_EDGE_NUM], int line_num)
{    
    ftime(&curTime);
    printf("have execute time:%ldms\n", (curTime.time-startTime.time)*1000 + (curTime.millitm - startTime.millitm));
    int i=4;
    while (strlen(topo[i++]) > 2);
    PRINT("readUserNodeInfo==================\n");
    for (int j=i; j<line_num; j++) {
        userNodeNum++;
        // PRINT("----------------line%d: %s", j, topo[j]);
        char tmpForUserNodeID[5] = {0,};       // 0~500
        char tmpForConNetNodeID[5] = {0,};     // 0~1000
        char tmpForBandwidth[5] = {0,};        // 0~100
        int tmpi=0;
        int preNum = tmpi;
        while(*((char *)topo[j]+tmpi) != ' ') {
            // PRINT("first string index:%d, char:%c\n", tmpi-preNum, *(topo[j]+tmpi));
            tmpForUserNodeID[tmpi-preNum] = *(topo[j]+tmpi);
            tmpi++;
        }
        tmpi++;
        preNum = tmpi;
        while(*((char *)topo[j]+tmpi) != ' ') {
            // PRINT("second string index:%d, char:%c\n", tmpi-preNum, *(topo[j]+tmpi));
            tmpForConNetNodeID[tmpi-preNum] = *(topo[j]+tmpi);
            tmpi++;
        }
        tmpi++;
        preNum = tmpi;
        while(*((char *)topo[j]+tmpi) != '\n') {
            // PRINT("third string index:%d, char:%c\n", tmpi-preNum, *(topo[j]+tmpi));
            tmpForBandwidth[tmpi-preNum] = *(topo[j]+tmpi);
            tmpi++;
        }
        
        // PRINT("UserNodeID is string:%s, int:%d\n", tmpForUserNodeID, atoi(tmpForUserNodeID));
        // PRINT("ConNetNodeID is string:%s, int:%d\n", tmpForConNetNodeID, atoi(tmpForConNetNodeID));
        // PRINT("Bandwidth is string:%s, int:%d\n", tmpForBandwidth, atoi(tmpForBandwidth));
         userNode[atoi(tmpForUserNodeID)].bandwidth = atoi(tmpForBandwidth);
         userNode[atoi(tmpForUserNodeID)].conNetNodeID = atoi(tmpForConNetNodeID);
    }
    // bubbleSortUserNode();
}

// 建立邻接顶点至邻接列表内
void createMLGraph(int networkNodeID, EdgePointer newEdge)
{
    EdgePointer pointer;       // 节点声明
    EdgePointer previous;      // 前一个节点

    previous = NULL;

    pointer = networkNode[networkNodeID].nextEdge;

    if (pointer == NULL) {
        networkNodeNum++;
    }

    while(pointer != NULL){
        // PRINT("search the end...\n");
        previous = pointer;

        if(pointer->networkNodeID1 == networkNodeID) {
            pointer = pointer->edge1;
        } else {
            pointer = pointer->edge2;
        }
    }
    if ( previous == NULL ){
        // PRINT("--------------networkNodeID:%d, previous==NULL\n", networkNodeID);
        networkNode[networkNodeID].nextEdge = newEdge;
    } else if ( previous->networkNodeID1 == networkNodeID) {
        previous->edge1 = newEdge;
        // PRINT("previous->edge1 = newEdge\n");
    } else {
        previous->edge2 = newEdge;
        // PRINT("previous->edge2 = newEdge\n");
    }
}

// 输出邻接列表内数据
void printNetworkNodeInfo(NetworkNodePointer networkNode)
{
    EdgePointer pointer;

    pointer = networkNode->nextEdge;
    int i = 0;
    while (pointer != NULL) {
        i++;
        PRINT("(%d, %d, %d, %d , %d)", pointer->networkNodeID1, pointer->networkNodeID2, pointer->bandwidth, pointer->costPerGB, pointer->flow);
        if ( networkNode->curNetworkNodeID == pointer->networkNodeID1) {
            pointer = pointer->edge1;
        } else {
            pointer = pointer->edge2;
        }
        // 如果在这里采用pointer访问数据，则可能此时pointer为空，所以会有segment core问题
    }
    PRINT("\n");
    tmp += i;
}

bool isDirectConnect(int networkNodeID1, int networkNodeID2, int addbandwidth, int newCost)
{    
    EdgePointer pointer, previous;
    pointer = networkNode[networkNodeID1].nextEdge;     // 这里借助了index=id
    int nextNetworkNodeID;
    while (pointer != NULL) {
        previous = pointer;
        if ( networkNodeID1 == pointer->networkNodeID1) {
            nextNetworkNodeID = pointer->networkNodeID2;
            pointer = pointer->edge1;
        } else {
            nextNetworkNodeID = pointer->networkNodeID1;
            pointer = pointer->edge2;
        }
        if(nextNetworkNodeID == networkNodeID2){
            previous->bandwidth += addbandwidth;
            if(previous->costPerGB < newCost){
                previous->costPerGB = newCost;
            }
            // printf("%d - %d is connected!\n", networkNodeID1, networkNodeID2);
            return true;
        }
    }
    return false;
}

void readNetworkNodeInfo(char * topo[MAX_EDGE_NUM], int line_num)
{
    PRINT("readNetworkNodeInfo==================\n");
    ftime(&curTime);
    printf("have execute time:%ldms\n", (curTime.time-startTime.time)*1000 + (curTime.millitm - startTime.millitm));
    for (int i=4; strlen(topo[i]) > 2; i++){
        // PRINT("----------------line%d: %s", i+1, topo[i]);
        char tmpForNetworkIDStart[5] = {0,};
        char tmpForNetworkIDEnd[5] = {0,};
        char tmpForBandwidth[5] = {0,};
        char tmpForCostPerGB[5] = {0,};

        int tmpi = 0;
        int preNum = tmpi;
        while(*((char *)topo[i]+tmpi) != ' ') {
            // PRINT("first string index:%d, char:%c\n", tmpi-preNum, *(topo[i]+tmpi));
            tmpForNetworkIDStart[tmpi-preNum] = *(topo[i]+tmpi);
            tmpi++;
        }
        tmpi++;
        preNum = tmpi;
        while(*((char *)topo[i]+tmpi) != ' ') {
            //  PRINT("second string index:%d, char:%c\n", tmpi-preNum, *(topo[i]+tmpi));
            tmpForNetworkIDEnd[tmpi-preNum] = *(topo[i]+tmpi);
            tmpi++;
        }
        tmpi++;
        preNum = tmpi;
        while(*((char *)topo[i]+tmpi) != ' ') {
            // PRINT("third string index:%d, char:%c\n", tmpi-preNum, *(topo[i]+tmpi));
            tmpForBandwidth[tmpi-preNum] = *(topo[i]+tmpi);
            tmpi++;
        }
        tmpi++;
        preNum = tmpi;
        while(*((char *)topo[i]+tmpi) != '\n') {
            // PRINT("fourth string index:%d, char:%c\n", tmpi-preNum, *(topo[i]+tmpi));
            tmpForCostPerGB[tmpi-preNum] = *(topo[i]+tmpi);
            tmpi++;
        }
        // PRINT("NetworkNodeIDStart is string:%s, int:%d\n", tmpForNetworkIDStart, atoi(tmpForNetworkIDStart));
        // PRINT("NetworkNodeIDEnd is string:%s, int:%d\n", tmpForNetworkIDEnd, atoi(tmpForNetworkIDEnd));
        // PRINT("Bandwidth is string:%s, int:%d\n", tmpForBandwidth, atoi(tmpForBandwidth));
        // PRINT("CostPerGB is string:%s, int:%d\n", tmpForCostPerGB, atoi(tmpForCostPerGB));
        // 创捷Edge结构
        if(isDirectConnect(atoi(tmpForNetworkIDStart), atoi(tmpForNetworkIDEnd), atoi(tmpForBandwidth), atoi(tmpForCostPerGB))){
            printf("%d - %d repeat edge!\n", atoi(tmpForNetworkIDStart), atoi(tmpForNetworkIDEnd));
            continue;
        }
        EdgePointer newEdge = (EdgePointer) malloc(sizeof(Edge));
        if (newEdge != NULL) {
            newEdge->networkNodeID1 = atoi(tmpForNetworkIDStart);
            newEdge->networkNodeID2 = atoi(tmpForNetworkIDEnd);
            newEdge->bandwidth = atoi(tmpForBandwidth);
            newEdge->costPerGB = atoi(tmpForCostPerGB);
            newEdge->flow = 0;      // 最初Edge上没流量;
            newEdge->edge1 = NULL;
            newEdge->edge2 = NULL;
            createMLGraph(atoi(tmpForNetworkIDStart), newEdge);
            createMLGraph(atoi(tmpForNetworkIDEnd), newEdge);
            networkLinkNum++;
        }
    }
}

#define MAX_8BIT 128  // 最大值100
int dijkstra(int networkNodeIDStart,int networkNodeIDEnd, int *preNetworkNodeID)
{
    bool isAccess[NETWORK_NODE_MAX_NUM];
    int distToStart[NETWORK_NODE_MAX_NUM];

    PRINT("==dijkstra\n");
    EdgePointer pointer = networkNode[networkNodeIDStart].nextEdge;
    EdgePointer previous = pointer;
    // 默认都没有被访问过
    memset(isAccess, 0, sizeof(bool)*NETWORK_NODE_MAX_NUM);
    // 默认所有节点到起始点的距离都为最大值
    for(int i=0; i<networkNodeNum; i++){
        distToStart[i] = MAX_8BIT;
        preNetworkNodeID[i] = MAX_8BIT;
    }

    while(pointer != NULL) {
        int nextNetworkNodeID;
        previous = pointer;
        if ( networkNodeIDStart == pointer->networkNodeID1) {
            // 当前结点相连的另外一个节点ID
            nextNetworkNodeID = pointer->networkNodeID2;
            // 下一条相连边
            pointer = pointer->edge1;
        } else {
            nextNetworkNodeID = pointer->networkNodeID1;
            pointer = pointer->edge2;
        }
        // 和start节点相连的节点和距离（成本）
        if(previous->bandwidth > previous->flow){
            PRINT("node:%d to node:%d bandwidth:%d > flow:%d\n", networkNodeIDStart, nextNetworkNodeID, previous->bandwidth, previous->flow);
            distToStart[nextNetworkNodeID] = previous->costPerGB;
            // 每一个起始点相连接点的前驱点都是总的起始节点
            // 如果bandwidth已经占完，相当于不链接
            preNetworkNodeID[nextNetworkNodeID] = networkNodeIDStart;
        } else {
            PRINT("======node:%d to node:%d bandwidth:%d <= flow:%d\n", networkNodeIDStart, nextNetworkNodeID, previous->bandwidth, previous->flow);
        }
        if(nextNetworkNodeID == networkNodeIDEnd)       // 找到就停止，速度又快了一倍
            return distToStart[networkNodeIDEnd];
        // PRINT("---the previous node of node%d is %d\n", nextNetworkNodeID, networkNodeIDStart); 
        // PRINT("the costPerGB to node:%d is %d\n", nextNetworkNodeID, previous->costPerGB); 
    }
    PRINT("start from node:%d\n", networkNodeIDStart);
    isAccess[networkNodeIDStart] = true;
    distToStart[networkNodeIDStart] = 0;

    int tmpi = 0;
    // PRINT("the dist:\n");
    // PRINT("index:");
    // for(int i=0; i<networkNodeNum; i++){
    //     PRINT("%3d\t", i);
    // }
    // PRINT("\n");
    while(tmpi<=networkNodeNum){
        int minDisToStart = MAX_8BIT;
        int minDisNetworkID = networkNodeIDStart;
        
        for(int i=0; i<networkNodeNum; i++){
            if(!isAccess[i] && minDisToStart > distToStart[i]){
                minDisNetworkID = i;
                minDisToStart = distToStart[i];
            }
        }
        // PRINT("%6d:",tmpi);
        // for(int i=0; i<networkNodeNum; i++){
        //     PRINT("%3d\t", i);
        // }
        // PRINT("\n");
        // PRINT("%6d:",tmpi);
        // for(int i=0; i<networkNodeNum; i++){
        //     PRINT("%3d\t", distToStart[i]);
        // }
        // PRINT("\n");
        // PRINT("the min dist at node:%d is %d\n", minDisNetworkID, minDisToStart); 
        // PRINT("start from node:%d\n", minDisNetworkID);
        isAccess[minDisNetworkID] = true;
        tmpi++;
        pointer = networkNode[minDisNetworkID].nextEdge;
        while(pointer != NULL){
            int nextNetworkNodeID;
            previous = pointer;
            if ( minDisNetworkID == pointer->networkNodeID1) {
                // 当前结点相连的另外一个节点ID
                nextNetworkNodeID = pointer->networkNodeID2;
                // 下一条相连边
                pointer = pointer->edge1;
            } else {
                nextNetworkNodeID = pointer->networkNodeID1;
                pointer = pointer->edge2;
            }
            // 如果这时的带宽和已用流量相当，则不用更新此距离
            PRINT("node:%d to node:%d bandwidth:%d, flow:%d\n", minDisNetworkID, nextNetworkNodeID, previous->bandwidth, previous->flow);
            if(previous->bandwidth <= previous->flow){
                // PRINT("WARNING: node:%d to node:%d bandwidth:%d, flow:%d\n", previous->networkNodeID1, previous->networkNodeID2, previous->bandwidth, previous->flow);
                continue;
            }
            // PRINT("judge node:%d\n", nextNetworkNodeID);
            // 如果当前节点距起始点距离+当前距下节点距离<下节点距起始点距离
            if(!isAccess[nextNetworkNodeID] && distToStart[minDisNetworkID] + previous->costPerGB < distToStart[nextNetworkNodeID]){
                // PRINT("W(%d)=%d > W(%d):%d+thisDist:%d=%d on %d\n", nextNetworkNodeID, distToStart[nextNetworkNodeID], minDisNetworkID, distToStart[minDisNetworkID], previous->costPerGB, distToStart[minDisNetworkID]+previous->costPerGB, minDisNetworkID);
                distToStart[nextNetworkNodeID] = distToStart[minDisNetworkID] + previous->costPerGB;
                preNetworkNodeID[nextNetworkNodeID] = minDisNetworkID;
                // PRINT("the previous node of node:%d is %d\n", nextNetworkNodeID, minDisNetworkID); 
            }
            // else if (!isAccess[nextNetworkNodeID]){
            //     PRINT("W(%d)=%d > W(%d):%d+thisDist:%d=%d on %d, no need update\n", nextNetworkNodeID, distToStart[nextNetworkNodeID], minDisNetworkID, distToStart[minDisNetworkID], previous->costPerGB, distToStart[minDisNetworkID]+previous->costPerGB, minDisNetworkID);
            // 
            // }
            if(nextNetworkNodeID == networkNodeIDEnd){       // 找到就停止，速度又快了一倍
                return distToStart[networkNodeIDEnd];
            }
        }
    }
    return distToStart[networkNodeIDEnd];
}

void getNetworkIDSeqOnMinDist(int *preNetworkNodeID,int networkNodeIDStart,  int networkNodeIDEnd, vector<int> &networkNodeIDSeq)
{
    int nextNetworkNodeID = networkNodeIDEnd;
    PRINT("==getNetworkIDSeqOnMinDist\n");
    while(nextNetworkNodeID != networkNodeIDStart){
        networkNodeIDSeq.push_back(nextNetworkNodeID);
        if(nextNetworkNodeID == 128){
            PRINT("maybe have some error\n");
        }
        nextNetworkNodeID = preNetworkNodeID[nextNetworkNodeID];
    }
    networkNodeIDSeq.push_back(networkNodeIDStart);
}

int getMinFlowOnMinDist(vector<int> networkNodeIDSeq)
{
    int networkNodeIDInSeq = networkNodeIDSeq[networkNodeIDSeq.size()-1];
    int nextNetworkNodeID, tmpi = 1, minFlow = MAX_8BIT;
    int curNetworkBandwidth = 0;
    PRINT("==getMinFlowOnMinDist\n");
    EdgePointer pointer = networkNode[networkNodeIDInSeq].nextEdge;
    while(pointer != NULL){
        curNetworkBandwidth = pointer->bandwidth - pointer->flow;
        // PRINT("curNetworkNodeID:%d\n", networkNodeIDInSeq);
        if (networkNodeIDInSeq  == pointer->networkNodeID1) {
            nextNetworkNodeID = pointer->networkNodeID2;
            pointer = pointer->edge1;
        } else {
            nextNetworkNodeID = pointer->networkNodeID1;
            pointer = pointer->edge2;
        }
        // PRINT("nextNetworkNodeID:%d\n", nextNetworkNodeID);
        if(nextNetworkNodeID == networkNodeIDSeq[networkNodeIDSeq.size()-tmpi-1]){
            // printf("the bandwidth node:%d to node:%d: %d, and minFlow = %d before updated\n", networkNodeIDInSeq, nextNetworkNodeID, curNetworkBandwidth, minFlow);
            if(minFlow > curNetworkBandwidth){
                minFlow = curNetworkBandwidth;       // 最小流
            }
            pointer = networkNode[nextNetworkNodeID].nextEdge;  // 下一个搜索起始点
            networkNodeIDInSeq = nextNetworkNodeID; // 起始点ID
            tmpi++;
            if(tmpi == (int)networkNodeIDSeq.size()){
                return minFlow;
            }
        }
    }
    return -1;
}

void updateFlow(vector<int> &networkNodeIDSeq, int minFlow)
{
    int networkNodeIDInSeq = networkNodeIDSeq[networkNodeIDSeq.size()-1];
    int nextNetworkNodeID, tmpi = 1;
    PRINT("==updateFlow\n");
    EdgePointer pointer = networkNode[networkNodeIDInSeq].nextEdge;
    EdgePointer previous = pointer;
    while(pointer != NULL){
        previous = pointer;
        if (networkNodeIDInSeq  == pointer->networkNodeID1) {
            nextNetworkNodeID = pointer->networkNodeID2;
            pointer = pointer->edge1;
        } else {
            nextNetworkNodeID = pointer->networkNodeID1;
            pointer = pointer->edge2;
        }
        // PRINT("nextNetworkNodeID:%d\n", nextNetworkNodeID);
        if(nextNetworkNodeID == networkNodeIDSeq[networkNodeIDSeq.size()-tmpi-1]){
            PRINT("the flow of node%d to node%d: %d change to %d, and bandwidth:%d\n", networkNodeIDInSeq, nextNetworkNodeID, previous->flow, minFlow+previous->flow, previous->bandwidth);
            previous->flow += minFlow;
            pointer = networkNode[nextNetworkNodeID].nextEdge;  // 下一个搜索起始点
            networkNodeIDInSeq = nextNetworkNodeID; // 起始点ID
            tmpi++;
            if(tmpi == (int)networkNodeIDSeq.size()){
                break;
            }
        }
    }

}

// return cannot arrive networknode connected to usernode
int calcFlowPath(int *serverID, int serverNum)
{
    int *tmpForPreNetworkNodeID, *preNetworkNodeID;    
    int iterationCount = 0;
    PRINT("==calcFlowPath\n");
    tmpForPreNetworkNodeID = (int *)malloc(sizeof(int)*networkNodeNum);
    preNetworkNodeID = (int *)malloc(sizeof(int)*networkNodeNum);
    for(int i=0; i<userNodeNum; i++){
        int networkNodeIDStart, networkNodeIDEnd = userNode[i].conNetNodeID;
        PRINT("======================================================to userNode[%d]:%d\n", userNode[i].curUserNodeID, networkNodeIDEnd);
        // 服务器到一个用户节点的计算
        int networkNodeProvider = serverID[0];
        int minFlow, allFlow = 0; 
        bool directConnectFlag = false;
        vector<int> networkNodeIDSeq;
        memset(tmpForPreNetworkNodeID, -1, sizeof(int)*networkNodeNum);
        memset(preNetworkNodeID, -1, sizeof(int)*networkNodeNum);

        while(allFlow < userNode[i].bandwidth){
            int minCost = MAXINT;
            int searchServerCount = 0, noPathcount = 0;;
            int charNum;
            iterationCount++;
            if(iterationCount > NETWORK_PATH_MAX_NUM){
                return MAXINT;
            }
            // 先看有没有直连
            for(int iForServerID=0; iForServerID<serverNum; iForServerID++){
                if(serverID[iForServerID] == userNode[i].conNetNodeID){
                    // printf("serverID:%d connect to userNode[%d]:%d directly\n", serverID[iForServerID], userNode[i].curUserNodeID, userNode[i].conNetNodeID);
                    directConnectFlag = true;
                    networkNodeProvider = serverID[iForServerID];
                    goto DirectConnect;
                    break;
                } else {
                    directConnectFlag = false;
                }
            }
            // 没有直连
            memset(preNetworkNodeID, -1, sizeof(int)*networkNodeNum);
            for(int iForServerID=0; iForServerID<serverNum; iForServerID++){
                int tmpForMinCost = MAXINT;
                searchServerCount++;
                // PRINT("=========serverID:%d to networkNodeID:%d\n", serverID[iForServerID], networkNodeIDEnd);
                networkNodeIDStart = serverID[iForServerID];
                memset(tmpForPreNetworkNodeID, -1, sizeof(int)*networkNodeNum);

                // 计算哪台服务器到用户点1的路径最优
                tmpForMinCost = dijkstra(networkNodeIDStart, networkNodeIDEnd, tmpForPreNetworkNodeID);
                // PRINT("serverID:%d, minCost:%d\n", serverID[iForServerID], tmpForMinCost);

                if(tmpForMinCost == 128){
                    noPathcount++;
                    // PRINT("WARNING: %d-%d have no path\n", networkNodeIDStart, networkNodeIDEnd);
                    continue;
                }
                // =====最终的结果
                if(minCost > tmpForMinCost){
                    networkNodeProvider = networkNodeIDStart;
                    minCost = tmpForMinCost;
                    memcpy(preNetworkNodeID, tmpForPreNetworkNodeID, sizeof(int)*networkNodeNum);
                    PRINT("select serverID:%d to userNode:%d, minCost:%d\n", serverID[iForServerID], networkNodeIDEnd, minCost);
                }
                break;  // 先不对比了，为了提高迭代次数
                // networkNodeIDSeq.clear();
                // getNetworkIDSeqOnMinDist(tmpForPreNetworkNodeID, networkNodeIDStart, networkNodeIDEnd, networkNodeIDSeq);
                // PRINT("previous node:\n");
                // for(int i=0; i<networkNodeNum; i++){
                //     PRINT("%3d\t", i);
                // }
                // PRINT("\n");
                // for(int i=0; i<networkNodeNum; i++){
                //     PRINT("%3d\t", preNetworkNodeID[i]);
                // }
                // PRINT("\n");
                // PRINT("server %d result: cur minCost:%d, %d-%d seq's minCost: %d and the seq:",searchServerCount, minCost, networkNodeIDStart, networkNodeIDEnd, tmpForMinCost);
                // for(int i=0; i<(int)networkNodeIDSeq.size(); i++){
                //     PRINT("%d\t", networkNodeIDSeq[networkNodeIDSeq.size()-1-i]);
                // }
                // PRINT("\n");
            }
            if(noPathcount >= serverNum){
                PRINT("ERROR: no server to the userNode[%d]:%d flow:%d  > allFlow:%d\n", i, userNode[i].conNetNodeID, userNode[i].bandwidth, allFlow);
                return userNode[i].conNetNodeID;
            }
            // PRINT("previous Node ID:\n");
            // for(int i=0; i<networkNodeNum; i++){
            //    PRINT("%3d\t", i);
            // }
            // PRINT("\n");
            // for(int i=0; i<networkNodeNum; i++){
            //    PRINT("%3d\t", preNetworkNodeID[i]);
            // }
            // PRINT("\n");
            // 获得start to end的最短路径上的ID序列
            networkNodeIDSeq.clear();
            getNetworkIDSeqOnMinDist(preNetworkNodeID, networkNodeProvider, networkNodeIDEnd, networkNodeIDSeq);
            // 获得该序列上点间最小的容量值
            minFlow = getMinFlowOnMinDist(networkNodeIDSeq);
            // 根据用户需求调整最小流
            if(allFlow + minFlow > userNode[i].bandwidth){
                PRINT("minFlow:%d + allFlow:%d > userNode[%d].bandwidth:%d, so update minFlow:%d and minCost:%d\n", minFlow, allFlow, i, userNode[i].bandwidth, userNode[i].bandwidth-allFlow, minCost);
                minFlow = userNode[i].bandwidth - allFlow;
                allFlow = userNode[i].bandwidth;
            } else {
                allFlow += minFlow;
            }
            PRINT("allCost:%d + addCost:%d = %d\n", allCost, minCost*minFlow, allCost+minCost*minFlow);
            // 根据最小流更新路径上的当前流量值
            updateFlow(networkNodeIDSeq, minFlow);

            PRINT("the %d-%d min flow: %d\n", networkNodeProvider, networkNodeIDEnd, minFlow);
            PRINT("select %d-%d seq:", networkNodeProvider, networkNodeIDEnd);
            // 将关注点写入内存缓冲区
            if((int)networkNodeIDSeq.size() > NETWORK_MAX_NUM_PER_PATH){
                return MAXINT;
            }
            for(int i=0; i<(int)networkNodeIDSeq.size(); i++){
                PRINT("%d\t",networkNodeIDSeq[networkNodeIDSeq.size()-1-i]);
                // itoa(networkNodeIDSeq[networkNodeIDSeq.size()-1-i], topoFileCurPointer, 10);
                charNum = sprintf(topoFileCurPointer, "%d ", networkNodeIDSeq[networkNodeIDSeq.size()-1-i]);
                topoFileCurPointer += charNum;
                // for(int i=0; i<10; i++){
                //     PRINT("%c", topo_file[i]);
                // }
                // PRINT("\n");
            }
            allCost += minCost*minFlow;
            PRINT("\n");
            // 如果成本太高，则忽略，时间缩短15ms左右
            if(allCost+costPerServer*serverNum > chromosome[chromosomeAllNum/2].cost){
                return MAXINT;
            }
            // 计算租用费用
DirectConnect:
            if(directConnectFlag){
                // printf("serverID:%d connect to userNode[%d]:%d directly, and the flow:%d\n", networkNodeProvider, userNode[i].curUserNodeID, userNode[i].conNetNodeID, userNode[i].bandwidth);
                minFlow = userNode[i].bandwidth;
                charNum = sprintf(topoFileCurPointer, "%d ", networkNodeProvider);
                topoFileCurPointer += charNum;

                minCost = 0;
            }
            // 加入用戶Node
            charNum = sprintf(topoFileCurPointer, "%d ", userNode[i].curUserNodeID);
            topoFileCurPointer += charNum;
            // 加入流量值
            charNum = sprintf(topoFileCurPointer, "%d", minFlow);
            topoFileCurPointer += charNum;
            *(topoFileCurPointer++) = '\n';
            PRINT("\n");
            
            // PRINT("start, end, band, cost, flow\n");
            // for (int i=0; i<NETWORK_NODE_MAX_NUM; i++) {
            //     if (networkNode[i].nextEdge == NULL) {
            //         break;
            //     }
            //     PRINT("NetworkNode %d: ", i);
            //     printNetworkNodeInfo(&networkNode[i]);
            // }
            PRINT("=======%d result:allFlow:%d, userNode[%d].bandwidth:%d\n", iterationCount, allFlow, userNode[i].curUserNodeID,  userNode[i].bandwidth);
            networkPathNum++;
            PRINT("--------------------networkPathNum:%d\n", networkPathNum);
            if(directConnectFlag)
                break;
        }
    }
    *(--topoFileCurPointer) = 0;
    return NOT_NETWORK_NODE_ID;
}
// 清除当前流，重新开始
void clearFlow(NetworkNodePointer networkNode)
{
    EdgePointer pointer;
    pointer = networkNode->nextEdge;
    while (pointer != NULL) {
        pointer->flow = 0;
        if ( networkNode->curNetworkNodeID == pointer->networkNodeID1) {
            pointer = pointer->edge1;
        } else {
            pointer = pointer->edge2;
        }
        // 如果在这里采用pointer访问数据，则可能此时pointer为空，所以会有segment core问题
    }
}

void initForRestart()
{
    for(int i=0; i<networkNodeNum; i++){
        clearFlow(&networkNode[i]);
    }
    networkPathNum = 0;
    allCost = 0;
    memset(topo_file, 0, sizeof(char)*NETWORK_NODE_MAX_NUM);
    strcpy(topo_file, "      \n\n");
    topoFileCurPointer = topo_file + 8;
}

void initialize()
{
    chromosomeAllNum = networkNodeNum + (int)networkNodeNum*CROSSOVER_PROBABILITY*2 + (int)networkNodeNum*VARIATION_PROBABILITY;
    chromoKeepNum = networkNodeNum;
    geneNumPerChromo = networkNodeNum;

    if(networkNodeNum >= 400)    filterCoefficient = 1;
    else if(300<=networkNodeNum && networkNodeNum<400) filterCoefficient = 1;
    else if(200<=networkNodeNum && networkNodeNum<300) filterCoefficient = 1;
    else if(100<=networkNodeNum && networkNodeNum<200) filterCoefficient = 1;
    else if(networkNodeNum < 100) filterCoefficient = 1;

    printf("==ga_initalize, chromosomeAllNum:%d, chromoKeepNum:%d, geneNumPerChromo:%d\n", chromosomeAllNum, chromoKeepNum, geneNumPerChromo);
    srand(time(0));
    for(int i=0; i<chromosomeAllNum; i++){
        bool *geneSeq = (bool *)malloc(sizeof(bool)*geneNumPerChromo);
        if( geneSeq != NULL){
            memset(geneSeq, 0, sizeof(bool)*(geneNumPerChromo+1));
            if( i==0 ){       // 第一个就是用户直连节点为服务器，保证有解
                for(int i=0; i<userNodeNum; i++){
                    geneSeq[userNode[i].conNetNodeID] = 1;
                }
            } else {
                // 这种方法比下面的方法快了一倍
                for(int tmpi=0; tmpi<=userNodeNum*filterCoefficient; tmpi++){
                    geneSeq[rand()%geneNumPerChromo] = rand()%2;
                }
                // int bit1Count = 0;
                // do{
                //     bit1Count = 0;
                //     for(int tmpi=0; tmpi<geneNumPerChromo; tmpi++){
                //         geneSeq[tmpi] = rand()%2;
                //         if(geneSeq[tmpi])
                //             bit1Count++;
                //         printf("%d\t", geneSeq[tmpi]);
                //     }
                //     printf("\n");
                // }while(bit1Count > userNodeNum);
            }
        }
        ChromosomePointer newChromo = (ChromosomePointer)malloc(sizeof(Chromosome));
        if(newChromo != NULL){
            newChromo->geneSeq = geneSeq;
            newChromo->cost = 0;
            newChromo->haveCalc = false;
            chromosome[i] = *newChromo;
        } else {
            PRINT("malloc failed\n");
        }
    }
    // ====================
    for(int i=0; i<chromosomeAllNum; i++){
        PRINT("chromosome:%d\t", i);
        for(int tmpi=0; tmpi<geneNumPerChromo; tmpi++){
            PRINT("%d\t", (int)chromosome[i].geneSeq[tmpi]);
        }
        PRINT("\n");
    }
}

void bubbleChromo()
{
    PRINT("==ga_bubbleChromo\n");
    for(int i=0; i<chromosomeAllNum; i++){
        for(int tmpi=i+1; tmpi<chromosomeAllNum; tmpi++){
            if(chromosome[i].cost > chromosome[tmpi].cost){
                Chromosome tmp = chromosome[tmpi];
                chromosome[tmpi] = chromosome[i];
                chromosome[i] = tmp;
            }
        }
    }
}

bool culProbability = true;
void fitness()
{
    printf("==ga_fitness\n");
    int serverID[NETWORK_NODE_MAX_NUM];
    int serverNum = 0;
    PRINT("chromosomeAllNum:%d\n", chromosomeAllNum);
    for(int i=0; i<chromosomeAllNum; i++){
        serverNum = 0;
        memset(serverID, -1, sizeof(NETWORK_NODE_MAX_NUM));
        initForRestart();
        for(int tmpi=0; tmpi<geneNumPerChromo; tmpi++){
            if(chromosome[i].geneSeq[tmpi]){
                serverID[serverNum++] = tmpi;
            }
        }
        // userNodeNum*的系数不好说是多少
        if(serverNum > userNodeNum*filterCoefficient){
            chromosome[i].cost = MAXINT-1;
        } else if(chromosome[i].haveCalc){
            PRINT("have calc:%d!\n", chromosome[i].cost);
            continue;
        }else if(NOT_NETWORK_NODE_ID == calcFlowPath(serverID, serverNum)){
            PRINT("cost:%d\n", allCost+costPerServer*serverNum);
            chromosome[i].cost = allCost+costPerServer*serverNum;
            // 这里排除相似基因，收敛速度会变慢，但更方便找全局最优
            // if(chromosome[i].cost == chromosome[0].cost){
            //     chromosome[i].cost = MAXINT-1;
                chromosome[i].haveCalc = true;
            // }
        } else {
            PRINT("have no path \n");
            chromosome[i].cost = MAXINT;
        }
    }
    culProbability = true;
    bubbleChromo();
    printf("=====================min_cost:%d\n", chromosome[0].cost);
}

int runnerGambleGetChromo()
{
    PRINT("==ga_runnerGambleGetChromo, culProbability:%d\n", culProbability);
    double costSum = 0, tmpCostSum = 0;
    if(culProbability){
        for(int i=0; i<chromoKeepNum; i++){
            costSum += chromosome[i].cost;
        }
        PRINT("costSum:%f\n", costSum);
        for(int i=0; i<chromoKeepNum; i++){
            if(i == 0){
                chromosome[i].probability_down = 0;
                chromosome[i].probability_up = (double)chromosome[i].cost / costSum;
                tmpCostSum = chromosome[i].cost;
            } else {
                chromosome[i].probability_down = chromosome[i-1].probability_up;
                tmpCostSum += chromosome[i].cost;
                chromosome[i].probability_up = (double)tmpCostSum / costSum;
            }
        }
        for(int i=0; i<chromoKeepNum; i++){
            PRINT("chromosome[%d].probability_up:%f, probability_down:%f\n", i, chromosome[i].probability_up, chromosome[i].probability_down);
        }
        culProbability = false;
    }
    double probability = (rand() % 100) / 100.0;
    if(probability == 0.0){
        PRINT("probability == 0, return index 0\n");
        return 0;
    }
    PRINT("rand probability:%f\n", probability);
    for(int i=0; i<chromoKeepNum; i++){
        if((chromosome[i].probability_down < probability) && (probability <= chromosome[i].probability_up)){
            return i;
        }
    }
    PRINT("ERROR: maybe you have an error!, cur probability:%f, distribution:\n", probability);
    for(int i=0; i<chromoKeepNum; i++){
        PRINT("chromosome[%d].probability_up:%f, probability_down:%f\n", i, chromosome[i].probability_up, chromosome[i].probability_down);
    }
    return -1;
}

void crossover()
{
    printf("==ga_crossover\n");
    ftime(&curTime);
    printf("have execute time:%ldms\n", (curTime.time-startTime.time)*1000 + (curTime.millitm - startTime.millitm));
    int chromosomeIndex1, chromosomeIndex2;
    int crossoverIndex = 0;
    int crossoverChildStartIndex = chromoKeepNum;
    for(int i=0; i<(int)chromoKeepNum*CROSSOVER_PROBABILITY; i++){
        chromosomeIndex1 = runnerGambleGetChromo();
        chromosomeIndex2 = runnerGambleGetChromo();
        PRINT("chromosomeIndex1:%d, chromosomeIndex2:%d\n", chromosomeIndex1, chromosomeIndex2);
        PRINT("index:%d, crossoverChildStartIndex:%d\n", i, crossoverChildStartIndex);
        // for(int tmpi=0; tmpi<geneNumPerChromo; tmpi++){
        //     PRINT("%d\t", (int)chromosome[crossoverChildStartIndex+i*2].geneSeq[tmpi]);
        // }
        // PRINT("\n");
        // for(int tmpi=0; tmpi<geneNumPerChromo; tmpi++){
        //     PRINT("%d\t", (int)chromosome[chromosomeIndex1].geneSeq[tmpi]);
        // }
        // PRINT("\n");
        memcpy(chromosome[crossoverChildStartIndex+i*2].geneSeq, chromosome[chromosomeIndex1].geneSeq, sizeof(bool)*geneNumPerChromo);
        memcpy(chromosome[crossoverChildStartIndex+i*2+1].geneSeq, chromosome[chromosomeIndex2].geneSeq, sizeof(bool)*geneNumPerChromo);
        chromosome[crossoverChildStartIndex+i*2].haveCalc = false;
        chromosome[crossoverChildStartIndex+i*2+1].haveCalc = false;
        for(int tmpi=0; tmpi<geneNumPerChromo*CROSSOVER_PROBABILITY_DOWN; tmpi++){
            crossoverIndex = rand() % geneNumPerChromo;
            bool tmp = chromosome[crossoverChildStartIndex+i*2].geneSeq[crossoverIndex];
            chromosome[crossoverChildStartIndex+i*2].geneSeq[crossoverIndex] = chromosome[crossoverChildStartIndex+i*2+1].geneSeq[crossoverIndex];
            chromosome[crossoverChildStartIndex+i*2+1].geneSeq[crossoverIndex] = tmp;
        }
    }
}

void mutation()
{
    printf("==ga_mutation\n");
    int mutationNum = (int)chromoKeepNum * VARIATION_PROBABILITY;   //计算变异的染色体数目
    int mutationChildStartIndex =  chromoKeepNum +  chromoKeepNum*CROSSOVER_PROBABILITY;
    for(int i = 0; i < mutationNum; i ++) {   
        int mutationChromoID = rand() % chromoKeepNum;                    //计算哪条染色体发生变异
        for(int tmpi=0; tmpi<VARIATION_PROBABILITY_DOWN; tmpi++){
            int mutationGenePlace = rand() % geneNumPerChromo;          //计算发生变异的染色体上需要变异的基因位点
            bool flag = chromosome[mutationChromoID].geneSeq[mutationGenePlace];  //对相应的基因位点进行变异
            memcpy(chromosome[mutationChildStartIndex+i].geneSeq, chromosome[mutationChromoID].geneSeq, sizeof(bool)*geneNumPerChromo);
            chromosome[mutationChildStartIndex+i].haveCalc = false;
            if(flag)
                chromosome[mutationChildStartIndex+i].geneSeq[mutationGenePlace] = 0;
            else
                chromosome[mutationChildStartIndex+i].geneSeq[mutationGenePlace] = 1;
        }
    } 
}

void ga()
{
    ftime(&curTime);
    initialize();
    // for(int i=0; i<ITERATION_NUM; i++){
    while(1){
        fitness();
        crossover();
        mutation();
        ftime(&curTime);
        if((curTime.time-startTime.time)*1000 + (curTime.millitm-startTime.millitm) > ITERATION_TIME){
            break;
        } else {
            ftime(&curTime);
            printf("have execute time:%ldms\n", (curTime.time-startTime.time)*1000 + (curTime.millitm - startTime.millitm));
        }
    }
    fitness();
    for(int i=0; i<3; i++){
        PRINT("min cost:%d\n", chromosome[i].cost);
    }
}

int getUserFlow(int userNodeID)
{
    for(int i=0; i<userNodeNum; i++){
        if(userNode[i].curUserNodeID == userNodeID){
            return userNode[i].bandwidth;
        }
    }
    return -1;
}

// 检查result是否符合结果规则
void checkResult(const char * const fileName)
{
    char *topo[NETWORK_PATH_MAX_NUM+2];
    int line_num;
    line_num = read_file(topo, NETWORK_PATH_MAX_NUM, fileName);

    printf("line num is :%d \n", line_num);
    if (line_num == 0)
    {
        printf("Please input valid topo file.\n");
        return;
    }
    if(line_num-2 != atoi(topo[0])){
        printf("path num != real num\n");
        return;
    }
    int tmpi = 0;
    int preNum = tmpi;
    vector<int> data[line_num];
    char curData[5];
    for(int i=2; i<line_num; i++){
        PRINT("new line %d\n", i);
        preNum = tmpi = 0;
        while(*((char *)topo[i]+tmpi) != '\n'){
            memset(curData, 0, sizeof(curData));
            while(*((char *)topo[i]+tmpi) != ' '){
                PRINT("%c\t", *(topo[i]+tmpi));
                curData[tmpi-preNum] = *(topo[i]+tmpi);
                tmpi++;
                if(*((char *)topo[i]+tmpi) == '\n'){
                    break;
                }
            }
            PRINT("new data %d\n", atoi(curData));
            data[i].push_back(atoi(curData));
            if(*((char *)topo[i]+tmpi) == '\n'){
                break;
            }
            tmpi++;
            preNum = tmpi;
        }
    }
    int curUserNodeID = -1, preUserNodeID = -1;
    int allFlowPerUser = data[2][data[2].size()-1];
    curUserNodeID = preUserNodeID = data[2][data[2].size()-2];
    for(int i=3; i<line_num; i++){
        curUserNodeID = data[i][data[i].size()-2];
        // printf("curUserNodeID:%d, preUserNodeID:%d\n", curUserNodeID, preUserNodeID);
        if(curUserNodeID == preUserNodeID){
            allFlowPerUser += data[i][data[i].size()-1];
        }else{
            // getUserFlow(preUserNodeID);
            if(allFlowPerUser == getUserFlow(preUserNodeID)){
                printf("check userNode %d is ok!\n", preUserNodeID);
            }else{
                printf("userNode %d flow is not correct!\n", preUserNodeID);
            }
            allFlowPerUser = data[i][data[i].size()-1];
            preUserNodeID = curUserNodeID;
        }
        if(i == line_num-1){
            if(allFlowPerUser == getUserFlow(preUserNodeID)){
                printf("last check userNode %d is ok!\n", preUserNodeID);
            }else{
                printf("userNode %d flow is not correct!\n", preUserNodeID);
            }
        }
    }
}


//你要完成的功能总入口
void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename)
{
    ftime(&startTime);
    init();

    // 读取每台服务器的成本数
    costPerServer = atoi(topo[2]);
    PRINT("costPerServer:%d\n", costPerServer);
    
    // 读取网络节点参数
    readNetworkNodeInfo(topo, line_num);
    // 读取用户节点信息
    readUserNodeInfo(topo, line_num);
    PRINT("=======================\n");
    printf("NetworkNodeNum:%d, NetworkLinkNum:%d, UserNodeNum:%d, AllCost:%d\n", networkNodeNum, networkLinkNum, userNodeNum, allCost);
    int i = 0;
    PRINT("start, end, band, cost, flow\n");
    for (i=0; i<NETWORK_NODE_MAX_NUM; i++) {
        if (networkNode[i].nextEdge == NULL) {
            break;
        }
        PRINT("NetworkNode %d: ", i);
        printNetworkNodeInfo(&networkNode[i]);
    }
    PRINT("LinkItemNum: %d, NetworkNodeNum = %d\n", tmp, i);
    //===========节点信息读入结构体完毕
    // 遗传算法
    ga();
    int serverID[NETWORK_NODE_MAX_NUM];
    int serverNum = 0;
    memset(serverID, -1, sizeof(NETWORK_NODE_MAX_NUM));
    initForRestart();

    for(int tmpi=0; tmpi<geneNumPerChromo; tmpi++){
        if(chromosome[0].geneSeq[tmpi]){
            serverID[serverNum++] = tmpi;
        }
    }
    if(NOT_NETWORK_NODE_ID == calcFlowPath(serverID, serverNum)){
        printf("\ncongratulation, have an answer!~_~\n");
        printf("PathNum:%d, RentCost:%d, ServerNum:%d, CostPerServer:%d, AllCost:%d\n", networkPathNum, allCost, serverNum, costPerServer,  allCost+costPerServer*serverNum);
        printf("and serverID:");
        for(int i=0; i<serverNum; i++){
            printf("%d\t", serverID[i]);
        }
        printf("\n");
        printf("write the networkPathNum:%d\n", networkPathNum);
        char tmp[6];
        int charNum = sprintf(tmp, "%d", networkPathNum);
        for(int i=0; i<charNum; i++){
            // PRINT("%c", tmp[i]);
            *(topo_file+i) = tmp[i];
        }
    } else {
        printf("ERROR: you are failed\n");
    }
    topoFileCurPointer = topo_file;
    
    ftime(&curTime);
    printf("END you have execute time:%ldms\n", (curTime.time-startTime.time)*1000 + (curTime.millitm - startTime.millitm));

    
	// 直接调用输出文件的方法输出到指定文件中(ps请注意格式的正确性，如果有解，第一行只有一个数据；第二行为空；第三行开始才是具体的数据，数据之间用一个空格分隔开)
	write_result((const char *)topoFileCurPointer, filename);

    // checkResult(filename);
}
