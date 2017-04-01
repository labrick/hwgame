#include "deploy.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <time.h>
#include <sys/timeb.h>
#include <assert.h>
#include <math.h>
#include <bitset>
#include<queue>

// #define _DEBUG
#undef _DEBUG

#ifdef _DEBUG
    #define PRINT printf
    #define ASSERT assert
#else
    #define PRINT(fmt, ...)
    #define ASSERT(x)
#endif

#define USER_NODE_MAX_NUM 500
#define NETWORK_NODE_MAX_NUM 1000
#define MAX_LINK_NUM_PER_NODE 20

// 答案限制条件
#define NETWORK_PATH_MAX_NUM 50000
#define NETWORK_MAX_NUM_PER_PATH 1000

#define PATH_SECCUSS -1
#define PATH_FAILED 0

#define MAXINT 0x7FFFFFFF

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
    int networkNodeID1; // 一种默认的规定：ID1等于网络节点头ID
    int networkNodeID2;
    struct Edge *edge1;
    struct Edge *edge2;
}Edge, *EdgePointer;
typedef struct NetworkNode{
    int curNetworkNodeID;
    EdgePointer nextEdge;
    vector<int> adjoinPoint;//存放邻接点
}NetworkNode, * NetworkNodePointer;
NetworkNode networkNode[NETWORK_NODE_MAX_NUM];

int networkNodeNum = 0;
int networkLinkNum = 0;
int userNodeNum = 0;
int costPerServer = 0;
int allCost = 0;
int allUserNeed = 0;    // 获取所有用户需求和方便判断MCMF是否成功
int superSourcePoint = 0;
int superCollectionPoint = 0;
bool startWriteToFile = false;

int tmp = 0;

#define ITERATION_KEEP_NUM 2000      // 最小成本保持代数
#define ITERATION_TIME (85*1000)    // ms
#define CROSSOVER_PROBABILITY 0.3   // 交叉染色体个数
#define CROSSOVER_PROBABILITY_DOWN 0.6  // 每次交叉的基因个数 
#define VARIATION_PROBABILITY 0.3   // 变异染色体个数
#define VARIATION_PROBABILITY_DOWN 0.4  // 每次变异的基因个数
#define SA_IN_GA_NUM    5
// 值越大，计算量越大，筛选的越少
// #define FILTER_COEFFICIENT 2.5
double filterCoefficient = 2.5;
int chromosomeBase = 10;
int chromosomeAllNum;
int chromoKeepNum;
int geneNumPerChromo;
double crossoverProbability = CROSSOVER_PROBABILITY;
typedef struct Chromosome{
    bool haveCalc;      // 已经计算的就不要重复计算了
    bitset<NETWORK_NODE_MAX_NUM> geneSeq;      // bool类型比char类型快了10ms左右
    // char *topoFile;
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
         allUserNeed += atoi(tmpForBandwidth);
    }
    // bubbleSortUserNode();
}

// 建立邻接顶点至邻接列表内
void createMLGraph(int networkNodeID, EdgePointer newEdge)
{
    EdgePointer pointer=NULL, previous=NULL;       // 节点声明

    pointer = networkNode[networkNodeID].nextEdge;
    // if (pointer == NULL) {
    //     networkNodeNum++;
    // }

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

void delMLGraphFirstEdge(int networkNodeID1)
{
    EdgePointer previous, pointer;       // 节点声明

    previous = pointer = networkNode[networkNodeID1].nextEdge;

    if(pointer != NULL){
        // PRINT("search the end...\n");
        previous = pointer;
        if(pointer->networkNodeID1 == networkNodeID1) {
            pointer = pointer->edge1;
        } else {
            pointer = pointer->edge2;
        }
        networkNode[networkNodeID1].nextEdge = pointer;
        free(previous);
    } else PRINT("WARNING: have no edge\n");
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

bool isDirectConnect(int networkNodeID1, int networkNodeID2)
{    
    EdgePointer pointer;
    pointer = networkNode[networkNodeID1].nextEdge;     // 这里借助了index=id
    int nextNetworkNodeID;
    while (pointer != NULL) {
        nextNetworkNodeID = pointer->networkNodeID2;
        pointer = pointer->edge1;
        if(nextNetworkNodeID == networkNodeID2){
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
        
        //save the adjoin node to vector
        int tmpSaveNetworkIDStart;//save the start and the end ID networknode
        int tmpSaveNetworkIDEnd;
        tmpSaveNetworkIDStart = atoi(tmpForNetworkIDStart);
        tmpSaveNetworkIDEnd = atoi(tmpForNetworkIDEnd);
        networkNode[tmpSaveNetworkIDStart].adjoinPoint.push_back(tmpSaveNetworkIDEnd);
        networkNode[tmpSaveNetworkIDEnd].adjoinPoint.push_back(tmpSaveNetworkIDStart);

        // 创捷Edge结构
        // if(isDirectConnect(atoi(tmpForNetworkIDStart), atoi(tmpForNetworkIDEnd), atoi(tmpForBandwidth), atoi(tmpForCostPerGB))){
        //     printf("%d - %d repeat edge!\n", atoi(tmpForNetworkIDStart), atoi(tmpForNetworkIDEnd));
        //     continue;
        // }
        EdgePointer newEdge = (EdgePointer) malloc(sizeof(Edge));
        if (newEdge != NULL) {
            newEdge->networkNodeID1 = tmpSaveNetworkIDStart;
            newEdge->networkNodeID2 = tmpSaveNetworkIDEnd;
            newEdge->bandwidth = atoi(tmpForBandwidth);
            newEdge->costPerGB = atoi(tmpForCostPerGB);
            newEdge->flow = 0;      // 最初Edge上没流量;
            newEdge->edge1 = NULL;
            newEdge->edge2 = NULL;
            createMLGraph(tmpSaveNetworkIDStart, newEdge);
        }
        newEdge = (EdgePointer) malloc(sizeof(Edge));
        if (newEdge != NULL) {
            newEdge->networkNodeID1 = tmpSaveNetworkIDEnd;
            newEdge->networkNodeID2 = tmpSaveNetworkIDStart;
            newEdge->bandwidth = atoi(tmpForBandwidth);
            newEdge->costPerGB = atoi(tmpForCostPerGB);
            newEdge->flow = 0;      // 最初Edge上没流量;
            newEdge->edge1 = NULL;
            newEdge->edge2 = NULL;
            createMLGraph(tmpSaveNetworkIDEnd, newEdge);
        }
        networkLinkNum++;
    }
    // 打印出所有点邻接表的信息
    // for(int i = 0; i < networkNodeNum-2; i ++)
    // {
    //     printf("networkNode[%d] :", i);
    //     for(unsigned int j = 0; j < networkNode[i].adjoinPoint.size(); j ++)
    //         printf("%d\t", networkNode[i].adjoinPoint[j]);
    //     printf("\n");
    // }
}

struct node{  
    int networkNodeID, distToStart;  
    node(int networkNodeID,int distToStart):networkNodeID(networkNodeID),distToStart(distToStart){}  
    inline bool operator<(const node &networkNode) const{  
        return distToStart > networkNode.distToStart;
    }         
};  
// start:networkNodeNum, End:networkNodeNum+1
int dijkstra(int networkNodeIDStart,int networkNodeIDEnd, int *preNetworkNodeID)
{
    bool isAccess[NETWORK_NODE_MAX_NUM];
    int distToStart[NETWORK_NODE_MAX_NUM];
    priority_queue<node> heap;
    heap.push(node(networkNodeIDStart, 0));

    ASSERT(networkNodeIDStart>=0 && networkNodeIDStart<=NETWORK_NODE_MAX_NUM);
    ASSERT(networkNodeIDEnd>=0 && networkNodeIDEnd<=NETWORK_NODE_MAX_NUM);

    PRINT("==dijkstra, start:%d, end:%d\n", networkNodeIDStart, networkNodeIDEnd);
    EdgePointer pointer = networkNode[networkNodeIDStart].nextEdge;
    EdgePointer previous = pointer;
    // 默认都没有被访问过
    memset(isAccess, 0, sizeof(bool)*NETWORK_NODE_MAX_NUM);
    // 默认所有节点到起始点的距离都为最大值
    for(int i=0; i<networkNodeNum; i++){
        distToStart[i] = MAXINT;
        preNetworkNodeID[i] = MAXINT;
    }

    // PRINT("start from node:%d\n", networkNodeIDStart);
    distToStart[networkNodeIDStart] = 0;
    int tmpi = 0;
    // PRINT("the dist:\n");
    // PRINT("index:");
    // for(int i=0; i<networkNodeNum; i++){
    //     PRINT("%5d\t", i);
    // }
    // PRINT("\n");
    while(tmpi<=networkNodeNum){
        int minDisNetworkID = networkNodeIDStart;
        
        do{
            if((int)heap.size() == 0){
                tmpi = networkNodeNum;
                break;
            }
            node tmpNode = heap.top();
            minDisNetworkID = tmpNode.networkNodeID;
            heap.pop();
            // printf("minDisNetworkID:%d, distToStart:%d, size:%d\n", minDisNetworkID, tmpNode.distToStart, (int)heap.size());
        }while(isAccess[minDisNetworkID] && !heap.empty());
        // PRINT("%6d:",tmpi);
        // for(int i=0; i<networkNodeNum; i++){
        //     PRINT("%5d\t", i);
        // }
        // PRINT("\n");
        // PRINT("%6d:",tmpi);
        // for(int i=0; i<networkNodeNum; i++){
        //     PRINT("%5d\t", distToStart[i]);
        // }
        // PRINT("\n");
        // PRINT("%6d:",tmpi);
        // for(int i=0; i<networkNodeNum; i++){
        //     PRINT("%5d\t", isAccess[i]);
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
            // 当前结点相连的另外一个节点ID
            nextNetworkNodeID = pointer->networkNodeID2;
            // 下一条相连边
            pointer = pointer->edge1;
            // 如果这时的带宽和已用流量相当，则不用更新此距离
            // printf("node:%d to node:%d bandwidth:%d, flow:%d\n", minDisNetworkID, nextNetworkNodeID, previous->bandwidth, previous->flow);
            if(previous->bandwidth <= previous->flow){
                // printf("WARNING: node:%d to node:%d bandwidth:%d, flow:%d\n", previous->networkNodeID1, previous->networkNodeID2, previous->bandwidth, previous->flow);
                continue;
            }
            // printf("judge node:%d\n", nextNetworkNodeID);
            // 如果当前节点距起始点距离+当前距下节点距离<下节点距起始点距离
            if(!isAccess[nextNetworkNodeID] && distToStart[minDisNetworkID] + previous->costPerGB < distToStart[nextNetworkNodeID]){
                // PRINT("W(%d)=%d > W(%d):%d+thisDist:%d=%d on %d\n", nextNetworkNodeID, distToStart[nextNetworkNodeID], minDisNetworkID, distToStart[minDisNetworkID], previous->costPerGB, distToStart[minDisNetworkID]+previous->costPerGB, minDisNetworkID);
                distToStart[nextNetworkNodeID] = distToStart[minDisNetworkID] + previous->costPerGB;
                heap.push(node(nextNetworkNodeID, distToStart[nextNetworkNodeID]));
                preNetworkNodeID[nextNetworkNodeID] = minDisNetworkID;
                // printf("the previous node of node:%d is %d\n", nextNetworkNodeID, minDisNetworkID); 
            }
            // else if (!isAccess[nextNetworkNodeID]){
            // PRINT("W(%d)=%d > W(%d):%d+thisDist:%d=%d on %d, no need update\n", nextNetworkNodeID, distToStart[nextNetworkNodeID], minDisNetworkID, distToStart[minDisNetworkID], previous->costPerGB, distToStart[minDisNetworkID]+previous->costPerGB, minDisNetworkID);
            // 
            // }
            if(nextNetworkNodeID == networkNodeIDEnd){       // 找到就停止，速度又快了一倍
                if(distToStart[networkNodeIDEnd] == 0)
                    return MAXINT-1;
                return distToStart[networkNodeIDEnd];
            }
        }
    }
    if(distToStart[networkNodeIDEnd] == 0)
        return MAXINT-1;
    return distToStart[networkNodeIDEnd];
}

void getNetworkIDSeqOnMinDist(int *preNetworkNodeID,int networkNodeIDStart,  int networkNodeIDEnd, vector<int> &networkNodeIDSeq)
{
    PRINT("==getNetworkIDSeqOnMinDist\n");
    PRINT("previous Node ID:\n");
    for(int i=0; i<networkNodeNum; i++){
       PRINT("%5d\t", i);
    }
    PRINT("\n");
    for(int i=0; i<networkNodeNum; i++){
       PRINT("%5d\t", preNetworkNodeID[i]);
    }
    PRINT("\n");
    int nextNetworkNodeID = networkNodeIDEnd;
    while(nextNetworkNodeID != networkNodeIDStart){
        networkNodeIDSeq.push_back(nextNetworkNodeID); 
        nextNetworkNodeID = preNetworkNodeID[nextNetworkNodeID];
    }
    networkNodeIDSeq.push_back(networkNodeIDStart);
}

int getMinFlowOnMinDist(vector<int> networkNodeIDSeq)
{
    int networkNodeIDInSeq = networkNodeIDSeq[networkNodeIDSeq.size()-1];
    int nextNetworkNodeID, tmpi = 1, minFlow = MAXINT;
    int curNetworkBandwidth = 0;
    PRINT("==getMinFlowOnMinDist\n");
    for(int i=0; i<(int)networkNodeIDSeq.size(); i++){
        PRINT("%d, ", networkNodeIDSeq[networkNodeIDSeq.size()-1-i]);
    }
    PRINT("\n");
    ASSERT(networkNodeIDSeq.size() >= 2);
    
    EdgePointer pointer = networkNode[networkNodeIDInSeq].nextEdge;
    while(pointer != NULL){
        curNetworkBandwidth = pointer->bandwidth - pointer->flow;
        // PRINT("curNetworkNodeID:%d\n", networkNodeIDInSeq);
        nextNetworkNodeID = pointer->networkNodeID2;
        pointer = pointer->edge1;
        // PRINT("nextNetworkNodeID:%d\n", nextNetworkNodeID);
        if(nextNetworkNodeID == networkNodeIDSeq[networkNodeIDSeq.size()-tmpi-1]){
            // printf("the bandwidth node:%d to node:%d: %d, and minFlow = %d before updated\n", networkNodeIDInSeq, nextNetworkNodeID, curNetworkBandwidth, minFlow);
            if(minFlow > curNetworkBandwidth){
                PRINT("%d-%d update minFlow:%d to %d\n", networkNodeIDInSeq, nextNetworkNodeID, minFlow, curNetworkBandwidth);
                minFlow = curNetworkBandwidth;       // 最小流
            }
            pointer = networkNode[nextNetworkNodeID].nextEdge;  // 下一个搜索起始点
            networkNodeIDInSeq = nextNetworkNodeID; // 起始点ID
            tmpi++;
            if(tmpi == (int)networkNodeIDSeq.size()){
                if(minFlow <= 0 || minFlow >= MAXINT)
                    printf("ERROR: getMinFlowOnMinDist, get minFlow failed\n");
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
    ASSERT(0<minFlow && minFlow<MAXINT);
    EdgePointer pointer = networkNode[networkNodeIDInSeq].nextEdge;
    EdgePointer previous = pointer;
    while(pointer != NULL){
        previous = pointer;
        nextNetworkNodeID = pointer->networkNodeID2;
        pointer = pointer->edge1;
        // PRINT("nextNetworkNodeID:%d\n", nextNetworkNodeID);
        if(nextNetworkNodeID == networkNodeIDSeq[networkNodeIDSeq.size()-tmpi-1]){
            PRINT("the flow of node:%d to node:%d: %d change to %d, and bandwidth:%d\n", networkNodeIDInSeq, nextNetworkNodeID, previous->flow, minFlow+previous->flow, previous->bandwidth);
            if(networkNodeIDInSeq != superSourcePoint)
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

void addSuperSourcePoint(int *serverID, int serverNum)
{
    for(int i=0; i<serverNum; i++){
        EdgePointer newEdge = (EdgePointer) malloc(sizeof(Edge));
        if (newEdge != NULL) {
            PRINT("add super edge %d to %d\n", superSourcePoint, serverID[i]);
            // 一种默认的规定：ID1等于节点头ID
            newEdge->networkNodeID1 = superSourcePoint;
            newEdge->networkNodeID2 = serverID[i];
            newEdge->bandwidth = MAXINT;
            newEdge->costPerGB = 0;
            newEdge->flow = 0;      // 最初Edge上没流量;
            newEdge->edge1 = NULL;
            newEdge->edge2 = NULL;
            createMLGraph(superSourcePoint, newEdge);
            networkLinkNum++;
        }
    }
    // PRINT("start, end, band, cost, flow\n");
    // for (int i=0; i<NETWORK_NODE_MAX_NUM; i++) {
    //     if (networkNode[i].nextEdge == NULL) {
    //         break;
    //     }
    //     PRINT("NetworkNode %d: ", i);
    //     printNetworkNodeInfo(&networkNode[i]);
    // }
}

int getUserNodeID(int conNetNodeID)
{
    PRINT("==getUserFlow, conNetNodeID:%d\n", conNetNodeID);
    for(int i=0; i<userNodeNum; i++){
        if(userNode[i].conNetNodeID == conNetNodeID)
            return userNode[i].curUserNodeID;
    }
    return -1;
}

void delSuperSourcePoint(int delCount)
{
    PRINT("==delMLGraphFirstEdge, delCount:%d\n", delCount);
    for(int i=0; i<delCount; i++){
        PRINT("del %d first edge\n", superSourcePoint);
        delMLGraphFirstEdge(superSourcePoint);
    }
}

int exeCalcFlowCount = 0;
// return cannot arrive networknode connected to usernode
int calcFlowPath(int *serverID, int serverNum)
{
    timeb timeStart, timeEnd;
    ftime(&timeStart);

    int *preNetworkNodeID;    
    int allFlow = 0, minCost, minFlow;
    int charNum;
    int iterationCount = 0;
    vector<int> networkNodeIDSeq;
    addSuperSourcePoint(serverID, serverNum);
   
    preNetworkNodeID = (int *)malloc(sizeof(int)*networkNodeNum);
    while((minCost = dijkstra(superSourcePoint, superCollectionPoint, preNetworkNodeID))){
        exeCalcFlowCount++;
        ftime(&timeEnd);
        // printf("calcFlowPath in dijkstra:%ldms\n", (timeEnd.time-startTime.time)*1000 + (timeEnd.millitm - startTime.millitm));
        if(minCost == MAXINT-1)       // 表示直连
            minCost = 0;
        else if(minCost == MAXINT){   // 表示已经没有了路径
            break;
        }
        iterationCount++;
        networkNodeIDSeq.clear();
        getNetworkIDSeqOnMinDist(preNetworkNodeID, superSourcePoint, superCollectionPoint, networkNodeIDSeq);
        // 获得该序列上点间最小的容量值
        minFlow = getMinFlowOnMinDist(networkNodeIDSeq);
        allFlow += minFlow;
        PRINT("allCost:%d + addCost:%d = %d\n", allCost, minCost*minFlow, allCost+minCost*minFlow);
        // 根据最小流更新路径上的当前流量值
        updateFlow(networkNodeIDSeq, minFlow);
        
        // printf("the %d-%d min flow: %d\n", superSourcePoint, superCollectionPoint, minFlow);
        if((int)networkNodeIDSeq.size() > NETWORK_MAX_NUM_PER_PATH){
            ftime(&timeEnd);
            printf("calcFlowPath: node per path too big END:%ldms\n", (timeEnd.time-startTime.time)*1000 + (timeEnd.millitm - startTime.millitm));
            delSuperSourcePoint(serverNum);
            return MAXINT;
        }
        allCost += minCost*minFlow;
        PRINT("select %d-%d seq:", superSourcePoint, superCollectionPoint);
        // 将关注点写入内存缓冲区
        // 舍弃超级原点和超级汇点
        for(int i=1; i<(int)networkNodeIDSeq.size()-1; i++){
            PRINT("%d\t",networkNodeIDSeq[networkNodeIDSeq.size()-1-i]);
            if(startWriteToFile){
                charNum = sprintf(topoFileCurPointer, "%d ", networkNodeIDSeq[networkNodeIDSeq.size()-1-i]);
                topoFileCurPointer += charNum;
            }
        }
        PRINT("\n");
        // 如果成本太高，则忽略，时间缩短15ms左右
        // if(allCost+costPerServer*serverNum > chromosome[chromosomeAllNum/2].cost){
        //     ftime(&timeEnd);
        //     printf("calcFlowPath: cost too big END:%ldms\n", (timeEnd.time-timeStart.time)*1000 + (timeEnd.millitm - timeStart.millitm));
        //     return MAXINT;
        // }
        if(startWriteToFile){
            // 加入用戶Node
            charNum = sprintf(topoFileCurPointer, "%d ", getUserNodeID(networkNodeIDSeq[1]));
            topoFileCurPointer += charNum;
            // 加入流量值
            charNum = sprintf(topoFileCurPointer, "%d", minFlow);
            topoFileCurPointer += charNum;
            *(topoFileCurPointer++) = '\n';
            PRINT("\n");
        }
        PRINT("=======%d result:allFlow:%d, allUserNeed:%d\n", iterationCount, allFlow, allUserNeed);
        networkPathNum++;
        PRINT("--------------------networkPathNum:%d\n", networkPathNum);

        if(allFlow == allUserNeed){
            delSuperSourcePoint(serverNum);
            return PATH_SECCUSS;
        }
    }
    // printf("MCMF OK, allFlow:%d\n", allFlow);

    ftime(&timeEnd);
    // printf("calcFlowPath: OK END:%ldms, calc path count:%d, the general:%lf\n", (timeEnd.time-startTime.time)*1000 + (timeEnd.millitm - startTime.millitm), exeCalcFlowCount, (double)((timeEnd.time-timeStart.time)*1000 + (timeEnd.millitm - timeStart.millitm)) / exeCalcFlowCount);
    delSuperSourcePoint(serverNum);
    return PATH_FAILED;
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
    PRINT("==initForRestart\n");
    for(int i=0; i<networkNodeNum; i++){
        clearFlow(&networkNode[i]);
    }
    networkPathNum = 0;
    allCost = 0;
    memset(topo_file, 0, sizeof(char)*NETWORK_NODE_MAX_NUM);
    strcpy(topo_file, "      \n\n");
    topoFileCurPointer = topo_file + 8;
}
typedef struct ServerInfo{
    int serverID;
    int serverFlow;
    int weigth;
}ServerInfo, *ServerInfoPointer;

void bubbleSort(ServerInfoPointer serverInfoPointer, int memNum)
{
    for(int i=0; i<memNum; i++){
        for(int tmpi=i+1; tmpi<memNum; tmpi++){
            if(serverInfoPointer[i].serverFlow*serverInfoPointer[i].weigth < serverInfoPointer[tmpi].serverFlow*serverInfoPointer[tmpi].weigth){
                ServerInfo tmp = serverInfoPointer[tmpi];
                serverInfoPointer[tmpi] = serverInfoPointer[i];
                serverInfoPointer[i] = tmp;
            }
        }
    }
}

void getServerID(int *serverID, int serverNum)
{
    vector<int> mustServerID, mustServerAllFlow;
    ServerInfo serverInfo[NETWORK_NODE_MAX_NUM];        // 为啥这里不能是数组指针
    int count = 0;
    PRINT("=====================getServerID\n");
    for(int i=0; i<NETWORK_NODE_MAX_NUM; i++){
        if (networkNode[i].nextEdge == NULL) {
            PRINT("Traversal over\n");
            break;
        } else {
            EdgePointer pointer = networkNode[i].nextEdge;
            int tmpForAllFlow = 0;
            while( pointer != NULL) {
                // PRINT("networkNode[%d]: %d + %d = %d\n", i, tmpForAllFlow, pointer->bandwidth, tmpForAllFlow + pointer->bandwidth);
                tmpForAllFlow += pointer->bandwidth;
                if ( networkNode[i].curNetworkNodeID == pointer->networkNodeID1) {
                    pointer = pointer->edge1;
                } else {
                    pointer = pointer->edge2;
                }
            }
            // 判断必选服务器
            // PRINT("networkNode[%d] flow is: %d\n", i, tmpForAllFlow);
            for(int tmpi=0; tmpi<userNodeNum; tmpi++){
                if((networkNode[i].curNetworkNodeID == userNode[i].conNetNodeID) && \
                        (tmpForAllFlow <= userNode[i].bandwidth)){
                    mustServerID.push_back(userNode[i].conNetNodeID);
                    mustServerAllFlow.push_back(tmpForAllFlow);
                }
            }
            ServerInfoPointer newServerInfo = (ServerInfoPointer)malloc(sizeof(ServerInfo)*networkNodeNum);
            if(newServerInfo != NULL){
                newServerInfo->serverID = networkNode[i].curNetworkNodeID;
                newServerInfo->serverFlow = tmpForAllFlow;
                newServerInfo->weigth = 100;
                serverInfo[count++] = *newServerInfo;
            }
        }   
        
    }
    bubbleSort(serverInfo, networkNodeNum);
    for(int i=0; i<serverNum; i++){
        if(serverInfo[i].weigth == 100){
            serverID[i] = serverInfo[i].serverID;
            serverInfo[i].weigth *= 0.8;
        }
        for(int tmpi=i; tmpi<networkNodeNum; tmpi++){
            if((serverInfo[tmpi].weigth == 100) && 
                    isDirectConnect(i, serverInfo[tmpi].serverID)){
                serverInfo[tmpi].weigth *= 0.8;
            }
        }
        bubbleSort(serverInfo, networkNodeNum);
    }
    for(int i=0; i<networkNodeNum; i++){
        PRINT("ID:%d, flow:%d, flow:%d*weight:%d=%d\n", serverInfo[i].serverID, serverInfo[i].serverFlow, serverInfo[i].serverFlow, serverInfo[i].weigth, serverInfo[i].serverFlow*serverInfo[i].weigth);
    }
    for(int i=0; i<(int)mustServerID.size(); i++){
        PRINT("serverID:%d is mustServerID\n", mustServerID[mustServerID.size()-1]);
        serverID[serverNum-i] = mustServerID[mustServerID.size()-1];
    }
    for(int i=0; i<serverNum; i++){
        PRINT("serverID:%d\tmaxFlow:%d\n", serverInfo[i].serverID, serverInfo[i].serverFlow);
        serverID[i] = serverInfo[i].serverID;
    }
}

bool firstExeFitness = true;
void createNewAnswer(int source, int newAnswer);
void initialize()
{
    chromosomeAllNum = chromosomeBase + (int)chromosomeBase*CROSSOVER_PROBABILITY*2 + (int)chromosomeBase*VARIATION_PROBABILITY + SA_IN_GA_NUM;
    chromoKeepNum = chromosomeBase;
    geneNumPerChromo = networkNodeNum-2;    // 去除两个超级节点
    firstExeFitness = true;

    if(networkNodeNum >= 400)    filterCoefficient = 1;
    else if(300<=networkNodeNum && networkNodeNum<400) filterCoefficient = 1;
    else if(200<=networkNodeNum && networkNodeNum<300) filterCoefficient = 1;
    else if(100<=networkNodeNum && networkNodeNum<200) filterCoefficient = 1;
    else if(networkNodeNum < 100) filterCoefficient = 1;

    printf("==ga_initalize, chromosomeAllNum:%d, chromoKeepNum:%d, geneNumPerChromo:%d\n", chromosomeAllNum, chromoKeepNum, geneNumPerChromo);
    srand(time(0));
    for(int i=0; i<chromosomeAllNum; i++){
        if( i==0 ){       // 第一个就是用户直连节点为服务器，保证有解
            for(int tmpi=0; tmpi<userNodeNum; tmpi++){
                chromosome[0].geneSeq.set(userNode[tmpi].conNetNodeID);
            }
        }else if(i<chromosomeAllNum-SA_IN_GA_NUM){
            // 这种方法比下面的方法快了一倍
            for(int tmpi=0; tmpi<=userNodeNum*filterCoefficient; tmpi++){
                chromosome[i].geneSeq.set(rand()%geneNumPerChromo);
            }
        }else{
            for(int tmpi=0; tmpi<SA_IN_GA_NUM; tmpi++){
                createNewAnswer(0, chromosomeAllNum-SA_IN_GA_NUM+tmpi);
            }
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

Chromosome tmpCache[NETWORK_NODE_MAX_NUM];
void merge(int low, int middle, int high)
{
    int i = low, j = middle+1;
    for(int k=low; k<=high; k++){
        printf("low=%d,high=%d\n",low,high);
        tmpCache[k] = chromosome[k];
        printf("tmpCache[%d].cost=%d\n",k,tmpCache[k].cost);
    }
    for(int k=low; k<=high; k++){
        if(i>middle) chromosome[k] = tmpCache[j++];
        else if(j>high) chromosome[k] = tmpCache[i++];
        else if(tmpCache[j].cost < tmpCache[i].cost) chromosome[k] = tmpCache[j++];
        else    chromosome[k] = tmpCache[i++];
        printf("chromosome[%d].cost=%d\n",k,chromosome[k].cost);
    }
}

void insertionSortChromo(int low, int high)
{
    int i,j;
    for(i=low+1; i<=high; i++){
//        int key = chromosome[i].cost;
//        printf("i=%d,key=%d\n",i,key);
        Chromosome tmp = chromosome[i];
        j = i-1;
        while(j>=low && chromosome[j].cost > tmp.cost){
           // chromosome[j+1] = chromosome[j];
            Chromosome tmp = chromosome[j+1];
            chromosome[j+1] = chromosome[j];
            chromosome[j] = tmp;
            j--;
        }
        chromosome[j+1] = tmp;
    }
    //for(int k=0; k<chromosomeAllNum; k++){
    //    printf("chromosome[%d]=%d\n",k,chromosome[k].cost);
    //}
}

void mergeSort(int low, int high)
{
    if(high-low <= 0) return;
//    printf("middle=%d",middle);
    else if(high-low < 50){
        insertionSortChromo(low, high);
    } else{
        int middle = (int)(low + (high-low)/2);
        mergeSort(low,middle); //左半边排序
        mergeSort(middle+1,high);  //右半边排序
        merge(low,middle,high);     //归并结果
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

bool isHaveComeOut(vector<int> &haveComeOut, int value)
{
    for(int i=0; i<(int)haveComeOut.size(); i++){
        if(haveComeOut[i] == value)
            return true;
    }
    return false;
}

bool culProbability = true;
int globalMinCost = MAXINT;
int minCostKeepCount = 0;
int gaExeCount = 0;
// 用来判断时间是不是到了，这里面的for执行的快
// false可以继续， true可以退出了
bool fitness()
{
    printf("==ga_fitness\n");
    int serverID[NETWORK_NODE_MAX_NUM];
    int serverNum = 0;
    PRINT("chromosomeAllNum:%d\n", chromosomeAllNum);
    vector<int> haveComeOut;
    int i = chromoKeepNum;
    if(firstExeFitness){
        firstExeFitness = false;
        i = 0;
    }
    for(; i<chromosomeAllNum; i++){
        // 这里排除相似基因，收敛速度会变慢，但更方便找全局最优 (1)
        // 这里做的还不够
        // if( (i != 0) && (chromosome[i].cost == chromosome[0].cost)){
        //     // printf("===========simile gene\n");
        //     chromosome[i].cost = MAXINT-1;
        //     continue;
        // }
            
        serverNum = 0;
        memset(serverID, -1, sizeof(NETWORK_NODE_MAX_NUM));
        initForRestart();
        for(int tmpi=0; tmpi<geneNumPerChromo; tmpi++){
            if(chromosome[i].geneSeq[tmpi]){
                serverID[serverNum++] = tmpi;
            }
        }
        // userNodeNum*的系数不好说是多少
        if(serverNum > (userNodeNum+SA_IN_GA_NUM)*filterCoefficient){
            printf("serverNum:%d > userNodeNum:%d * %f = %f\n", serverNum, userNodeNum, filterCoefficient, userNodeNum*filterCoefficient);
        } else if(chromosome[i].haveCalc){
            printf("chromosome[%d], have calc:%d!\n", i, chromosome[i].cost);
            if(!isHaveComeOut(haveComeOut, chromosome[i].cost)){
                haveComeOut.push_back(chromosome[i].cost);
            }
            continue;
        }else if(PATH_SECCUSS == calcFlowPath(serverID, serverNum)){
            chromosome[i].cost = allCost+costPerServer*serverNum;
            chromosome[i].haveCalc = true;
            if(chromoKeepNum<=i && i<(chromoKeepNum+chromoKeepNum*crossoverProbability*2))
                printf("crossover, chromosome[%d].cost:%d\n", i, chromosome[i].cost);
            else if((chromoKeepNum+chromoKeepNum*crossoverProbability*2)<=i && i<(chromosomeAllNum-SA_IN_GA_NUM))
                printf("mutation,  chromosome[%d].cost:%d\n", i, chromosome[i].cost);
            else if((chromosomeAllNum-SA_IN_GA_NUM)<=i && i<chromosomeAllNum)
                printf("sa,        chromosome[%d].cost:%d\n", i, chromosome[i].cost);
            else
                printf("original,  chromosome[%d].cost:%d\n", i, chromosome[i].cost);
            if(isHaveComeOut(haveComeOut, chromosome[i].cost)){
                // printf("cost:%d have come out\n", chromosome[i].cost);
                chromosome[i].cost = MAXINT;
            } else {
                haveComeOut.push_back(chromosome[i].cost);
            }
        } else {
            PRINT("have no path \n");
            chromosome[i].cost = MAXINT;
        }
        ftime(&curTime);
        if((curTime.time-startTime.time)*1000 + (curTime.millitm-startTime.millitm) > ITERATION_TIME){
            printf("saExeCount:%d\n", gaExeCount);
            printf("time out, ready to end it!\n");
            mergeSort(0,chromosomeAllNum-1);
            // bubbleChromo();
            return true;
        }
    }
    culProbability = true;
    mergeSort(0,chromosomeAllNum-1);
    // bubbleChromo();
    
    return false;
}

int runnerGambleGetChromo()
{
    PRINT("==ga_runnerGambleGetChromo, culProbability:%d\n", culProbability);
    double costSum = 0, tmpCostSum = 0;
    if(culProbability){
        for(int i=0; i<chromoKeepNum; i++){
            costSum += 1.0/chromosome[i].cost;
        }
        // printf("costSum:%f\n", costSum);
        for(int i=0; i<chromoKeepNum; i++){
            if(i == 0){
                chromosome[i].probability_down = 0;
                chromosome[i].probability_up = (double)(1.0/chromosome[i].cost) / costSum;
                tmpCostSum = 1.0/chromosome[i].cost;
            } else {
                chromosome[i].probability_down = chromosome[i-1].probability_up;
                tmpCostSum += 1.0/chromosome[i].cost;
                chromosome[i].probability_up = (double)tmpCostSum / costSum;
            }
        }
        // for(int i=0; i<chromoKeepNum; i++){
        //     PRINT("chromosome[%d].probability_up:%f, probability_down:%f\n", i, chromosome[i].probability_up, chromosome[i].probability_down);
        // }
        culProbability = false;
    }
    double probability = (rand() % 100) / 100.0;
    PRINT("rand probability:%f\n", probability);
    if(probability == 0.0){
        PRINT("probability == 0, return index 0\n");
        return 0;
    }
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

double ga_Tdelta = 0.995;//温度的下降率
int gaIterationCount = 0;
void crossover()
{
    printf("==ga_crossover, crossoverProbability:%f\n", crossoverProbability);
    gaIterationCount++;
    int chromosomeIndex1, chromosomeIndex2;
    int crossoverIndex = 0;
    int crossoverChildStartIndex = chromoKeepNum;
    crossoverProbability = crossoverProbability * ga_Tdelta;
    if(crossoverProbability < 0.1)  crossoverProbability = 0.1;
    for(int i=0; i<(int)chromoKeepNum*crossoverProbability; i++){
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
        
        chromosome[crossoverChildStartIndex+i*2].geneSeq = chromosome[chromosomeIndex1].geneSeq;
        chromosome[crossoverChildStartIndex+i*2+1].geneSeq = chromosome[chromosomeIndex2].geneSeq;
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
    int mutationChildStartIndex =  chromoKeepNum +  chromoKeepNum*crossoverProbability*2;
    int mutationNum = chromosomeAllNum - mutationChildStartIndex;   //计算变异的染色体数目
    for(int i = 0; i < mutationNum; i ++) {   
        int mutationChromoID = rand() % chromoKeepNum;                    //计算哪条染色体发生变异
        // for(int tmpi=0; tmpi<10*VARIATION_PROBABILITY_DOWN; tmpi++){
        for(int tmpi=0; tmpi<3; tmpi++){
            int mutationGenePlace = rand() % geneNumPerChromo;          //计算发生变异的染色体上需要变异的基因位点
            chromosome[mutationChildStartIndex+i].geneSeq = chromosome[mutationChromoID].geneSeq;
            chromosome[mutationChildStartIndex+i].haveCalc = false;
            // chromosome[mutationChildStartIndex+i].geneSeq.flip(mutationGenePlace);
            chromosome[mutationChildStartIndex+i].geneSeq.flip(mutationGenePlace);
        }
    } 
}

void createNewAnswer(int source, int newAnswer);
void ga()
{
    printf("==ga\n");
    initialize();
    // for(int i=0; i<ITERATION_NUM; i++){
    while(1){
        gaExeCount++;
        ftime(&curTime);
        printf("have execute time:%ldms\n", (curTime.time-startTime.time)*1000 + (curTime.millitm - startTime.millitm));
        if(fitness())   return;
        crossover();
        mutation();
        printf("===============gaExeCount:%d, min_cost:%d, keep %d times!\n", gaExeCount, chromosome[0].cost, minCostKeepCount);
        for(int i=0; i<SA_IN_GA_NUM; i++)
            createNewAnswer(0, chromosomeAllNum-SA_IN_GA_NUM+i);
        if(globalMinCost != chromosome[0].cost){
            globalMinCost = chromosome[0].cost;
            minCostKeepCount = 0;
        }else{
            minCostKeepCount++;
            printf("minCost keep %d times, ready to end it!\n", minCostKeepCount);
            if(minCostKeepCount > ITERATION_KEEP_NUM){
                break;
            }
        }
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

void addSuperCollectionPoint()
{
    for(int i=0; i<userNodeNum; i++){
        EdgePointer newEdge = (EdgePointer) malloc(sizeof(Edge));
        if (newEdge != NULL) {
            // PRINT("add super edge %d to %d\n", userNode[i].conNetNodeID, superCollectionPoint);
            newEdge->networkNodeID1 = userNode[i].conNetNodeID;
            newEdge->networkNodeID2 = superCollectionPoint;
            newEdge->bandwidth = userNode[i].bandwidth;
            newEdge->costPerGB = 0;
            newEdge->flow = 0;      // 最初Edge上没流量;
            newEdge->edge1 = NULL;
            newEdge->edge2 = NULL;
            createMLGraph(userNode[i].conNetNodeID, newEdge);
            networkLinkNum++;
        }
    }
}

void delSuperCollectionPoint()
{
    for(int i=0; i<userNodeNum; i++){
        // PRINT("del %d first edge\n", userNode[i].conNetNodeID);
        delMLGraphFirstEdge(userNode[i].conNetNodeID);
    }
}

//模拟退火算法
double T = 20;//初始化温度
double Tmin = 1e-8;//温度的下界
int k = 100;//迭代的次数
double Tdelta = 0.99999;//温度的下降率
//这里延用遗传算法的数据结构，将每个网络节点保存为一个二进制位，初始状态随机选择一个网络节点和上一个状态比较
void initializeSa(){
    int serverID[NETWORK_NODE_MAX_NUM];
    int serverNum = 0;
    printf("进入退火算法初始化阶段\n");
    geneNumPerChromo = networkNodeNum - 2;//实际的网络节点数
    // 0:保存历史最优值
    // 创建两条一模一样的染色体，其中一条需要进行变异(从消费节点)
    // 创建两条一模一样的染色体，其中一条需要进行变异(从流量最大点)
    for(int i = 0; i < 5; i ++){
        if( i<5 ){
            // 第一个就是用户直连节点为服务器，保证有解
            for(int tmpi=0; tmpi<userNodeNum; tmpi++){
                chromosome[i].geneSeq.set(userNode[tmpi].conNetNodeID);
            }
        }
        // if( i>2 ){
        //     // 第一个就是用户直连节点为服务器，保证有解
        //     getServerID(serverID, userNodeNum/2);
        //     for(int tmpi=0; tmpi<userNodeNum/2; tmpi++){
        //         chromosome[i].geneSeq.set(serverID[tmpi]);
        //     }
        // }
    }

    serverNum = 0;
    for(int tmpi=0; tmpi<geneNumPerChromo; tmpi++){
        if(chromosome[1].geneSeq[tmpi]){
            serverID[serverNum++] = tmpi;
        }
    }
    if(PATH_SECCUSS == calcFlowPath(serverID , serverNum)){
        printf("userNode:chromosome[%d].cost:%d\n", 0, allCost+costPerServer*serverNum);
        chromosome[1].cost = allCost+costPerServer*serverNum;
        chromosome[3].cost = chromosome[1].cost;
        chromosome[0].cost = chromosome[1].cost;
    }

    // serverNum = 0;
    // for(int tmpi=0; tmpi<geneNumPerChromo; tmpi++){
    //     if(chromosome[3].geneSeq[tmpi]){
    //         serverID[serverNum++] = tmpi;
    //     }
    // }
    // initForRestart();
    // if(PATH_SECCUSS == calcFlowPath(serverID , serverNum)){
    //     printf("maxFlow:chromosome[%d].cost:%d\n", 0, allCost+costPerServer*serverNum);
    //     chromosome[3].cost = allCost+costPerServer*serverNum;
    //     if(chromosome[3].cost < chromosome[0].cost){
    //         chromosome[0].geneSeq = chromosome[3].geneSeq;
    //         chromosome[0].cost = allCost+costPerServer*serverNum;
    //     }
    // } else chromosome[3].cost = MAXINT;
} 

void createNewAnswer(int source, int newAnswer){
    //状态转移需要稳定的状态开始转移到它的邻接点，作为新的状态
    //从直连的点中选择5~10个点进行处理，这几个点中选择他们的邻接点来作为新状态
    printf("==createNewAnswer\n");
    int serverID[NETWORK_NODE_MAX_NUM];
    int serverNum = 0;
    int statusChangeNum = 1;//稳态中的需要改变状态的点个数
    int mutationGenePlace;//稳态中需要改变状态的点的位置
    int tempIndex;//临时记录选定的点的index
    chromosome[newAnswer].geneSeq = chromosome[source].geneSeq;
    for(int tmpi=0; tmpi<geneNumPerChromo; tmpi++){
        if(chromosome[newAnswer].geneSeq[tmpi]){
            serverID[serverNum++] = tmpi;
        }
    }
    int conNetNodeIDRand = 0;
    double addProbabilityUp = 0.1;
    double delProbabilityUp = 0.9;
    double moveProbabilityUp = 1;
    for(int i = 0; i < statusChangeNum; i ++){
        double randProbability = (rand() % 100) / 100.0;
        mutationGenePlace = serverID[rand() % serverNum];
        if(randProbability<addProbabilityUp){
            chromosome[newAnswer].geneSeq.set(rand()%geneNumPerChromo);
            printf("[%d] randProbability:%f, add serverID %d\n", source, randProbability, mutationGenePlace);
        }else if(addProbabilityUp<randProbability && randProbability<delProbabilityUp){
            chromosome[newAnswer].geneSeq.reset(mutationGenePlace);
            printf("[%d] randProbability:%f, del serverID %d\n", source, randProbability, mutationGenePlace);
        }else if(delProbabilityUp<randProbability && randProbability<moveProbabilityUp){
            chromosome[newAnswer].geneSeq[mutationGenePlace] = 0;
            conNetNodeIDRand = rand() % networkNode[mutationGenePlace].adjoinPoint.size();
            tempIndex = networkNode[mutationGenePlace].adjoinPoint[conNetNodeIDRand];
            chromosome[newAnswer].geneSeq[tempIndex] = 1;
            printf("[%d] move serverID from %d to %d\n", source, mutationGenePlace, tempIndex);
        }
    }
}

int saExeCount = 0;
void judgeNewAnswer(int source, int newAnswer){
    printf("==judgeNewAnswer\n");
    int serverID[NETWORK_NODE_MAX_NUM];
    int serverNum = 0;
    for(int tmpi=0; tmpi<geneNumPerChromo; tmpi++){
        if(chromosome[newAnswer].geneSeq[tmpi]){
            serverID[serverNum++] = tmpi;
        }
    }

    initForRestart();
    if(PATH_SECCUSS == calcFlowPath(serverID, serverNum)){
        // printf("chromosome[%d].cost:%d\n", 2, allCost+costPerServer*serverNum);
        chromosome[newAnswer].cost = allCost+costPerServer*serverNum;
    } else {
        // printf("have no path \n");
        chromosome[newAnswer].cost = MAXINT;
    }
    printf("chromosome[0].cost:%d, chromosome[%d].cost:%d, chromosome[%d].cost:%d\n",chromosome[0].cost, source, chromosome[source].cost, newAnswer, chromosome[newAnswer].cost);
    if(chromosome[source].cost > chromosome[newAnswer].cost){
        chromosome[source].geneSeq = chromosome[newAnswer].geneSeq;
        chromosome[source].cost = chromosome[newAnswer].cost;

        if(chromosome[source].cost < chromosome[0].cost){
            chromosome[0].geneSeq = chromosome[source].geneSeq;
            chromosome[0].cost = chromosome[source].cost;
        }
    } else if(chromosome[source].cost < chromosome[newAnswer].cost){
        int subCost = chromosome[source].cost - chromosome[newAnswer].cost;
        double accetpProbability = exp(subCost / T);
        if(accetpProbability > (rand()%100 / 100.0)){
            printf("accept new answer\n");
            chromosome[source].geneSeq = chromosome[newAnswer].geneSeq;
            chromosome[source].cost = chromosome[newAnswer].cost;
        }else printf("deny new answer\n");
    }
}

void sa(){
    initializeSa();
    //unsigned int r = 0;
    printf("origianl T:%f\n", T);
    while(T > Tmin){
        saExeCount++;
        createNewAnswer(1, 2);
        judgeNewAnswer(1, 2);

        createNewAnswer(3, 4);
        judgeNewAnswer(3, 4);
        
        ftime(&curTime);
        printf("current T:%f\n", T);
        //r ++;
        printf("===============saExeCount:%d, min_cost:%d, keep %d times!\n", saExeCount, chromosome[0].cost, minCostKeepCount);
        if((curTime.time-startTime.time)*1000 + (curTime.millitm-startTime.millitm) > ITERATION_TIME){
            printf("time out, ready to end it!\n");
            break;
        }
        if(globalMinCost != chromosome[0].cost){
            globalMinCost = chromosome[0].cost;
            minCostKeepCount = 0;
        }else{
            minCostKeepCount++;
            printf("minCost keep %d times, ready to end it!\n", minCostKeepCount);
            if(minCostKeepCount > ITERATION_KEEP_NUM){
                break;
            }
            // if(minCostKeepCount > 200)
            //     T = T * (2-Tdelta);
            // else   
            //  T = T * Tdelta;     // 最优值不变化就变化温度
        }
                T = T * Tdelta;     // 最优值不变化就变化温度
    }
    //printf("次数为：%d\n", r);
}

//你要完成的功能总入口
void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename)
{
    int serverID[NETWORK_NODE_MAX_NUM];
    int serverNum = 0;
    int i = 0;
    ftime(&startTime);
    init();

    // 读取每台服务器的成本数
    costPerServer = atoi(topo[2]);
    PRINT("costPerServer:%d\n", costPerServer);
    networkNodeNum = atoi(topo[0]);
    PRINT("networkNodeNum:%d\n", networkNodeNum);
    superSourcePoint = networkNodeNum;
    superCollectionPoint = networkNodeNum+1;
    networkNodeNum += 2;
    
    // 先读取用户节点信息并将超级汇点边首先加入邻接表，优先到达汇点
    readUserNodeInfo(topo, line_num);
    addSuperCollectionPoint();
    // 读取网络节点参数
    readNetworkNodeInfo(topo, line_num);
    // networkNodeNum+1 = superCollectionPoint
    PRINT("=======================\n");
    printf("NetworkNodeNum:%d, NetworkLinkNum:%d, UserNodeNum:%d, AllCost:%d\n", networkNodeNum, networkLinkNum, userNodeNum, allCost);
    // if(networkNodeNum > 500){
    //     memset(serverID, -1, sizeof(NETWORK_NODE_MAX_NUM));
    //     startWriteToFile = true;
    //     initForRestart();

    //     for(int tmpi=0; tmpi<userNodeNum; tmpi++){
    //         serverID[serverNum++] = userNode[tmpi].conNetNodeID;
    //     }
    //     goto UP_500;
    // }
    //
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
    startWriteToFile = true;
    initForRestart();
    memset(serverID, -1, sizeof(NETWORK_NODE_MAX_NUM));
    if(0<networkNodeNum && networkNodeNum<200){
        chromosomeBase  = 200;
        // 遗传算法
        ga();
        startWriteToFile = true;
        initForRestart();
        memset(serverID, -1, sizeof(NETWORK_NODE_MAX_NUM));
        for(int tmpi=0; tmpi<geneNumPerChromo; tmpi++){
            if(chromosome[0].geneSeq[tmpi]){
                serverID[serverNum++] = tmpi;
            }
        }
    }else if(200<networkNodeNum && networkNodeNum <500){
        chromosomeBase  = 100;
        // 遗传算法
        ga();
        startWriteToFile = true;
        initForRestart();
        memset(serverID, -1, sizeof(NETWORK_NODE_MAX_NUM));
        for(int tmpi=0; tmpi<geneNumPerChromo; tmpi++){
            if(chromosome[0].geneSeq[tmpi]){
                serverID[serverNum++] = tmpi;
            }
        }
    }else{
        ga();
        startWriteToFile = true;
        initForRestart();
        memset(serverID, -1, sizeof(NETWORK_NODE_MAX_NUM));
        for(int tmpi=0; tmpi<geneNumPerChromo; tmpi++){
            if(chromosome[0].geneSeq[tmpi]){
                serverID[serverNum++] = tmpi;
            }
        }
    }

    if(PATH_SECCUSS == calcFlowPath(serverID, serverNum)){
        printf("%d", allCost+costPerServer*serverNum);
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
    
    ftime(&curTime);
    printf("END you have execute time:%ldms\n", (curTime.time-startTime.time)*1000 + (curTime.millitm - startTime.millitm));

    *(--topoFileCurPointer) = 0;    // 去除多余的空行
    
	// 直接调用输出文件的方法输出到指定文件中(ps请注意格式的正确性，如果有解，第一行只有一个数据；第二行为空；第三行开始才是具体的数据，数据之间用一个空格分隔开)
	write_result((const char *)topo_file, filename);

    // checkResult(filename);
}
