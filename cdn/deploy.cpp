#include "deploy.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#define _DEBUG

#ifndef _DEBUG
#define printf
#endif

#define USER_NODE_MAX_NUM 500
#define NETWORK_NODE_MAX_NUM 1000
#define MAX_LINK_NUM_PER_NODE 20

#define NETWORK_PATH_MAX_NUM 50000
#define NETWORK_MAX_NUM_PER_PATH 1000

#define NOT_NETWORK_NODE_ID -1

// 网络路径数量
int networkPathNum = 0;
char topo_file_master[NETWORK_PATH_MAX_NUM*NETWORK_MAX_NUM_PER_PATH]; 
char topo_file_maxFlow[NETWORK_PATH_MAX_NUM*NETWORK_MAX_NUM_PER_PATH]; 
char *topoFileCurPointer = topo_file_master;

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
int allCost = 0;

int tmp = 0;

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
    strcpy(topo_file_master, "      \n\n");
    strcpy(topo_file_maxFlow, "      \n\n");
    topoFileCurPointer = topo_file_master + 8;
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
    int i=4;
    while (strlen(topo[i++]) > 2);
    printf("readUserNodeInfo==================\n");
    for (int j=i; j<line_num; j++) {
        userNodeNum++;
        // printf("----------------line%d: %s", j, topo[j]);
        char tmpForUserNodeID[5] = {0,};       // 0~500
        char tmpForConNetNodeID[5] = {0,};     // 0~1000
        char tmpForBandwidth[5] = {0,};        // 0~100
        int tmpi=0;
        int preNum = tmpi;
        while(*((char *)topo[j]+tmpi) != ' ') {
            // printf("first string index:%d, char:%c\n", tmpi-preNum, *(topo[j]+tmpi));
            tmpForUserNodeID[tmpi-preNum] = *(topo[j]+tmpi);
            tmpi++;
        }
        tmpi++;
        preNum = tmpi;
        while(*((char *)topo[j]+tmpi) != ' ') {
            // printf("second string index:%d, char:%c\n", tmpi-preNum, *(topo[j]+tmpi));
            tmpForConNetNodeID[tmpi-preNum] = *(topo[j]+tmpi);
            tmpi++;
        }
        tmpi++;
        preNum = tmpi;
        while(*((char *)topo[j]+tmpi) != '\n') {
            // printf("third string index:%d, char:%c\n", tmpi-preNum, *(topo[j]+tmpi));
            tmpForBandwidth[tmpi-preNum] = *(topo[j]+tmpi);
            tmpi++;
        }
        
        // printf("UserNodeID is string:%s, int:%d\n", tmpForUserNodeID, atoi(tmpForUserNodeID));
        // printf("ConNetNodeID is string:%s, int:%d\n", tmpForConNetNodeID, atoi(tmpForConNetNodeID));
        // printf("Bandwidth is string:%s, int:%d\n", tmpForBandwidth, atoi(tmpForBandwidth));
         userNode[atoi(tmpForUserNodeID)].bandwidth = atoi(tmpForBandwidth);
         userNode[atoi(tmpForUserNodeID)].conNetNodeID = atoi(tmpForConNetNodeID);
    }
    bubbleSortUserNode();
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
        // printf("search the end...\n");
        previous = pointer;

        if(pointer->networkNodeID1 == networkNodeID) {
            pointer = pointer->edge1;
        } else {
            pointer = pointer->edge2;
        }
    }
    if ( previous == NULL ){
        // printf("--------------networkNodeID:%d, previous==NULL\n", networkNodeID);
        networkNode[networkNodeID].nextEdge = newEdge;
    } else if ( previous->networkNodeID1 == networkNodeID) {
        previous->edge1 = newEdge;
        // printf("previous->edge1 = newEdge\n");
    } else {
        previous->edge2 = newEdge;
        // printf("previous->edge2 = newEdge\n");
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
        printf("(%d, %d, %d, %d , %d)", pointer->networkNodeID1, pointer->networkNodeID2, pointer->bandwidth, pointer->costPerGB, pointer->flow);
        if ( networkNode->curNetworkNodeID == pointer->networkNodeID1) {
            pointer = pointer->edge1;
        } else {
            pointer = pointer->edge2;
        }
        // 如果在这里采用pointer访问数据，则可能此时pointer为空，所以会有segment core问题
    }
    printf("\n");
    tmp += i;
}

void readNetworkNodeInfo(char * topo[MAX_EDGE_NUM], int line_num)
{
    printf("readNetworkNodeInfo==================\n");
    for (int i=4; strlen(topo[i]) > 2; i++){
        // printf("----------------line%d: %s", i+1, topo[i]);
        char tmpForNetworkIDStart[5] = {0,};
        char tmpForNetworkIDEnd[5] = {0,};
        char tmpForBandwidth[5] = {0,};
        char tmpForCostPerGB[5] = {0,};

        int tmpi = 0;
        int preNum = tmpi;
        while(*((char *)topo[i]+tmpi) != ' ') {
            // printf("first string index:%d, char:%c\n", tmpi-preNum, *(topo[i]+tmpi));
            tmpForNetworkIDStart[tmpi-preNum] = *(topo[i]+tmpi);
            tmpi++;
        }
        tmpi++;
        preNum = tmpi;
        while(*((char *)topo[i]+tmpi) != ' ') {
            //  printf("second string index:%d, char:%c\n", tmpi-preNum, *(topo[i]+tmpi));
            tmpForNetworkIDEnd[tmpi-preNum] = *(topo[i]+tmpi);
            tmpi++;
        }
        tmpi++;
        preNum = tmpi;
        while(*((char *)topo[i]+tmpi) != ' ') {
            // printf("third string index:%d, char:%c\n", tmpi-preNum, *(topo[i]+tmpi));
            tmpForBandwidth[tmpi-preNum] = *(topo[i]+tmpi);
            tmpi++;
        }
        tmpi++;
        preNum = tmpi;
        while(*((char *)topo[i]+tmpi) != '\n') {
            // printf("fourth string index:%d, char:%c\n", tmpi-preNum, *(topo[i]+tmpi));
            tmpForCostPerGB[tmpi-preNum] = *(topo[i]+tmpi);
            tmpi++;
        }
        // printf("NetworkNodeIDStart is string:%s, int:%d\n", tmpForNetworkIDStart, atoi(tmpForNetworkIDStart));
        // printf("NetworkNodeIDEnd is string:%s, int:%d\n", tmpForNetworkIDEnd, atoi(tmpForNetworkIDEnd));
        // printf("Bandwidth is string:%s, int:%d\n", tmpForBandwidth, atoi(tmpForBandwidth));
        // printf("CostPerGB is string:%s, int:%d\n", tmpForCostPerGB, atoi(tmpForCostPerGB));
        // 创捷Edge结构
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

bool isDirectConnect(int networkNodeID1, int networkNodeID2)
{    
    EdgePointer pointer;
    pointer = networkNode[networkNodeID1].nextEdge;
    int nextNetworkNodeID;
    while (pointer != NULL) {
        if ( networkNode->curNetworkNodeID == pointer->networkNodeID1) {
            nextNetworkNodeID = pointer->networkNodeID2;
            pointer = pointer->edge1;
        } else {
            nextNetworkNodeID = pointer->networkNodeID1;
            pointer = pointer->edge2;
        }
        if(nextNetworkNodeID == networkNodeID2){
            printf("%d - %d is connected!\n", networkNodeID1, networkNodeID2);
            return true;
        }
    }
    return false;
}

void getServerID(int *serverID, int serverNum)
{
    vector<int> mustServerID, mustServerAllFlow;
    ServerInfo serverInfo[NETWORK_NODE_MAX_NUM];        // 为啥这里不能是数组指针
    int count = 0;
    printf("=====================getServerID\n");
    for(int i=0; i<NETWORK_NODE_MAX_NUM; i++){
        if (networkNode[i].nextEdge == NULL) {
            printf("Traversal over\n");
            break;
        } else {
            EdgePointer pointer = networkNode[i].nextEdge;
            int tmpForAllFlow = 0;
            while( pointer != NULL) {
                // printf("networkNode[%d]: %d + %d = %d\n", i, tmpForAllFlow, pointer->bandwidth, tmpForAllFlow + pointer->bandwidth);
                tmpForAllFlow += pointer->bandwidth;
                if ( networkNode[i].curNetworkNodeID == pointer->networkNodeID1) {
                    pointer = pointer->edge1;
                } else {
                    pointer = pointer->edge2;
                }
            }
            // 判断必选服务器
            // printf("networkNode[%d] flow is: %d\n", i, tmpForAllFlow);
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
        printf("ID:%d, flow:%d, flow:%d*weight:%d=%d\n", serverInfo[i].serverID, serverInfo[i].serverFlow, serverInfo[i].serverFlow, serverInfo[i].weigth, serverInfo[i].serverFlow*serverInfo[i].weigth);
    }
    for(int i=0; i<(int)mustServerID.size(); i++){
        printf("serverID:%d is mustServerID\n", mustServerID[mustServerID.size()-1]);
        serverID[serverNum-i] = mustServerID[mustServerID.size()-1];
    }
    for(int i=0; i<serverNum; i++){
        printf("serverID:%d\tmaxFlow:%d\n", serverInfo[i].serverID, serverInfo[i].serverFlow);
        serverID[i] = serverInfo[i].serverID;
    }
}

#define MAXINT 128  // 最大值100
int dijkstra(int networkNodeIDStart,int networkNodeIDEnd, int *preNetworkNodeID)
{
    bool isAccess[NETWORK_NODE_MAX_NUM];
    int distToStart[NETWORK_NODE_MAX_NUM];

    printf("==dijkstra\n");
    EdgePointer pointer = networkNode[networkNodeIDStart].nextEdge;
    EdgePointer previous = pointer;
    // 默认都没有被访问过
    memset(isAccess, 0, sizeof(bool)*NETWORK_NODE_MAX_NUM);
    // 默认所有节点到起始点的距离都为最大值
    for(int i=0; i<networkNodeNum; i++){
        distToStart[i] = MAXINT;
        preNetworkNodeID[i] = MAXINT;
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
            printf("node:%d to node:%d bandwidth:%d > flow:%d\n", networkNodeIDStart, nextNetworkNodeID, previous->bandwidth, previous->flow);
            distToStart[nextNetworkNodeID] = previous->costPerGB;
        }
        preNetworkNodeID[nextNetworkNodeID] = networkNodeIDStart;
        // printf("---the previous node of node%d is %d\n", nextNetworkNodeID, networkNodeIDStart); 
        // printf("the costPerGB to node:%d is %d\n", nextNetworkNodeID, previous->costPerGB); 
    }
    // printf("start from node%d\n", networkNodeIDStart);
    isAccess[networkNodeIDStart] = true;
    distToStart[networkNodeIDStart] = 0;

    int tmpi = 0;
    // printf("the dist:\n");
    // printf("index:");
    // for(int i=0; i<networkNodeNum; i++){
    //     printf("%3d\t", i);
    // }
    // printf("\n");
    while(tmpi<=networkNodeNum){
        int minDisToStart = MAXINT;
        int minDisNetworkID = networkNodeIDStart;
        
        for(int i=0; i<networkNodeNum; i++){
            if(!isAccess[i] && minDisToStart > distToStart[i]){
                minDisNetworkID = i;
                minDisToStart = distToStart[i];
            }
        }
        // printf("%6d:",tmpi);
        // for(int i=0; i<networkNodeNum; i++){
        //     printf("%3d\t", distToStart[i]);
        // }
        // printf("\n");
        // printf("the min dist at node:%d is %d\n", minDisNetworkID, minDisToStart); 
        // printf("start from node:%d\n", minDisNetworkID);
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
            // printf("node:%d to node:%d bandwidth:%d, flow:%d\n", minDisNetworkID, nextNetworkNodeID, previous->bandwidth, previous->flow);
            if(previous->bandwidth <= previous->flow){
                // printf("WARNING: node:%d to node:%d bandwidth:%d, flow:%d\n", previous->networkNodeID1, previous->networkNodeID2, previous->bandwidth, previous->flow);
                continue;
            }
            // printf("judge node:%d\n", nextNetworkNodeID);
            // 如果当前节点距起始点距离+当前距下节点距离<下节点距起始点距离
            if(!isAccess[nextNetworkNodeID] && distToStart[minDisNetworkID] + previous->costPerGB < distToStart[nextNetworkNodeID]){
                // printf("W(%d)=%d > W(%d):%d+thisDist:%d=%d on %d\n", nextNetworkNodeID, distToStart[nextNetworkNodeID], minDisNetworkID, distToStart[minDisNetworkID], previous->costPerGB, distToStart[minDisNetworkID]+previous->costPerGB, minDisNetworkID);
                distToStart[nextNetworkNodeID] = distToStart[minDisNetworkID] + previous->costPerGB;
                preNetworkNodeID[nextNetworkNodeID] = minDisNetworkID;
                // printf("the previous node of node:%d is %d\n", nextNetworkNodeID, minDisNetworkID); 
            }
            // else if (!isAccess[nextNetworkNodeID]){
            //     printf("W(%d)=%d > W(%d):%d+thisDist:%d=%d on %d, no need update\n", nextNetworkNodeID, distToStart[nextNetworkNodeID], minDisNetworkID, distToStart[minDisNetworkID], previous->costPerGB, distToStart[minDisNetworkID]+previous->costPerGB, minDisNetworkID);
            // 
            // }
        }
    }
    return distToStart[networkNodeIDEnd];
}

void getNetworkIDSeqOnMinDist(int *preNetworkNodeID,int networkNodeIDStart,  int networkNodeIDEnd, vector<int> &networkNodeIDSeq)
{
    int nextNetworkNodeID = networkNodeIDEnd;
    printf("==getNetworkIDSeqOnMinDist\n");
    while(nextNetworkNodeID != networkNodeIDStart){
        networkNodeIDSeq.push_back(nextNetworkNodeID);
        if(nextNetworkNodeID == 128){
            printf("maybe have some error\n");
        }
        nextNetworkNodeID = preNetworkNodeID[nextNetworkNodeID];
    }
    networkNodeIDSeq.push_back(networkNodeIDStart);
}

int getMinFlowOnMinDist(vector<int> networkNodeIDSeq)
{
    int networkNodeIDInSeq = networkNodeIDSeq[networkNodeIDSeq.size()-1];
    int nextNetworkNodeID, tmpi = 1, minFlow = MAXINT;
    int curNetworkBandwidth = 0;
    printf("==getMinFlowOnMinDist\n");
    EdgePointer pointer = networkNode[networkNodeIDInSeq].nextEdge;
    while(pointer != NULL){
        curNetworkBandwidth = pointer->bandwidth - pointer->flow;
        // printf("curNetworkNodeID:%d\n", networkNodeIDInSeq);
        if (networkNodeIDInSeq  == pointer->networkNodeID1) {
            nextNetworkNodeID = pointer->networkNodeID2;
            pointer = pointer->edge1;
        } else {
            nextNetworkNodeID = pointer->networkNodeID1;
            pointer = pointer->edge2;
        }
        // printf("nextNetworkNodeID:%d\n", nextNetworkNodeID);
        if(nextNetworkNodeID == networkNodeIDSeq[networkNodeIDSeq.size()-tmpi-1]){
            // printf("the bandwidth node%d to node%d: %d, and minFlow = %d before updated\n", networkNodeIDInSeq, nextNetworkNodeID, curNetworkBandwidth, minFlow);
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
    printf("==updateFlow\n");
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
        // printf("nextNetworkNodeID:%d\n", nextNetworkNodeID);
        if(nextNetworkNodeID == networkNodeIDSeq[networkNodeIDSeq.size()-tmpi-1]){
            printf("the flow of node%d to node%d: %d change to %d, and bandwidth:%d\n", networkNodeIDInSeq, nextNetworkNodeID, previous->flow, minFlow+previous->flow, previous->bandwidth);
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

int getNetworkNodeIDLinkNum(NetworkNodePointer networkNode)
{
    EdgePointer pointer;
    pointer = networkNode->nextEdge;
    int networkNodeIDLinkNum = 0;
    printf("==getNetworkNodeIDLinkNum\n");
    while (pointer != NULL) {
        networkNodeIDLinkNum++;
        if ( networkNode->curNetworkNodeID == pointer->networkNodeID1) {
            pointer = pointer->edge1;
        } else {
            pointer = pointer->edge2;
        }
        // 如果在这里采用pointer访问数据，则可能此时pointer为空，所以会有segment core问题
    }
    printf("networkNode:%d, LinkNum:%d", networkNode->curNetworkNodeID, networkNodeIDLinkNum);
    return networkNodeIDLinkNum;
}

int getLinkCost(int networkNodeID1, int networkNodeID2)
{
    EdgePointer pointer, previous;
    pointer = networkNode[networkNodeID1].nextEdge;
    int nextNetworkNodeID = 0;
    printf("==getLinkCost:NodeID1:%d, NodeID2:%d\n", networkNodeID1, networkNodeID2);
    while(pointer != NULL){
        previous = pointer;
        if ( networkNodeID1 == pointer->networkNodeID1) {
            nextNetworkNodeID = pointer->networkNodeID2;
            pointer = pointer->edge1;
        } else {
            nextNetworkNodeID = pointer->networkNodeID1;
            pointer = pointer->edge2;
        }
        if(nextNetworkNodeID == networkNodeID2){
            printf("costPerGB between %d-%d:%d\n", networkNodeID1, networkNodeID2, previous->costPerGB);
            return previous->costPerGB;
        }
    }
    return -1;
}

// return cannot arrive networknode connected to usernode
int calcFlowPath(int *serverID, int serverNum)
{
    int *tmpForPreNetworkNodeID, *preNetworkNodeID;    
    printf("==calcFlowPath\n");
    tmpForPreNetworkNodeID = (int *)malloc(sizeof(int)*networkNodeNum);
    preNetworkNodeID = (int *)malloc(sizeof(int)*networkNodeNum);
    for(int i=0; i<userNodeNum; i++){
        int networkNodeIDStart, networkNodeIDEnd = userNode[i].conNetNodeID;
        printf("======================================================to userNode[%d]:%d\n", userNode[i].curUserNodeID, networkNodeIDEnd);
        // 服务器到一个用户节点的计算
        int iterationCount = 0;
        int networkNodeProvider = serverID[0];
        int minFlow, allFlow = 0; 
        bool directConnectFlag = false;
        vector<int> networkNodeIDSeq;
        memset(tmpForPreNetworkNodeID, -1, sizeof(int)*networkNodeNum);
        memset(preNetworkNodeID, -1, sizeof(int)*networkNodeNum);

        while(allFlow < userNode[i].bandwidth){
            int minCost = 32768;
            int searchServerCount = 0, noPathcount = 0;;
            int charNum;
            iterationCount++;
            memset(preNetworkNodeID, -1, sizeof(int)*networkNodeNum);
            for(int iForServerID=0; iForServerID<serverNum; iForServerID++){
                int tmpForMinCost = 32768;
                searchServerCount++;
                printf("=========serverID:%d to networkNodeID:%d\n", serverID[iForServerID], networkNodeIDEnd);
                networkNodeIDStart = serverID[iForServerID];
                memset(tmpForPreNetworkNodeID, -1, sizeof(int)*networkNodeNum);

                if(networkNodeIDStart == userNode[i].conNetNodeID){
                    directConnectFlag = true;
                    networkNodeProvider = networkNodeIDStart;
                    goto DirectConnect;
                    break;
                } else {
                    directConnectFlag = false;
                }

                // 计算哪台服务器到用户点1的路径最优
                // 先获得最优，再获取次优，根据ncost+link小者为合适(n为经过的节点数，所经过节点周围链路总和)
                tmpForMinCost = dijkstra(networkNodeIDStart, networkNodeIDEnd, tmpForPreNetworkNodeID);
                printf("serverID:%d, minCost:%d\n", serverID[iForServerID], tmpForMinCost);

                if(tmpForMinCost == 128){
                    noPathcount++;
                    printf("WARNING: %d-%d have no path\n", networkNodeIDStart, networkNodeIDEnd);
                    continue;
                }
                // =====最终的结果
                if(minCost > tmpForMinCost){
                    networkNodeProvider = networkNodeIDStart;
                    minCost = tmpForMinCost;
                    memcpy(preNetworkNodeID, tmpForPreNetworkNodeID, sizeof(int)*networkNodeNum);
                    // printf("select serverID:%d to userNode:%d, minCost:%d\n", serverID[iForServerID], networkNodeIDEnd, minCost);
                }
                networkNodeIDSeq.clear();
                getNetworkIDSeqOnMinDist(tmpForPreNetworkNodeID, networkNodeIDStart, networkNodeIDEnd, networkNodeIDSeq);
                // printf("previous node:\n");
                // for(int i=0; i<networkNodeNum; i++){
                //     printf("%3d\t", i);
                // }
                // printf("\n");
                // for(int i=0; i<networkNodeNum; i++){
                //     printf("%3d\t", preNetworkNodeID[i]);
                // }
                // printf("\n");
                printf("server %d result: cur minCost:%d, %d-%d seq's minCost: %d and the seq:",searchServerCount, minCost, networkNodeIDStart, networkNodeIDEnd, tmpForMinCost);
                for(int i=0; i<(int)networkNodeIDSeq.size(); i++){
                    printf("%d\t", networkNodeIDSeq[networkNodeIDSeq.size()-1-i]);
                }
                printf("\n");
            }
            if(noPathcount >= serverNum){
                printf("ERROR: no server to the userNode[%d]:%d flow:%d  > allFlow:%d\n", i, userNode[i].conNetNodeID, userNode[i].bandwidth, allFlow);
                return userNode[i].conNetNodeID;
            }
            // printf("previous Node ID:\n");
            // for(int i=0; i<networkNodeNum; i++){
            //    printf("%3d\t", i);
            // }
            // printf("\n");
            // for(int i=0; i<networkNodeNum; i++){
            //    printf("%3d\t", preNetworkNodeID[i]);
            // }
            // printf("\n");
            // 获得start to end的最短路径上的ID序列
            networkNodeIDSeq.clear();
            getNetworkIDSeqOnMinDist(preNetworkNodeID, networkNodeProvider, networkNodeIDEnd, networkNodeIDSeq);
            // 获得该序列上点间最小的容量值
            minFlow = getMinFlowOnMinDist(networkNodeIDSeq);
            // 根据用户需求调整最小流
            if(allFlow + minFlow > userNode[i].bandwidth){
                printf("minFlow:%d + allFlow:%d > userNode[%d].bandwidth:%d, so update minFlow:%d and minCost:%d\n", minFlow, allFlow, i, userNode[i].bandwidth, userNode[i].bandwidth-allFlow, minCost);
                minFlow = userNode[i].bandwidth - allFlow;
                allFlow = userNode[i].bandwidth;
            } else {
                allFlow += minFlow;
            }
            printf("allCost:%d + addCost:%d = %d", allCost, minCost*minFlow, allCost+minCost*minFlow);
            // 根据最小流更新路径上的当前流量值
            updateFlow(networkNodeIDSeq, minFlow);

            printf("the %d-%d min flow: %d\n", networkNodeProvider, networkNodeIDEnd, minFlow);
            printf("select %d-%d seq:", networkNodeProvider, networkNodeIDEnd);
            // 将关注点写入内存缓冲区
            for(int i=0; i<(int)networkNodeIDSeq.size(); i++){
                printf("%d\t",networkNodeIDSeq[networkNodeIDSeq.size()-1-i]);
                // itoa(networkNodeIDSeq[networkNodeIDSeq.size()-1-i], topoFileCurPointer, 10);
                charNum = sprintf(topoFileCurPointer, "%d ", (char)networkNodeIDSeq[networkNodeIDSeq.size()-1-i]);
                topoFileCurPointer += charNum;
                // for(int i=0; i<10; i++){
                //     printf("%c", topo_file[i]);
                // }
                // printf("\n");
            }
            allCost += minCost*minFlow;
            printf("\n");
            // 计算租用费用
DirectConnect:
            if(directConnectFlag){
                printf("serverID:%d connect to userNode[%d]:%d directly, and the flow:%d\n", networkNodeProvider, userNode[i].curUserNodeID, userNode[i].conNetNodeID, userNode[i].bandwidth);
                minFlow = userNode[i].bandwidth;
                charNum = sprintf(topoFileCurPointer, "%d ", (char)networkNodeProvider);
                topoFileCurPointer += charNum;

                minCost = 0;
            }
            // 加入用戶Node
            charNum = sprintf(topoFileCurPointer, "%d ", (char)userNode[i].curUserNodeID);
            topoFileCurPointer += charNum;
            // 加入流量值
            charNum = sprintf(topoFileCurPointer, "%d", (char)minFlow);
            topoFileCurPointer += charNum;
            *(topoFileCurPointer++) = '\n';
            printf("\n");
            
            // printf("start, end, band, cost, flow\n");
            // for (int i=0; i<NETWORK_NODE_MAX_NUM; i++) {
            //     if (networkNode[i].nextEdge == NULL) {
            //         break;
            //     }
            //     printf("NetworkNode %d: ", i);
            //     printNetworkNodeInfo(&networkNode[i]);
            // }
            printf("=======%d result:allFlow:%d, userNode[%d].bandwidth:%d\n", iterationCount, allFlow, userNode[i].curUserNodeID,  userNode[i].bandwidth);
            networkPathNum++;
            printf("--------------------networkPathNum:%d\n", networkPathNum);
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
}

int addNoPathMethod()
{
    // 找到视频服务器位置
    // 1. 用户节点相连周围线路流量总和<用户需求为服务器
    // 2. 最大流量前serverNum名为服务器，服务器直接相连的*0.8->*0.6...
    int serverNum = 3;
    int *serverID = (int *)malloc(sizeof(int)*serverNum);
    // 获取初始的serverNum服务器位置
    getServerID(serverID, serverNum);
    printf("max flow networkNodeID:");
    for (int i=0; i<serverNum; i++) {
        printf("%d, ", serverID[i]);
    }
    printf("\n");
    int addServerID = -1;
    int *tmpForServerID = (int *)malloc(sizeof(int)*serverNum);
    do{
        if(addServerID != NOT_NETWORK_NODE_ID){
            printf("have no way to %d, add it to serverID\n", addServerID);
            memcpy(tmpForServerID, serverID, sizeof(int)*serverNum);
            initForRestart();
            free(serverID);
            serverNum++;
            serverID = (int *)malloc(sizeof(int)*serverNum);
            memcpy(serverID, tmpForServerID, sizeof(int)*serverNum);
            free(tmpForServerID);
            tmpForServerID = (int *)malloc(sizeof(int)*serverNum);  // 让备份部分空间也+1

            serverID[serverNum-1] = addServerID;

            memset(topo_file_master, 0, sizeof(char)*NETWORK_NODE_MAX_NUM);
            strcpy(topo_file_master, "      \n\n");
            topoFileCurPointer = topo_file_master + 8;
        }

        printf("server location ID:");
        for (int i=0; i<serverNum; i++) {
            printf("%d, ", serverID[i]);
        }
        printf("\n");
        // max flow min cost
        // addServerID = calcFlowPath(serverID, serverNum);
    }while(NOT_NETWORK_NODE_ID != (addServerID = calcFlowPath(serverID, serverNum)));

    printf("congratulation, have an answer!~_~\n");
    printf("PathNum:%d, RentCost:%d, ServerNum:%d, AllCost:%d\n", networkPathNum, allCost, serverNum, allCost+300*serverNum);
    printf("and serverID:");
    for(int i=0; i<serverNum; i++){
        printf("%d\t", serverID[i]);
    }
    printf("\n");
    char tmp[6];
    int charNum = sprintf(tmp, "%d", networkPathNum);
    for(int i=0; i<charNum; i++){
        *(topo_file_master+i) = tmp[i];
    }
    free(tmpForServerID);
    return (allCost+300*serverNum);
}

int maxFlowServerMethod()
{  
    int noPathServerID = 0;
    // 找到视频服务器位置
    // 1. 用户节点相连周围线路流量总和<用户需求为服务器
    // 2. 最大流量前serverNum名为服务器，服务器直接相连的*0.8->*0.6...
    int serverNum = 3;
    int *serverID = (int *)malloc(sizeof(int)*serverNum);
    topoFileCurPointer = topo_file_maxFlow;
    do{
        if(noPathServerID != 0){
            printf("have no way to %d, research serverID\n", noPathServerID);
            initForRestart();
            free(serverID);
            serverNum++;
            serverID = (int *)malloc(sizeof(int)*serverNum);

            memset(topo_file_maxFlow, 0, sizeof(char)*NETWORK_NODE_MAX_NUM);
            strcpy(topo_file_maxFlow, "      \n\n");
            topoFileCurPointer = topo_file_maxFlow + 8;
        }
        getServerID(serverID, serverNum);
        printf("max flow networkNodeID:");
        for (int i=0; i<serverNum; i++) {
            printf("%d, ", serverID[i]);
        }
        printf("\n");
    }while(NOT_NETWORK_NODE_ID != (noPathServerID = calcFlowPath(serverID, serverNum)));

    printf("congratulation, have an answer!~_~\n");
    printf("PathNum:%d, RentCost:%d, ServerNum:%d, AllCost:%d\n", networkPathNum, allCost, serverNum, allCost+300*serverNum);
    printf("and serverID:");
    for(int i=0; i<serverNum; i++){
        printf("%d\t", serverID[i]);
    }
    printf("\n");
    char tmp[6];
    int charNum = sprintf(tmp, "%d", networkPathNum);
    for(int i=0; i<charNum; i++){
        *(topo_file_maxFlow+i) = tmp[i];
    }
    return (allCost+300*serverNum);
}
//你要完成的功能总入口
void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename)
{
    // print_time("deploy_server");
    init();

    // 读取网络节点参数
    readNetworkNodeInfo(topo, line_num);
    // 读取用户节点信息
    readUserNodeInfo(topo, line_num);
    printf("=======================\n");
    printf("NetworkNodeNum:%d, NetworkLinkNum:%d, UserNodeNum:%d, AllCost:%d\n", networkNodeNum, networkLinkNum, userNodeNum, allCost);
    int i = 0;
    printf("start, end, band, cost, flow\n");
    for (i=0; i<NETWORK_NODE_MAX_NUM; i++) {
        if (networkNode[i].nextEdge == NULL) {
            break;
        }
        printf("NetworkNode %d: ", i);
        printNetworkNodeInfo(&networkNode[i]);
    }
    printf("LinkItemNum: %d, NetworkNodeNum = %d\n", tmp, i);
    //===========节点信息读入结构体完毕
    // 遗传算法进化
    //
    // 第一种方法：采用过不去就是服务器master
    int firstMethodCost = 100000;
    initForRestart();
    firstMethodCost = addNoPathMethod();
    // 第二种方法：采用支持流量最大为服务器maxFlow
    int secondMethodCost = 10000;
    initForRestart();
    secondMethodCost = maxFlowServerMethod();
    // printf("max flow networkNodeID:");
    // for (int i=0; i<serverNum; i++) {
    //     printf("%d, ", serverID[i]);
    // }
    // printf("\n");

    printf("firstMethodCost:%d, secondMethodCost:%d\n", firstMethodCost, secondMethodCost);
    if( firstMethodCost > secondMethodCost){
        printf("select maxFlow\n");
        topoFileCurPointer = topo_file_maxFlow;
    } else {
        printf("select master\n");
        topoFileCurPointer = topo_file_master;
    }
    
	// 直接调用输出文件的方法输出到指定文件中(ps请注意格式的正确性，如果有解，第一行只有一个数据；第二行为空；第三行开始才是具体的数据，数据之间用一个空格分隔开)
	write_result((const char *)topoFileCurPointer, filename);
}
