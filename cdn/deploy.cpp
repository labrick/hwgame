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

#define SERVER_NUM  3

#define NETWORK_PATH_MAX_NUM 50000
#define NETWORK_MAX_NUM_PER_PATH 1000

// 网络路径数量
int networkPathNum = 0;
char topo_file[NETWORK_PATH_MAX_NUM*NETWORK_MAX_NUM_PER_PATH]; 
char *topoFileCurPointer = topo_file;

using namespace std;

typedef struct UserNode{
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
    strcpy(topo_file, "      \n\n");
    topoFileCurPointer = topo_file + 8;
}

// 用户需求量的大的排在前面以便优先处理
void bubbleSortUserNode()
{
    for(int i=0; i<userNodeNum; i++){
        for(int j=i; j<(userNodeNum-1); j++){
            if(userNode[j].bandwidth < userNode[j+1].bandwidth){
                UserNode tmp = userNode[j+1];
                userNode[j+1] = userNode[j];
                userNode[j] = tmp;
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

void bubbleSort(int *arrayHead, int *relationArray, int memNum)
{
    for(int i=0; i<memNum; i++){
        for(int j=i; j<(memNum-1); j++){
            if(arrayHead[j] > arrayHead[j+1]){
                int tmp = arrayHead[j+1];
                arrayHead[j+1] = arrayHead[j];
                arrayHead[j] = tmp;
                tmp = relationArray[j+1];
                relationArray[j+1] = relationArray[j];
                relationArray[j] = tmp;
            }
        }
    }
}

void getServerID(int *serverID, int serverNum)
{
    int *allFlow = (int *)malloc(sizeof(int)*serverNum);
    vector<int> mustServerID, mustServerAllFlow;
    memset(allFlow, 0, sizeof(int)*serverNum);
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
            printf("networkNode[%d] flow is: %d\n", i, tmpForAllFlow);
            for(int tmpi=0; tmpi<userNodeNum; tmpi++){
                if((networkNode[i].curNetworkNodeID == userNode[i].conNetNodeID) && \
                        (tmpForAllFlow <= userNode[i].bandwidth)){
                    mustServerID.push_back(userNode[i].conNetNodeID);
                    mustServerAllFlow.push_back(tmpForAllFlow);
                }
            }
            if (tmpForAllFlow > *allFlow) {
                // printf("networkNode[%d] is the max flow: %d\n", i, tmpForAllFlow);
                *allFlow = tmpForAllFlow;
                *serverID = i;
                // printf("before bubbleSort:\n");
                // for(int i=0; i<serverNum; i++){
                //     printf("%d, ", allFlow[i]);
                // }
                bubbleSort(allFlow, serverID, serverNum);
                // printf("\nafter bubbleSort:\n");
            }
        }   
    }
    for(int i=0; i<(int)mustServerID.size(); i++){
        printf("serverID:%d is mustServerID\n", mustServerID[mustServerID.size()-1]);
        serverID[serverNum-i] = mustServerID[mustServerID.size()-1];
    }
    for(int i=0; i<serverNum; i++){
        printf("serverID:%d\tmaxFlow:%d\n", serverID[i], allFlow[i]);
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
                printf("WARNING: node:%d to node:%d bandwidth:%d, flow:%d\n", previous->networkNodeID1, previous->networkNodeID2, previous->bandwidth, previous->flow);
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

void calcFlowPath(int *serverID, int serverNum)
{
    int *tmpForPreNetworkNodeID, *preNetworkNodeID;    
    printf("==calcFlowPath\n");
    tmpForPreNetworkNodeID = (int *)malloc(sizeof(int)*networkNodeNum);
    preNetworkNodeID = (int *)malloc(sizeof(int)*networkNodeNum);
    for(int i=0; i<userNodeNum; i++){
        int networkNodeIDStart, networkNodeIDEnd = userNode[i].conNetNodeID;
        printf("======================================================to userNode[%d]:%d\n", i, networkNodeIDEnd);
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
                printf("no server to the userNode[%d]:%d flow:%d  > allFlow:%d\n", i, userNode[i].conNetNodeID, userNode[i].bandwidth, allFlow);
                break;
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
            printf("\n");
DirectConnect:
            if(directConnectFlag){
                printf("serverID:%d connect to userNode[%d]:%d directly, and the flow:%d\n", networkNodeProvider, i, userNode[i].conNetNodeID, userNode[i].bandwidth);
                minFlow = userNode[i].bandwidth;
                charNum = sprintf(topoFileCurPointer, "%d ", (char)networkNodeProvider);
                topoFileCurPointer += charNum;
            }
            // 加入用戶Node
            charNum = sprintf(topoFileCurPointer, "%d ", (char)i);
            topoFileCurPointer += charNum;
            // 加入流量值
            charNum = sprintf(topoFileCurPointer, "%d ", (char)minFlow);
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
            printf("=======%d result:allFlow:%d, userNode[%d].bandwidth:%d\n", iterationCount, allFlow, i,  userNode[i].bandwidth);
            networkPathNum++;
            printf("--------------------networkPathNum:%d\n", networkPathNum);
            if(directConnectFlag)
                break;
        }
    }
}

//你要完成的功能总入口
void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename)
{
    // print_time("deploy_server");
    init();

	// 需要输出的内容
    // char * topo_file = (char *)"17\n\n0 8 0 20\n21 8 0 20\n9 11 1 13\n21 22 2 20\n23 22 2 8\n1 3 3 11\n24 3 3 17\n27 3 3 26\n24 3 3 10\n18 17 4 11\n1 19 5 26\n1 16 6 15\n15 13 7 13\n4 5 8 18\n2 25 9 15\n0 7 10 10\n23 24 11 23";
    
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
    // 找到视频服务器位置：根据最大流量输出来计算
    int serverID[SERVER_NUM]={-1, };
    getServerID(serverID, SERVER_NUM);
//     serverID[0] = 0;
//     serverID[1] = 1;
//     serverID[2] = 21;
    printf("max flow networkNodeID:");
    for (int i=0; i<SERVER_NUM; i++) {
        printf("%d, ", serverID[i]);
    }
    printf("\n");
    // max flow min cost
    calcFlowPath(serverID, SERVER_NUM);
    
    
	// 直接调用输出文件的方法输出到指定文件中(ps请注意格式的正确性，如果有解，第一行只有一个数据；第二行为空；第三行开始才是具体的数据，数据之间用一个空格分隔开)
    char tmp[6];
    int charNum = sprintf(tmp, "%d", networkPathNum);
    for(int i=0; i<charNum; i++){
        *(topo_file+i) = tmp[i];
    }
	write_result((const char *)topo_file, filename);
}
