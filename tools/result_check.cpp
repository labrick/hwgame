#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <errno.h>
#include <unistd.h>
#include <signal.h>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

int read_file(char buff[1024][1024], char * filename){
    ifstream fp(filename);
    if(!fp.is_open()){
        printf("Fail to open file %s.\n", filename);
        return 0;
    }
    printf("Open file %s is ok.\n", filename);
    unsigned int cnt =0;
    while(!fp.eof()){
        fp.getline(buff[cnt],100);
        //cout << buff[cnt] << endl;
        cnt++;
    }
    return cnt;
}

struct edge{
    int S;
    int D;
    int width;
};

int main(int argc, char *argv[])
{
    int line_num;
    char *topo_file = argv[1];
    char *result_file = argv[2];
    char buff[1024][1024];
    char resultBuff[1024][1024];
    line_num = read_file(buff, topo_file);
    printf("There are %d lines in %s.\n", line_num, topo_file);
    //for(int i = 0; i<=line_num; i++)
    //    cout<<buff[i]<<endl;
    //string a=buff[0];
    //cout << a << endl;
    //
    //
    vector <edge> Saler;
    vector <edge> Eater;
    int nodeNum = 0;
    int edgeNum = 0;
    int eaterNum = 0;
    for(int i=0; i<line_num; i++){
        const char *d = " \t";
        char *p;
        int num=0;
        int a[100];

        if(i==1 || i==3)
            continue;
        if(i == 0){
            p = strtok(buff[0],d);
            while(p){
                a[num] = atoi(p);
            //    printf("topo_file:i=%d,p=%s,num= %d\n",i, p, num);
                p=strtok(NULL,d);
                ++num;
            }
            nodeNum = a[0]; 
            edgeNum = a[1];
            eaterNum = a[2];
            continue;
        }

        p = strtok(buff[i],d);
        while(p){
           // printf("topo_file:i=%d,p=%s,num= %d\n",i, p, num);
            a[num] = atoi(p);
            p=strtok(NULL,d);
            ++num;
        }
        if(i == 4+edgeNum || i > 4+edgeNum+eaterNum)
            continue;
        printf("topo_file:i=%d,num=%d\n",i,num);
        if(num == 4){
            edge thisEdge = {a[0], a[1], a[2]};
            edge thatEdge = {a[1], a[0], a[2]};
            Saler.push_back(thisEdge);
            Saler.push_back(thatEdge);
            continue;
        }
        if(num == 3){
            edge thisEater = {a[0], a[1], a[2]};
            Eater.push_back(thisEater);
        }        
        //for(int i=0; i<num; i++)
        //    cout << a[i] << endl;
    }

    int salerSize = Saler.size();
    int eaterSize = Eater.size();
    //cout << Saler.size() << endl;
    //cout << Eater.size() << endl;
    
    line_num = read_file(resultBuff,result_file);
    printf("There are %d lines in %s.\n", line_num, result_file);
    printf("Result use edges:\n");
    int resultEdge = 0;
    //for(int k=0; k<=line_num; k++){
    for(int k=0; k<line_num; k++){
        const char *d = " \t";
        char *p;
        int num=0;
        int b[100];
        if (k==0){
            p = strtok(resultBuff[k],d);
            //printf("p=%s\n",p);
            while(p){
                printf("p=%s\n",p);
                b[num] = atoi(p);
                //printf("%s\n",p)
                p = strtok(NULL,d);
                num++;
                printf("n=%d\n",num);
            }
            //for(int i=0; i<num; i++)
            //    cout << b[i] << endl;
            resultEdge = b[0];
            continue;
        }
        if (k==1)
            continue;
        if (k > 2+resultEdge)
            continue;
        p = strtok(resultBuff[k],d);
        while(p){
            b[num] = atoi(p);
            p = strtok(NULL,d);
            num++;
        }
        printf("result_file:k=%d,num=%d\n",k,num);
        if(num != 3)
            for(int i = 0; i < num-3; i++){
                for(int j = 0; j < salerSize; j++){
                    if(Saler[j].S == b[i] && Saler[j].D == b[i+1]){
                        int temp = Saler[j].width;
                        Saler[j].width -= b[num -1];
                        if(Saler[j].width < 0){
                            printf("failed line %d, location %d,%d to %d, width %d, need %d\n",k+1,i+1,b[i], b[i+1],Saler[j].width,b[num-1]);
                            //return 0;
                        }
                       //  printf("%d to %d, with_change:%d to %d\n", Saler[j].S, Saler[j].D, temp, Saler[j].width);
                    }        
                }
            }
        for(int j = 0; j < eaterSize; j++){
            if(Eater[j].S == b[num-2])
                Eater[j].width -= b[num-1];
        }
    }
    for(int j = 0; j < eaterSize; j++){
        if(Eater[j].width != 0){
            printf("eater%d is false, width%d\n",j,Eater[j].width);
            //return 0;
        }
    }
    // printf("All edges:\n");
    // for(int i=0; i < salerSize; i++){
    //     printf("%d to %d, width:%d\n",Saler[i].S, Saler[i].D,Saler[i].width);
    // }
    printf("True\n");
    return 0;

}


