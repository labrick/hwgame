#!/bin/bash

##################参数测试脚本使用说明####################################
#1.首先需要配置好各个文件夹的地址
#2.设置好温度
#3.设置好温度下降率
#4.设置case的数值


#######################（1）文件夹地址####################################
casePlace=/mnt/hdd2/yanan/huawei/case_example/second/advance
workPlace=/mnt/hdd2/yanan/huawei/zheng_test
sourceFile=$workPlace/gcc
outDir=$workPlace/result
mkdir -p $outDir

#######################（2）温度参数设置##################################
sa_T_original=20
sa_T_end=200
sa_T_unit=20
((T_count=($sa_T_end-$sa_T_original)/$sa_T_unit - 1))

#######################（3）温度变化率设置################################
sa_T_list=(0.9 0.99 0.999 0.9999 0.99999)
sa_Tdelta_original=0.9
# sa_Tdelta_end=0.99999
# sa_Tdelta_unit=0.05
((Tdelta_count=5))
#######################（4）随机加点概率##################################
addProbabilityUp=0.1

#######################（5）随机删点概率##################################
delProbabilityUp=0.9

#######################（6）case的数值####################################
case_count=0 #比实际case数小1


#ga_Tdelta=0.995
#ga_Tdelta_unit=0.05
#ga_Tdelta_line=936

############################ 定位数据在源程序中的行号 ####################################
sa_T_line=`grep "double T = " -n $sourceFile/cdn/deploy.cpp | grep -o '^[0-9]\{1,9\}' `
sa_Tdelta_line=`grep "double Tdelta = " -n $sourceFile/cdn/deploy.cpp | grep -o '^[0-9]\{1,9\}' `
sa_addP_line=`grep "double addProbabilityUp = " -n $sourceFile/cdn/deploy.cpp | grep -o '^[0-9]\{1,9\}' `
sa_delP_line=`grep "double delProbabilityUp = " -n $sourceFile/cdn/deploy.cpp | grep -o '^[0-9]\{1,9\}' `

#echo "T: $sa_T_line"
#echo "Tdelta: $sa_Tdelta_line"
#echo "addProbabilityUp: $sa_addP_line"
#echo "delProbabilityUp: $sa_delP_line"


############################ 正式程序开始启动 ####################################


#i的索引号就是case的索引号，每个case都是9个
for i in $(seq 0 $case_count)
do
    caseSum=0
    finalFile=$outDir/final_case$i.txt
    echo "case$i" > $finalFile
    echo "addProbabilityUp: $addProbabilityUp" > $finalFile
    echo "delProbabilityUp: $delProbabilityUp" > $finalFile

    echo "***********************************************************************" >> $finalFile
    echo -e "\t\t\t\t\t\t\t1\t\t2\t\t3\t\t4\t\t5\t\tav" >> $finalFile
    echo "-----------------------------------------------------------------------" >> $finalFile
    #k用来更新T-delta
    for sa_Tdelta in 0.9 0.99 0.999 0.9999 0.99999       #$(seq 0 $Tdelta_count)
    do
        #k用来更新T
        sa_T=$sa_T_original
        for T in $(seq 0 $T_count) 
        do
            sed -i "${sa_T_line}s/[0-9]\{1,9\}/$sa_T/"  $sourceFile/cdn/deploy.cpp
            sed -i "${sa_Tdelta_line}s/[0-9].[0-9]\{1,9\}/$sa_Tdelta/"  $sourceFile/cdn/deploy.cpp
            sed -i "${sa_addP_line}s/[0-9].[0-9]\{1,9\}/$addProbabilityUp/"  $sourceFile/cdn/deploy.cpp
            sed -i "${sa_delP_line}s/[0-9].[0-9]\{1,9\}/$delProbabilityUp/"  $sourceFile/cdn/deploy.cpp
            #sed -i "${ga_Tdelta_line}s/[0-9].[0-9]\{1,5\}/$ga_Tdelta/"  $sourceFile/cdn/deploy.cpp
            
            #j的索引号是每个case循环的次数，这里设置为5次
            #修改变量完成之后，开始编译程序
            $sourceFile/build.sh
            for j in $(seq 0 1)
            do
                echo "case$i 第 $j 次计算" 
                #判断上一次执行成功以后，进行下次参数替换，然后继续执行程序
                #对指定的参数进行修改
                
                echo "T: $sa_T"
                echo "Tdelta: $sa_Tdelta"
                echo "addProbabilityUp: $addProbabilityUp"
                echo "delProbabilityUp: $delProbabilityUp"

                $sourceFile/bin/cdn $casePlace/case$j.txt  $outDir/result$i-D$sa_Tdelta-T$sa_T-$j.txt > $outDir/log_case$i-D$sa_Tdelta-T$sa_T-$j.txt 
                
                #从打印信息中提取出最小cost
                awk '/congratulation, have an answer!~_~/ { print x }; { x = $0 }'  $workPlace/result/log_case$i-D$sa_Tdelta-T$sa_T-$j.txt | grep -o '[0-9]\{1,9\}' > $outDir/temp.txt
                tmp[$j]=`grep -o '[0-9]\{1,9\}' $outDir/temp.txt`
                echo ${tmp[$j]}
            done
            #计算5次开销的平均值
            ((caseSum=${tmp[0]}+${tmp[1]}))
            tmp[2]=0
            tmp[3]=0
            tmp[4]=0
            #((caseSum=${tmp[0]}+${tmp[1]}+${tmp[2]}+${tmp[3]}+${tmp[4]})) 
            caseNum=2
            ((caseAverage=$caseSum/$caseNum))
            echo $caseAverage
           

            echo -e "T: $sa_T | delta: $sa_Tdelta | \t${tmp[0]}\t${tmp[1]}\t\t\t\t\t\t\t\t$caseAverage" >> $finalFile
            #echo -e "T: $sa_T | delta: $sa_Tdelta | \t${tmp[0]}\t${tmp[1]}\t${tmp[2]}\t${tmp[3]}\t${tmp[4]}\t$caseAverage" >> $finalFile
            #echo -e "ga_Tdelta:       $ga_Tdelta     | \t\t${tmp[0]}\t${tmp[1]}\t${tmp[2]}\t${tmp[3]}\t${tmp[4]}\t$caseAverage" >> $finalFile
            
            #计算完了每次的开销之后更新T或者T-delta进入下一次计算     
            sa_T=`expr $sa_T + $sa_T_unit`
            #sa_Tdelta=`expr $sa_Tdelta - $sa_Tdelta_unit`
        done 
        echo "-----------------------------------------------------------------------" >> $finalFile
    done
done

