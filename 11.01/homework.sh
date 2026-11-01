#!/bin/bash

# 进入目标文件夹（注意：判断 -f / -d 是基于当前目录的）
cd bash_homework

# 每次运行先清空结果文件，避免重复追加
> filenames.txt
> dirname.txt

# 用命令替换把 ls 的结果存到变量
CUR_DIR=`ls`

# 逐个遍历
for val in $CUR_DIR
do
    if [ -f $val ]; then
        # >> 表示追加写入，> 会覆盖
        echo $val >> filenames.txt
    elif [ -d $val ]; then
        echo $val >> dirname.txt
    fi
done

exit 0
