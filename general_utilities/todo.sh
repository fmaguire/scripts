#!/bin/zsh

usage="todo {add,remove} 'task'"
todo_file=/home/fin/.zsh/todo.txt

if [ $1 = "add" ]; then
    echo $2 >> $todo_file
elif [ $1 = "remove" ]; then
    sed -i "/${2}/d" $todo_file
fi

cat $todo_file
