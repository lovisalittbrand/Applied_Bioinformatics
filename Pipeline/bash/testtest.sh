#!/bin/bash

file="data/cog_list.txt"
vec=()
while read LINE  
do  
  echo $LINE
  vec+=($LINE)
    
done <$file

echo "heeeej"
echo ${vec[@]} 
