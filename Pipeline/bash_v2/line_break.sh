#!/bin/bash
#Script that removes line breaks in FASTA-files

input=$1
output=$2

awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}' $input > $output
