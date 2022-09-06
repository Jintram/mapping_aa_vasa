#!/bin/bash

step=$1

if [[ "${step}" == *"(1)"* || "${step}" == "all" ]]; then
	echo "yes"
else
	echo "no" 
fi





