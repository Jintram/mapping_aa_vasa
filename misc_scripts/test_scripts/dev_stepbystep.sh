#!/bin/bash


step=$1

if [ "$step" = "" ] || [ "$step" = "all" ]; then
  echo "Step was not set, setting to all"
  step="all"
fi

echo "Executing step $step"

if [ "$step" = "1" ] || [ "$step" = "all" ]; then
  echo "Performing step 1"
fi

if [ "$step" = "2" ] || [ "$step" = "all" ]; then
  echo "Performing step 2"
fi

if [ "$step" = "3" ] || [ "$step" = "all" ]; then
  echo "Performing step 3"
fi

