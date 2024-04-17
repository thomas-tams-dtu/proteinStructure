#!/bin/bash

# Specify the directory name
alphamissense_dir="data/alphamissense/P28698"

# Check if the directory exists
if [ ! -d "$alphamissense_dir" ]; then
    # If the directory does not exist, create it
    mkdir -p "$alphamissense_dir"
    echo "Directory created: $alphamissense_dir"
else
    # If the directory already exists, display a message
    echo "Directory '$alphamissense_dir' already exists."
fi

echo "Moving to $alphamissense_dir"
cd $alphamissense_dir

alphamissense_cmd="nohup /home/ctools/alphamissense/alphamissense P28698 &"
echo "Running $alphamissense_cmd"
$alphamissense_cmd
