#!/bin/bash
# Example helper script for the pipeline

echo "Running custom analysis script..."

# Example: simple statistics
if [ -f "$1" ]; then
    echo "File: $1"
    echo "Lines: $(wc -l < "$1")"
fi

echo "Done!"