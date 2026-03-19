#!/bin/bash

# Only add to safe.directory if it's not already there
if ! git config --global --get-all safe.directory | grep -q "$(pwd)"; then
    git config --global --add safe.directory "$(pwd)"
fi

# This command tells Git: "When I commit, convert to LF. When I check out, leave it alone."
git config --global core.autocrlf input

# Get the current username to identify the computer
HOST_NAME=$(hostname)

echo "Detected machine: $HOST_NAME"

if [ "$HOST_NAME" == "DESKTOP-4BM2SLT" ]; then
    echo "--- Desktop detected ---"
    # Create the desktop venv if it doesn't exist
    # (Update the path below to match your desktop's Python path!)
    if [ ! -d ".venv_desktop" ]; then
        echo "Creating desktop virtual environment..."
        python -m venv .venv_desktop
    fi
    source .venv_desktop/Scripts/activate
else
    echo "--- Laptop detected ---"
    # Create the laptop venv if it doesn't exist
    if [ ! -d ".venv_laptop" ]; then
        echo "Creating laptop virtual environment..."
        "/c/Users/jaume/AppData/Local/Programs/Python/Python314/python.exe" -m venv .venv_laptop
    fi
    source .venv_laptop/Scripts/activate
fi

# Ensure pip is ready and show current version
python --version