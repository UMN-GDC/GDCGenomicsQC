#!/bin/bash

# Check if the number of arguments is correct
if [ $# -ne 3 ]; then
    echo "Usage: $0 <file> <line_number> <new_content>"
    exit 1
fi

file=$1
line_number=$2
new_content=$3

# Check if the file exists
if [ ! -f "$file" ]; then
    echo "Error: File $file does not exist!"
    exit 1
fi

# Check if the line number is valid
if ! [[ "$line_number" =~ ^[0-9]+$ ]]; then
    echo "Error: Invalid line number!"
    exit 1
fi

# Check if the line number is within the file range
file_lines=$(wc -l < "$file")
if [ "$line_number" -gt "$file_lines" ]; then
    echo "Error: Line number is greater than the number of lines in the file!"
    exit 1
fi

# Temporarily store the file content in an array
mapfile -t lines < "$file"

# Replace the content of the specified line with the new content
lines["$((line_number-1))"]=$new_content

# Write the modified content back to the file
printf "%s\n" "${lines[@]}" > "$file"

echo "Content of line $line_number in $file replaced with: $new_content"


# Replace the content of the specified line with the new content
# sed -i "${line_number}s/.*/$new_content/" "$file"

# echo "Content of line $line_number in $file replaced with: $new_content"
