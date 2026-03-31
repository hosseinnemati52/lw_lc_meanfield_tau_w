ROOT_DIR="${1:-.}"

find "$ROOT_DIR" -type d -name "GD_temp" | while read -r dir; do
    echo "Processing: $dir"
    
    parent_dir=$(dirname "$dir")
    zip_name="${parent_dir}/GD_temp.zip"

    # Zip from inside the parent directory
    ( cd "$parent_dir" && zip -r "GD_temp.zip" "GD_temp" )
    
    # Check the integrity of the zip file
    if unzip -t "$zip_name" > /dev/null; then
        echo "Zip file is valid. Deleting original folder: $dir"
        rm -rf "$dir"
    else
        echo "Zip integrity check failed for: $zip_name"
        echo "Skipping deletion of original folder: $dir"
    fi
done
