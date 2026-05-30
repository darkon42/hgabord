#!/bin/bash

# --- Configuration ---
# 1. Define the predefined source directory you want to zip
SOURCE_DIR="/home/chris/Dropbox/workspace/hgabord" 

# 2. Define the predefined destination directory
DEST_DIR="/home/chris/archives/hgabord"       

# --- Execution ---
# 3. Generate the timestamp (Format: YYYY-MM-DD_HH-MM-SS)
TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")

# 4. Set the name of the zip file
ZIP_FILENAME="hgabord_${TIMESTAMP}.zip"

# 5. Ensure the destination directory exists (creates it if it doesn't)
mkdir -p "$DEST_DIR"

# 6. Navigate to the parent directory of SOURCE_DIR to avoid zipping the full absolute path
# Using dirname and basename keeps the zip file clean (just the folder itself)
PARENT_DIR=$(dirname "$SOURCE_DIR")
TARGET_FOLDER=$(basename "$SOURCE_DIR")

cd "$PARENT_DIR" || exit

# 7. Create the zip file (-r makes it recursive)
##echo "Zipping '$TARGET_FOLDER' into '$ZIP_FILENAME'..."
zip -rq "$ZIP_FILENAME" "$TARGET_FOLDER"

# 8. Move the zip file to the predefined location
##echo "Moving '$ZIP_FILENAME' to '$DEST_DIR'..."
mv "$ZIP_FILENAME" "$DEST_DIR/"

echo "Success! Archive saved at: ${DEST_DIR}/${ZIP_FILENAME}"
