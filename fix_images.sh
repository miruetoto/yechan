#\!/bin/bash

# For 2000-01-10-(메모) DGX staion.md
file="posts/Yechan-md/2000-01-10-(메모) DGX staion.md"

# Create backup
cp "$file" "$file.backup"

# Read the file and replace images one by one
sed -e 's/첨부파일들\/이미지파일\.png/첨부파일들\/img_7.png/1' \
    -e 's/첨부파일들\/이미지파일\.png/첨부파일들\/img_8.png/1' \
    -e 's/첨부파일들\/이미지파일\.png/첨부파일들\/img_4.png/1' \
    -e 's/첨부파일들\/이미지파일\.png/첨부파일들\/img_14.png/1' \
    -e 's/첨부파일들\/이미지파일\.png/첨부파일들\/img_5.png/1' \
    -e 's/첨부파일들\/이미지파일\.png/첨부파일들\/img_11.png/1' \
    -e 's/첨부파일들\/이미지파일\.png/첨부파일들\/img_12.png/1' \
    -e 's/첨부파일들\/이미지파일\.png/첨부파일들\/img_15.png/1' \
    -e 's/첨부파일들\/이미지파일\.png/첨부파일들\/img_10.png/1' \
    -e 's/첨부파일들\/이미지파일\.png/첨부파일들\/img_6.png/1' \
    -e 's/첨부파일들\/이미지파일\.png/첨부파일들\/img_9.png/1' \
    -e 's/첨부파일들\/이미지파일\.png/첨부파일들\/img_13.png/1' \
    "$file.backup" > "$file"

rm "$file.backup"
