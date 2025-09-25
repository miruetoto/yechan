#!/bin/bash

# Yechan Blog 자동 동기화 및 배포 스크립트
# macOS 전용

set -e  # 에러 발생 시 스크립트 중단

# 스크립트가 있는 디렉토리로 이동
cd "$(dirname "$0")"

# macOS 체크
if [[ "$OSTYPE" != "darwin"* ]]; then
    echo "Error: This script is designed to run on macOS only."
    echo "Current OS: $OSTYPE"
    exit 1
fi

echo "Yechan Blog 동기화 및 배포 시작..."

# 1. 백업용 커밋
echo "현재 상태 백업 커밋 중..."
if [ -n "$(git status --porcelain)" ]; then
    git add .
    BACKUP_MSG="Backup before sync - $(date '+%Y-%m-%d %H:%M:%S')"
    git commit -m "$BACKUP_MSG"
    echo "백업 커밋 완료"
else
    echo "변경사항이 없어 백업 커밋을 건너뜁니다."
fi

# 2. Yechan-md 내용 삭제
echo "기존 Yechan-md 내용 삭제 중..."
TARGET_PATH="./posts/Yechan-md"
if [ -d "$TARGET_PATH" ]; then
    rm -rf "$TARGET_PATH"/*
    echo "기존 내용 삭제 완료"
fi

# 3. Obsidian 파일들 복사
echo "Obsidian 파일 복사 중..."
OBSIDIAN_PATH="/Users/cgb/Library/Mobile Documents/iCloud~md~obsidian/Documents/TopoNotes/Yechan"

if [ -d "$OBSIDIAN_PATH" ]; then
    # indexes 폴더를 제외하고 모든 파일 복사
    for item in "$OBSIDIAN_PATH"/*; do
        if [ "$(basename "$item")" != "indexes" ]; then
            cp -R "$item" "$TARGET_PATH"/
        fi
    done
    echo "파일 복사 완료 (indexes 폴더 제외)"
else
    echo "Error: Obsidian 폴더를 찾을 수 없습니다: $OBSIDIAN_PATH"
    exit 1
fi

# 4. Quarto 렌더링
echo "Quarto 렌더링 중..."
quarto render

echo "렌더링 완료"

# 5. Git 커밋 및 푸시
echo "Git 커밋 및 푸시 중..."

# 변경사항 확인
if [ -n "$(git status --porcelain)" ]; then
    # 변경사항이 있는 경우
    git add .
    
    # 현재 시간으로 커밋 메시지 생성
    COMMIT_MSG="Auto sync from Obsidian - $(date '+%Y-%m-%d %H:%M:%S')"
    
    git commit -m "$COMMIT_MSG

Generated with [Claude Code](https://claude.ai/code)

Co-Authored-By: Claude <noreply@anthropic.com>"
    
    git push origin main
    
    echo "블로그 업데이트 완료!"
    echo "변경사항이 GitHub Pages에 반영되기까지 몇 분 소요될 수 있습니다."
else
    echo "변경사항이 없습니다. 푸시를 건너뜁니다."
fi

echo "모든 작업이 완료되었습니다!"