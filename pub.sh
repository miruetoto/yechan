#!/bin/bash
set -e

cd ~/Dropbox/01-rsch/9999-Yechan
source .venv/bin/activate

# 이미지 처리 스크립트 실행 전 현재 상태 백업 커밋
# (url_encode_images.py 또는 cleanup.py가 파일 변경/삭제를 일으켜도
#  이 시점으로 git reset --hard HEAD~1 하면 복구 가능)
git add -A
git diff --cached --quiet || git commit -m "backup: before image processing"

uv run python url_encode_images.py
uv run python cleanup.py
quarto render

git add -A
git diff --cached --quiet || git commit -m "."
git push origin main
