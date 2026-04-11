#!/bin/bash
set -e

# 원격 서버에서 실행
sshpass -p '123qwe' ssh -p 43052 root@210.117.173.186 << 'ENDSSH'
cd ~/Dropbox/01-rsch/9999-Yechan

# 이미지 참조를 URL 인코딩 형식으로 변환 (옵시디언 호환성)
uv run python url_encode_images.py

# 이미지 정리
uv run python cleanup.py

# Quarto 렌더링
quarto render

# Git 커밋 및 푸시
git add -A
git commit -m "."
git push origin main
ENDSSH
