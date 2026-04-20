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

# push는 186 서버(Linux)에서 실행한다.
# 이유: macOS APFS는 한글 파일명을 NFD(자모 분해)로 저장하는 반면
#       GitHub Pages(Linux)는 NFC로 요청하므로, 맥에서 직접 push하면
#       리포 바이트가 NFD가 되어 배포 시 이미지·경로가 404로 깨진다.
#       Dropbox가 맥↔186 동기화 과정에서 파일명을 NFC로 정규화하므로
#       186에서 push하면 리포에 NFC 바이트가 올라간다.
# SSH 인증은 키 기반(~/.ssh/id_ed25519). 스크립트에 비밀번호 없음.
SSH_HOST="210.117.173.186"
REPO_PATH="~/Dropbox/01-rsch/9999-Yechan"
LOCAL_HEAD=$(git rev-parse HEAD)

echo "Local HEAD: $LOCAL_HEAD"
echo "Waiting for Dropbox sync to $SSH_HOST..."

# 186의 HEAD가 로컬과 일치할 때까지 폴링 (최대 ~60초)
for i in $(seq 1 30); do
  REMOTE_HEAD=$(ssh -o BatchMode=yes "$SSH_HOST" "cd $REPO_PATH && git rev-parse HEAD" 2>/dev/null || echo "")
  if [ "$REMOTE_HEAD" = "$LOCAL_HEAD" ]; then
    echo "Sync confirmed after ${i}x2s"
    break
  fi
  sleep 2
done

if [ "$REMOTE_HEAD" != "$LOCAL_HEAD" ]; then
  echo "ERROR: Dropbox sync not complete after 60s"
  echo "  local : $LOCAL_HEAD"
  echo "  186   : ${REMOTE_HEAD:-<none>}"
  echo "직접 186 터미널에서 'git push origin main'을 실행해 주세요."
  exit 1
fi

# 186에서 push
ssh -o BatchMode=yes "$SSH_HOST" "cd $REPO_PATH && git push origin main"
echo "Pushed from $SSH_HOST"
