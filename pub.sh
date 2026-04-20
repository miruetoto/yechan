#!/bin/bash
# ============================================================
# Yechan 블로그 배포 스크립트 (맥 로컬 실행)
# ------------------------------------------------------------
# 맥에서 render·commit까지 하고, push만 Linux 서버(182)로 위임.
#
# 왜 push만 원격?
#   macOS APFS는 한글 파일명을 NFD(자모 분해)로 저장한다.
#   맥에서 push하면 NFD 바이트가 GitHub 리포에 그대로 올라가고,
#   GitHub Pages(Linux)는 NFC로 요청하므로 URL↔파일 매칭 실패 → 404.
#   Dropbox가 맥→182 동기화 시 파일명을 NFC로 정규화하므로,
#   182에서 push하면 리포에 NFC 바이트가 기록됨.
# ============================================================
set -e

cd ~/Dropbox/01-rsch/9999-Yechan
[ -f .venv/bin/activate ] && source .venv/bin/activate

# 1단계: 백업 커밋 — 실패 시 복원 지점
git add -A
git diff --cached --quiet || git commit -m "backup: before image processing"
BACKUP_COMMIT=$(git rev-parse HEAD)
echo "Backup commit: $BACKUP_COMMIT"

# 에러/중단/종료 시 백업 시점으로 자동 복원
rollback_on_failure() {
  local exit_code=$?
  echo ""
  echo "================================================================"
  echo "pub.sh 실패 (exit $exit_code) — 백업 시점으로 자동 복원합니다"
  echo "  target: $BACKUP_COMMIT"
  echo "================================================================"
  git reset --hard "$BACKUP_COMMIT" || true
  git clean -fd || true
  echo "복원 완료. HEAD=$(git rev-parse --short HEAD)"
}
trap rollback_on_failure ERR INT TERM

# 2단계: 이미지 정리 + 블로그 렌더 (맥에서)
uv run python cleanup.py
quarto render

# 3단계: 최종 커밋 (맥에서)
git add -A
git diff --cached --quiet || git commit -m "."
LOCAL_HEAD=$(git rev-parse HEAD)
echo "Local HEAD: $LOCAL_HEAD"

# 4단계: Dropbox가 .git을 182로 동기화할 때까지 대기
SSH_HOST="210.117.173.182"
REPO_PATH='~/Dropbox/01-rsch/9999-Yechan'

echo "Waiting for Dropbox sync to $SSH_HOST..."
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
  echo "  182   : ${REMOTE_HEAD:-<none>}"
  echo "직접 182에서 'git push origin main' 수동 실행하세요."
  exit 1
fi

# 5단계: 182에서 push만 실행 (read + network)
ssh -o BatchMode=yes "$SSH_HOST" "cd $REPO_PATH && git push origin main"
echo "Pushed from $SSH_HOST"

# 성공: 복원 트랩 해제
trap - ERR INT TERM
