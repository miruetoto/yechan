#!/bin/bash
set -e

cd ~/Dropbox/01-rsch/9999-Yechan

# ============================================================
# 1단계 (로컬/맥): 백업 커밋 — 실패 시 복원 지점
# ============================================================
git add -A
git diff --cached --quiet || git commit -m "backup: before image processing"
BACKUP_COMMIT=$(git rev-parse HEAD)
echo "Backup commit: $BACKUP_COMMIT"

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

# ============================================================
# 2단계: Dropbox가 백업 커밋을 186으로 전파할 때까지 대기
# ============================================================
SSH_HOST="210.117.173.186"
REPO_PATH='~/Dropbox/01-rsch/9999-Yechan'

echo "Waiting for Dropbox sync to $SSH_HOST..."
for i in $(seq 1 30); do
  REMOTE_HEAD=$(ssh -o BatchMode=yes "$SSH_HOST" "cd $REPO_PATH && git rev-parse HEAD" 2>/dev/null || echo "")
  if [ "$REMOTE_HEAD" = "$BACKUP_COMMIT" ]; then
    echo "Sync confirmed after ${i}x2s"
    break
  fi
  sleep 2
done

if [ "$REMOTE_HEAD" != "$BACKUP_COMMIT" ]; then
  echo "ERROR: Dropbox sync not complete after 60s"
  echo "  local  : $BACKUP_COMMIT"
  echo "  186    : ${REMOTE_HEAD:-<none>}"
  exit 1
fi

# ============================================================
# 3단계 (원격/186): 이미지 정리 + render + commit + push
# ------------------------------------------------------------
# 왜 186에서?
#   - macOS APFS는 한글 파일명을 NFD로 강제 저장하는 반면
#     GitHub Pages(Linux)는 NFC로 요청 → 맥에서 만든 파일이 NFD로
#     push되면 배포 후 이미지 404.
#   - 리눅스에서 파일 생성·rename·commit·push를 하면 전 과정이 NFC.
# 인증: SSH 공개키(~/.ssh/id_ed25519). 스크립트에 비번 없음.
# ============================================================
ssh -o BatchMode=yes "$SSH_HOST" bash -l <<'REMOTE'
set -e
# bash -l (login shell)로 실행해야 .bash_profile/.profile/.bashrc 체인이 로드되어
# uv·quarto PATH가 잡힌다. non-interactive SSH는 기본적으로 rc 파일을 건너뜀.
# (fallback으로 .bashrc도 직접 source — guard로 early return하면 무시되니 무해)
[ -f ~/.bashrc ] && source ~/.bashrc

cd ~/Dropbox/01-rsch/9999-Yechan
[ -f .venv/bin/activate ] && source .venv/bin/activate

uv run python cleanup.py
quarto render

git add -A
git diff --cached --quiet || git commit -m "."
git push origin main
REMOTE

echo "Done — rendered and pushed from $SSH_HOST"

# 성공: 트랩 해제
trap - ERR INT TERM
