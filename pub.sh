#!/bin/bash
# ============================================================
# Yechan 블로그 배포 스크립트 (186 서버 전용)
# ------------------------------------------------------------
# 실행 위치: 210.117.173.186 (Linux)
# 실행 방법: 사용자가 ssh로 접속해서 `bash pub.sh` 직접 실행
#
# 왜 맥이 아닌 186에서 돌리나?
#   - macOS APFS는 한글 파일명을 NFD(자모 분해)로 저장하는데
#     GitHub Pages(Linux)는 NFC로 요청함 → 맥에서 push하면 이미지 404
#   - 리눅스에서 파일 생성·cleanup·render·commit·push를 전부 하면 NFC
#
# 왜 ssh 자동화를 안 쓰나?
#   - 맥에서 ssh로 원격 실행하면 양쪽 git 조작으로 Dropbox `.git`
#     동기화 충돌(`conflicted copy`)이 발생함. 전부 한 호스트에서.
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

# 2단계: 이미지 정리 + 블로그 렌더
uv run python cleanup.py
quarto render

# 3단계: 최종 커밋 + push
git add -A
git diff --cached --quiet || git commit -m "."
git push origin main

echo "Done — rendered and pushed"

# 성공: 복원 트랩 해제
trap - ERR INT TERM
