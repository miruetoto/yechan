#!/bin/bash
# ============================================================
# Yechan 블로그 배포 스크립트 (맥 로컬 실행)
# ------------------------------------------------------------
# 맥에서 render·commit·push 전부 수행.
#
# 한글 파일명 NFC 처리:
#   macOS APFS는 파일명을 NFD(자모 분해)로 저장. 그대로 push하면
#   GitHub Pages(Linux)의 NFC 요청과 매칭 실패 → 404.
#   해결: 이 리포에 `git config core.precomposeunicode true` 설정.
#   이러면 맥 git이 `git add`/commit 시점에 NFD→NFC 자동 정규화하여
#   커밋 오브젝트의 파일명 바이트가 NFC로 기록된다.
# ============================================================
set -e

cd ~/Dropbox/01-rsch/9999-Yechan
[ -f .venv/bin/activate ] && source .venv/bin/activate

# 선제 가드: precomposeunicode가 true가 아니면 NFD로 커밋될 위험
PRECOMPOSE=$(git config --get core.precomposeunicode || echo "")
if [ "$PRECOMPOSE" != "true" ]; then
  echo "ERROR: core.precomposeunicode is not 'true' (got: '${PRECOMPOSE:-<unset>}')"
  echo "실행: git config core.precomposeunicode true"
  echo "이 설정 없이 push하면 한글 파일명이 NFD로 기록되어 GitHub Pages에서 404."
  exit 1
fi

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

# 2단계: 이미지 정리 + index.qmd 자동생성 + 블로그 렌더 (맥에서)
uv run python cleanup.py
uv run python gen_index.py
rm -rf docs
quarto render
uv run python cleanup.py --postrender

# 3단계: 최종 커밋 (맥에서)
git add -A
git diff --cached --quiet || git commit -m "."
LOCAL_HEAD=$(git rev-parse HEAD)
echo "Local HEAD: $LOCAL_HEAD"

# 4단계: push
git push origin main
echo "Pushed from local."

# 성공: 복원 트랩 해제
trap - ERR INT TERM
