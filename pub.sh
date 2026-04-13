#!/bin/bash
set -e

run_local() {
  cd ~/Dropbox/01-rsch/9999-Yechan
  source .venv/bin/activate
  uv run python url_encode_images.py
  uv run python cleanup.py
  quarto render
  git add -A
  git commit -m "."
  git push origin main
}

if [ "$(whoami)" = "root" ]; then
  # 186 서버(컨테이너) 내부에서 직접 실행
  run_local
else
  # 로컬(Mac)에서 실행 → 186으로 SSH
  sshpass -p '123qwe' ssh -tt -p 7749 root@210.117.173.186 << 'ENDSSH'
set -e
cd ~/Dropbox/01-rsch/9999-Yechan
source .venv/bin/activate
uv run python url_encode_images.py
uv run python cleanup.py
quarto render
git add -A
git commit -m "."
git push origin main
exit
ENDSSH
fi
