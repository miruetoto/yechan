#!/usr/bin/env python3
import re
from pathlib import Path

# Posts 디렉토리 경로
posts_dir = Path("/root/Dropbox/01-rsch/9999-Yechan/Posts")

# 모든 .md와 .qmd 파일 찾기
files = list(posts_dir.glob("*.md")) + list(posts_dir.glob("*.qmd"))

updated_count = 0

for file_path in files:
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()

        original_content = content

        # title 필드에서 " ▷ "를 "▷"로 변환
        content = re.sub(r'(title:\s*[^\n]*?)\s+▷\s+', r'\1▷', content)

        # 변경사항이 있으면 파일 저장
        if content != original_content:
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(content)
            print(f"✓ {file_path.name}")
            updated_count += 1

    except Exception as e:
        print(f"✗ Error: {file_path.name} - {e}")

print()
if updated_count > 0:
    print(f"✓ {updated_count}개 파일의 제목 간격 수정 완료")
else:
    print("수정할 제목이 없습니다.")
