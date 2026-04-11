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

        # title에서 ▷가 있는 경우, 마지막 ▷까지 볼드체로
        # 예: title: 연구 ▷ HST ▷ 눈이... → title: **연구 ▷ HST ▷** 눈이...
        def add_bold(match):
            title_content = match.group(1)

            # 이미 볼드체가 있으면 건너뛰기
            if '**' in title_content:
                return match.group(0)

            # 마지막 ▷의 위치 찾기
            last_arrow_pos = title_content.rfind(' ▷ ')

            if last_arrow_pos != -1:
                # 마지막 ▷ 이후의 텍스트
                after_last_arrow = title_content[last_arrow_pos + 3:].strip()
                # 마지막 ▷까지 (포함)
                before_and_arrow = title_content[:last_arrow_pos + 3].strip()

                return f"title: **{before_and_arrow}** {after_last_arrow}"
            else:
                return match.group(0)

        content = re.sub(r'^title:\s*(.+)$', add_bold, content, flags=re.MULTILINE)

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
    print(f"✓ {updated_count}개 파일의 제목에 볼드체 적용 완료")
else:
    print("수정할 제목이 없습니다.")
