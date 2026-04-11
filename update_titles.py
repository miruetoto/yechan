#!/usr/bin/env python3
import os
import re
from pathlib import Path

# Posts 디렉토리 경로
posts_dir = Path("/root/Dropbox/01-rsch/9999-Yechan/Posts")

# 모든 .md와 .qmd 파일 찾기
files = list(posts_dir.glob("*.md")) + list(posts_dir.glob("*.qmd"))

# 변환할 패턴
patterns = [
    (r"title: \(공부\) (.+)", r"title: 공부 ▷ \1"),
    (r"title: \(리뷰\) (.+)", r"title: 리뷰 ▷ \1"),
    (r"title: \(리비전\) (.+)", r"title: 리비전 ▷ \1"),
    (r"title: \(책공부\) (.+)", r"title: 책공부 ▷ \1"),
    (r"title: \(연구\) (.+)", r"title: 연구 ▷ \1"),
    (r"title: \(메모\) (.+)", r"title: 메모 ▷ \1"),
    (r"title: \(코드\) (.+)", r"title: 코드 ▷ \1"),
    (r"title: '\(메모\) (.+)'", r"title: '메모 ▷ \1'"),
]

# 각 파일 처리
for file_path in files:
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()

        original_content = content

        # 각 패턴 적용
        for pattern, replacement in patterns:
            content = re.sub(pattern, replacement, content)

        # 변경사항이 있으면 파일 저장
        if content != original_content:
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(content)
            print(f"✓ Updated: {file_path.name}")
        else:
            print(f"- Skipped: {file_path.name} (already updated or no match)")

    except Exception as e:
        print(f"✗ Error processing {file_path.name}: {e}")

print("\n완료!")
