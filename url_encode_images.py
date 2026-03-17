#!/usr/bin/env python3
"""
이미지 참조를 URL 인코딩 형식으로 자동 변환하는 스크립트
옵시디언 호환성을 위해 attachments/ 경로의 모든 이미지 참조를 URL 인코딩 형식으로 변환합니다.
"""

import re
import urllib.parse
from pathlib import Path


def url_encode_filename(match):
    """이미지 참조에서 파일명을 추출하여 URL 인코딩"""
    full_match = match.group(0)
    prefix = match.group(1)  # ![...](attachments/
    filename = match.group(2)  # 파일명

    # 이미 URL 인코딩된 경우 건너뛰기 (% 포함 여부로 판단)
    if '%' in filename:
        return full_match

    # 파일명만 URL 인코딩 (확장자 포함)
    encoded_filename = urllib.parse.quote(filename)

    return f"{prefix}{encoded_filename})"


def convert_markdown_file(file_path):
    """마크다운 파일의 이미지 참조를 URL 인코딩으로 변환"""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()

        # attachments/ 경로의 이미지 참조 패턴 매칭
        # ![...](attachments/파일명.확장자) 형식
        pattern = r'(!\[.*?\]\(attachments/)([^)]+)(\))'

        original_content = content
        content = re.sub(pattern, url_encode_filename, content)

        # 변경사항이 있는 경우에만 파일 저장
        if content != original_content:
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(content)
            return True
        return False

    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return False


def main():
    """Posts 디렉토리의 모든 마크다운 파일 처리"""
    posts_dir = Path(__file__).parent / 'Posts'

    if not posts_dir.exists():
        print(f"Posts directory not found: {posts_dir}")
        return

    converted_count = 0
    total_count = 0

    # .md 파일만 처리
    for md_file in posts_dir.glob('*.md'):
        total_count += 1
        if convert_markdown_file(md_file):
            converted_count += 1
            print(f"✓ Converted: {md_file.name}")

    print(f"\nProcessed {total_count} files, converted {converted_count} files")


if __name__ == '__main__':
    main()
