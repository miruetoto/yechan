#!/usr/bin/env python3
"""
포스트 파일명 정리 도구
- frontmatter의 title과 date를 기반으로 파일명 자동 변경
- 형식: YYMMDD_카테고리_제목.ext
"""

import re
import json
from pathlib import Path
from datetime import datetime

# 색상 코드
class Colors:
    RED = '\033[0;31m'
    GREEN = '\033[0;32m'
    YELLOW = '\033[1;33m'
    BLUE = '\033[0;34m'
    NC = '\033[0m'

def print_color(text, color=Colors.NC):
    print(f"{color}{text}{Colors.NC}")

def parse_frontmatter(content: str) -> dict:
    """YAML frontmatter 파싱"""
    match = re.match(r'^---\s*\n(.*?)\n---', content, re.DOTALL)
    if not match:
        return {}

    frontmatter = {}
    for line in match.group(1).split('\n'):
        if ':' in line:
            key, value = line.split(':', 1)
            key = key.strip()
            value = value.strip().strip('"').strip("'")
            frontmatter[key] = value

    return frontmatter

def parse_frontmatter_from_ipynb(file_path: Path) -> dict:
    """Jupyter notebook에서 frontmatter 파싱"""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            notebook = json.load(f)

        # 첫 번째 raw 또는 markdown 셀에서 frontmatter 찾기
        for cell in notebook.get('cells', []):
            source = ''.join(cell.get('source', []))
            if source.startswith('---'):
                return parse_frontmatter(source)
    except Exception:
        pass

    return {}

def parse_date(date_str: str) -> str:
    """날짜 문자열을 YYMMDD 형식으로 변환"""
    if not date_str:
        return None

    # 다양한 날짜 형식 처리
    formats = [
        "%m/%d/%Y",      # 12/23/2022
        "%Y-%m-%d",      # 2022-12-23
        "%Y/%m/%d",      # 2022/12/23
        "%d/%m/%Y",      # 23/12/2022
        "%B %d, %Y",     # December 23, 2022
        "%b %d, %Y",     # Dec 23, 2022
    ]

    for fmt in formats:
        try:
            dt = datetime.strptime(date_str, fmt)
            return dt.strftime("%y%m%d")
        except ValueError:
            continue

    # YYMMDD 형식인지 확인
    if re.match(r'^\d{6}$', date_str):
        return date_str

    return None

def parse_title(title: str) -> tuple:
    """제목에서 카테고리와 실제 제목 추출

    예: "**(리뷰) EbayesThresh**" -> ("리뷰", "EbayesThresh")
    예: "**(공부) 토폴로지**" -> ("공부", "토폴로지")
    """
    if not title:
        return None, None

    # ** 제거
    title = title.replace('**', '').strip()

    # (카테고리) 제목 패턴
    match = re.match(r'^\(([^)]+)\)\s*(.+)$', title)
    if match:
        category = match.group(1).strip()
        actual_title = match.group(2).strip()
        return category, actual_title

    # 카테고리 없는 경우
    return None, title

def sanitize_filename(name: str) -> str:
    """파일명에 사용할 수 없는 문자 제거"""
    # 특수문자 제거 (일부 허용)
    name = re.sub(r'[<>:"/\\|?*]', '', name)
    # 연속 공백 제거
    name = re.sub(r'\s+', ' ', name)
    return name.strip()

def get_new_filename(file_path: Path) -> str:
    """파일의 새 이름 계산"""
    ext = file_path.suffix

    # frontmatter 읽기
    if ext == '.ipynb':
        frontmatter = parse_frontmatter_from_ipynb(file_path)
    else:
        try:
            content = file_path.read_text(encoding='utf-8')
            frontmatter = parse_frontmatter(content)
        except Exception:
            return None

    if not frontmatter:
        return None

    # 날짜 파싱
    date_str = frontmatter.get('date', '')
    date_prefix = parse_date(date_str)
    if not date_prefix:
        return None

    # 제목 파싱
    title = frontmatter.get('title', '')
    category, actual_title = parse_title(title)
    if not actual_title:
        return None

    # 새 파일명 생성
    actual_title = sanitize_filename(actual_title)

    if category:
        new_name = f"{date_prefix}_{category}_{actual_title}{ext}"
    else:
        new_name = f"{date_prefix}_{actual_title}{ext}"

    return new_name

def main():
    print("=" * 50)
    print_color("포스트 파일명 정리 도구", Colors.BLUE)
    print("=" * 50)
    print()

    current_dir = Path.cwd()
    posts_dir = current_dir / "posts"

    if not posts_dir.exists():
        print_color("에러: posts 디렉토리를 찾을 수 없습니다.", Colors.RED)
        return 1

    # 모든 포스트 파일 검색
    patterns = ['*.md', '*.qmd', '*.ipynb']
    all_files = []
    for pattern in patterns:
        all_files.extend(posts_dir.glob(pattern))

    all_files = sorted(all_files)

    print_color(f"총 {len(all_files)}개 파일 검사 중...", Colors.GREEN)
    print()

    # 변경 계획 수립
    rename_plan = []  # [(old_path, new_name), ...]

    for file_path in all_files:
        new_name = get_new_filename(file_path)
        if new_name and new_name != file_path.name:
            rename_plan.append((file_path, new_name))

    if not rename_plan:
        print_color("✓ 변경할 파일이 없습니다.", Colors.GREEN)
        return 0

    # 변경 내용 출력
    print_color(f"변경 예정: {len(rename_plan)}개 파일", Colors.YELLOW)
    print("-" * 50)

    for old_path, new_name in rename_plan:
        print(f"  {old_path.name}")
        print(f"  → {new_name}")
        print()

    # 확인
    print("-" * 50)
    try:
        confirm = input("변경하시겠습니까? (y/N): ").strip().lower()
    except (EOFError, KeyboardInterrupt):
        print()
        print_color("취소되었습니다.", Colors.YELLOW)
        return 0

    if confirm != 'y':
        print_color("취소되었습니다.", Colors.YELLOW)
        return 0

    # 실제 변경
    print()
    renamed = 0
    for old_path, new_name in rename_plan:
        new_path = old_path.parent / new_name

        # 새 파일명이 이미 존재하는지 확인
        if new_path.exists():
            print_color(f"  스킵 (이미 존재): {new_name}", Colors.YELLOW)
            continue

        try:
            old_path.rename(new_path)
            print(f"  ✓ {old_path.name} → {new_name}")
            renamed += 1
        except Exception as e:
            print_color(f"  실패: {old_path.name} - {e}", Colors.RED)

    print()
    print_color(f"✓ {renamed}개 파일 이름 변경 완료", Colors.GREEN)

    print()
    print("=" * 50)
    print_color("완료", Colors.BLUE)
    print("=" * 50)

    return 0

if __name__ == "__main__":
    exit(main())
