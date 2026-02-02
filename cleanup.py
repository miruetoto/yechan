#!/usr/bin/env python3
"""
블로그 관리 도구
1. 포스트 파일명 정리 (frontmatter 기반)
2. 미사용 이미지 정리 (삭제)
3. 이미지 파일명 정리 (파일명_NN 형식)
"""

import re
import json
import shutil
from pathlib import Path
from datetime import datetime
from urllib.parse import unquote, quote

# 색상 코드
class Colors:
    RED = '\033[0;31m'
    GREEN = '\033[0;32m'
    YELLOW = '\033[1;33m'
    BLUE = '\033[0;34m'
    NC = '\033[0m'

def print_color(text, color=Colors.NC):
    print(f"{color}{text}{Colors.NC}")

# ============================================================
# 포스트 파일명 정리 관련 함수
# ============================================================

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

def rename_posts(posts_dir: Path):
    """포스트 파일명을 frontmatter 기반으로 정리"""
    print_color("포스트 파일명 정리 중...", Colors.GREEN)
    print()

    # 모든 포스트 파일 검색
    patterns = ['*.md', '*.qmd', '*.ipynb']
    all_files = []
    for pattern in patterns:
        all_files.extend(posts_dir.glob(pattern))

    all_files = sorted(all_files)

    print(f"  총 {len(all_files)}개 파일 검사 중...")
    print()

    # 변경 계획 수립
    rename_plan = []  # [(old_path, new_name), ...]

    for file_path in all_files:
        new_name = get_new_filename(file_path)
        if new_name and new_name != file_path.name:
            rename_plan.append((file_path, new_name))

    if not rename_plan:
        print_color("✓ 변경할 파일이 없습니다.", Colors.GREEN)
        return

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
        return

    if confirm != 'y':
        print_color("취소되었습니다.", Colors.YELLOW)
        return

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

# ============================================================
# 이미지 정리 관련 함수
# ============================================================

def find_image_references_with_source(file_path: Path) -> dict:
    """파일에서 이미지 참조 찾기 (이미지 -> 소스파일 매핑)"""
    image_extensions = r'\.(png|jpg|jpeg|gif|webp|svg|PNG|JPG|JPEG|GIF|WEBP|SVG)'
    patterns = [
        rf'\!\[[^\]]*\]\(([^)]*{image_extensions})\)',
        rf'<img[^>]+src=["\']([^"\']*{image_extensions})["\']',
        rf'((?:\.?/)?attachments/[^\s\)\]"\']*{image_extensions})',
    ]

    references = {}
    try:
        content = file_path.read_text(encoding='utf-8', errors='ignore')
        for pattern in patterns:
            matches = re.findall(pattern, content, re.IGNORECASE)
            for match in matches:
                if isinstance(match, tuple):
                    match = match[0]
                filename = Path(match).name
                filename = unquote(filename)
                if filename:
                    references[filename] = file_path
    except Exception:
        pass

    return references

def find_image_references_ordered(file_path: Path) -> list:
    """파일에서 이미지 참조를 등장 순서대로 찾기 (중복 제거)"""
    image_pattern = r'\!\[[^\]]*\]\(([^)]*\.(?:png|jpg|jpeg|gif|webp|svg))\)'

    references = []
    seen = set()
    try:
        content = file_path.read_text(encoding='utf-8', errors='ignore')
        matches = re.findall(image_pattern, content, re.IGNORECASE)
        for match in matches:
            filename = Path(match).name
            filename = unquote(filename)
            if filename and filename not in seen:
                seen.add(filename)
                references.append(filename)
    except Exception:
        pass

    return references

def get_all_attachments(attachments_dir: Path) -> set:
    """attachments 폴더의 모든 이미지 파일 목록"""
    image_extensions = {'.png', '.jpg', '.jpeg', '.gif', '.webp', '.svg'}
    images = set()

    if not attachments_dir.exists():
        return images

    for file in attachments_dir.iterdir():
        if file.is_file() and file.suffix.lower() in image_extensions:
            images.add(file.name)

    return images

def cleanup_unused_images(posts_dir: Path, attachments_dir: Path, current_dir: Path):
    """미사용 이미지 정리"""
    print_color("1. MD/QMD 파일에서 참조된 이미지 수집 중...", Colors.GREEN)

    used_images = set()
    image_sources = {}

    for pattern in ['**/*.md', '**/*.qmd']:
        for file_path in posts_dir.glob(pattern):
            refs = find_image_references_with_source(file_path)
            used_images.update(refs.keys())
            image_sources.update(refs)

    index_qmd = current_dir / "index.qmd"
    if index_qmd.exists():
        refs = find_image_references_with_source(index_qmd)
        used_images.update(refs.keys())
        image_sources.update(refs)

    print(f"  ✓ 참조된 이미지: {len(used_images)}개")
    print()

    print_color("2. attachments 폴더 스캔 중...", Colors.GREEN)
    all_images = get_all_attachments(attachments_dir)
    print(f"  ✓ 전체 이미지: {len(all_images)}개")
    print()

    unused_images = all_images - used_images
    missing_images = used_images - all_images
    used_count = len(used_images & all_images)
    usage_rate = used_count * 100 // len(all_images) if all_images else 0

    print("=" * 50)
    print_color("분석 결과", Colors.BLUE)
    print("=" * 50)
    print(f"전체: {Colors.YELLOW}{len(all_images)}개{Colors.NC} | "
          f"사용: {Colors.GREEN}{used_count}개{Colors.NC} | "
          f"미사용: {Colors.RED}{len(unused_images)}개{Colors.NC} | "
          f"사용률: {usage_rate}%")

    if missing_images:
        print()
        print_color(f"⚠ 누락된 이미지 {len(missing_images)}개", Colors.YELLOW)
        for img in sorted(missing_images)[:5]:
            source = image_sources.get(img)
            if source:
                print(f"  ! {img}")
                print(f"    → {source.name}")
            else:
                print(f"  ! {img}")
        if len(missing_images) > 5:
            print(f"  ... 외 {len(missing_images) - 5}개")

    if unused_images:
        print()
        print_color("3. 미사용 이미지 삭제 중...", Colors.GREEN)

        sorted_unused = sorted(unused_images)
        deleted = 0
        for img in sorted_unused:
            file_path = attachments_dir / img
            try:
                file_path.unlink()
                print(f"  삭제: {img}")
                deleted += 1
            except Exception as e:
                print_color(f"  실패: {img} - {e}", Colors.RED)

        print()
        print_color(f"✓ {deleted}개 파일 삭제 완료", Colors.GREEN)
    else:
        print()
        print_color("✓ 모든 이미지가 사용 중입니다!", Colors.GREEN)

def rename_images_for_file(md_path: Path, attachments_dir: Path) -> list:
    """MD 파일의 이미지들을 파일명 기반으로 리네임"""
    # MD 파일명에서 확장자 제거
    base_name = md_path.stem  # e.g., "250910_책공부_The Book of Why"

    # 이미 정리된 파일명 패턴 (파일명_NN.ext)
    escaped_base = re.escape(base_name)
    already_renamed_pattern = re.compile(rf'^{escaped_base}_\d{{2}}\.\w+$')

    # 등장 순서대로 이미지 찾기
    images = find_image_references_ordered(md_path)
    if not images:
        return []

    # 이미 정리된 이미지는 제외
    images_to_rename = []
    for img in images:
        if not already_renamed_pattern.match(img):
            images_to_rename.append(img)

    if not images_to_rename:
        return []

    # 기존에 해당 파일명으로 된 파일들 확인 (번호 이어서 부여)
    existing_nums = []
    for f in attachments_dir.iterdir():
        m = re.match(rf'^{escaped_base}_(\d{{2}})\.\w+$', f.name)
        if m:
            existing_nums.append(int(m.group(1)))

    next_num = max(existing_nums) + 1 if existing_nums else 1

    # 리네임 계획 수립
    rename_plan = []  # [(old_name, new_name), ...]

    for old_name in images_to_rename:
        old_path = attachments_dir / old_name
        if not old_path.exists():
            continue

        ext = old_path.suffix.lower()
        new_name = f"{base_name}_{next_num:02d}{ext}"

        # 새 이름이 이미 존재하면 번호 증가
        while (attachments_dir / new_name).exists():
            next_num += 1
            new_name = f"{base_name}_{next_num:02d}{ext}"

        rename_plan.append((old_name, new_name))
        next_num += 1

    if not rename_plan:
        return []

    # MD 파일 내용 업데이트
    content = md_path.read_text(encoding='utf-8')
    for old_name, new_name in rename_plan:
        # URL 인코딩된 버전도 처리
        old_encoded = quote(old_name)
        new_encoded = quote(new_name)

        content = content.replace(f"attachments/{old_encoded}", f"attachments/{new_encoded}")
        content = content.replace(f"attachments/{old_name}", f"attachments/{new_name}")

    md_path.write_text(content, encoding='utf-8')

    # 실제 파일 리네임
    for old_name, new_name in rename_plan:
        old_path = attachments_dir / old_name
        new_path = attachments_dir / new_name
        if old_path.exists():
            shutil.move(str(old_path), str(new_path))

    return rename_plan

def rename_images(posts_dir: Path, attachments_dir: Path):
    """이미지 파일명을 MD 파일명 기반으로 정리"""
    print_color("이미지 파일명 정리 중...", Colors.GREEN)
    print()

    total_renamed = 0

    # 모든 MD/QMD 파일 처리
    all_files = list(posts_dir.glob('**/*.md')) + list(posts_dir.glob('**/*.qmd'))
    for md_path in sorted(all_files):
        renamed = rename_images_for_file(md_path, attachments_dir)
        if renamed:
            print(f"  {md_path.name}")
            for old_name, new_name in renamed:
                print(f"    {old_name} → {new_name}")
            total_renamed += len(renamed)

    if total_renamed > 0:
        print()
        print_color(f"✓ {total_renamed}개 이미지 파일명 변경 완료", Colors.GREEN)
    else:
        print_color("✓ 변경할 이미지가 없습니다 (이미 정리됨)", Colors.GREEN)

# ============================================================
# 메인 함수
# ============================================================

def main():
    print("=" * 50)
    print_color("블로그 관리 도구", Colors.BLUE)
    print("=" * 50)
    print()

    current_dir = Path.cwd()
    print_color(f"작업 디렉토리: {current_dir}", Colors.YELLOW)
    print()

    posts_dir = current_dir / "posts"
    attachments_dir = posts_dir / "attachments"

    if not posts_dir.exists():
        print_color("에러: posts 디렉토리를 찾을 수 없습니다.", Colors.RED)
        return 1

    # 모든 작업 실행
    rename_posts(posts_dir)
    print()
    print("-" * 50)
    print()
    if attachments_dir.exists():
        cleanup_unused_images(posts_dir, attachments_dir, current_dir)
        print()
        print("-" * 50)
        print()
        rename_images(posts_dir, attachments_dir)
    else:
        print_color("경고: posts/attachments 디렉토리가 없어 이미지 정리를 건너뜁니다.", Colors.YELLOW)

    print()
    print("=" * 50)
    print_color("완료", Colors.BLUE)
    print("=" * 50)

    return 0

if __name__ == "__main__":
    exit(main())
