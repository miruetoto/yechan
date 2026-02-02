#!/usr/bin/env python3
"""
이미지 파일 정리 프로그램
1. posts/attachments/ 폴더의 미사용 이미지 자동 삭제
2. 이미지 파일명을 MD 파일명 기준으로 정리 (파일명_01.png 형식)
"""

import re
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

def get_date_prefix(md_filename: str) -> str:
    """MD 파일명에서 날짜 prefix 추출 (YYMMDD)"""
    match = re.match(r'^(\d{6})_', md_filename)
    if match:
        return match.group(1)
    return None

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

        # 리포트 생성
        print()
        print_color("4. 리포트 생성 중...", Colors.GREEN)

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        report_file = current_dir / f"이미지정리_리포트_{timestamp}.txt"

        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("=" * 50 + "\n")
            f.write("이미지 파일 정리 리포트\n")
            f.write(f"생성: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("=" * 50 + "\n\n")

            f.write("【요약】\n")
            f.write(f"  • 전체: {len(all_images)}개\n")
            f.write(f"  • 사용: {used_count}개\n")
            f.write(f"  • 삭제: {deleted}개\n")
            f.write(f"  • 사용률: {usage_rate}%\n\n")

            f.write("【삭제된 이미지】\n")
            for img in sorted_unused:
                f.write(f"  ✗ {img}\n")

            if missing_images:
                f.write("\n【누락된 이미지】\n")
                for img in sorted(missing_images):
                    f.write(f"  ! {img}\n")

        print_color(f"✓ 리포트: {report_file.name}", Colors.GREEN)
    else:
        print()
        print_color("✓ 모든 이미지가 사용 중입니다!", Colors.GREEN)

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

def main():
    print("=" * 50)
    print_color("이미지 파일 정리 프로그램", Colors.BLUE)
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

    if not attachments_dir.exists():
        print_color("에러: posts/attachments 디렉토리를 찾을 수 없습니다.", Colors.RED)
        return 1

    # 메뉴 선택
    print("작업 선택:")
    print("  1) 미사용 이미지 정리 (삭제)")
    print("  2) 이미지 파일명 정리 (파일명_NN 형식)")
    print("  3) 모두 실행 (1 + 2)")
    print()

    try:
        choice = input("선택 (1/2/3): ").strip()
    except (EOFError, KeyboardInterrupt):
        print()
        print_color("취소되었습니다.", Colors.YELLOW)
        return 0

    print()

    if choice == '1':
        cleanup_unused_images(posts_dir, attachments_dir, current_dir)
    elif choice == '2':
        rename_images(posts_dir, attachments_dir)
    elif choice == '3':
        cleanup_unused_images(posts_dir, attachments_dir, current_dir)
        print()
        print("-" * 50)
        print()
        rename_images(posts_dir, attachments_dir)
    else:
        print_color("잘못된 선택입니다.", Colors.RED)
        return 1

    print()
    print("=" * 50)
    print_color("완료", Colors.BLUE)
    print("=" * 50)

    return 0

if __name__ == "__main__":
    exit(main())
