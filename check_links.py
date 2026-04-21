#!/usr/bin/env python3
"""
문서별 이미지 링크 상태 체크.
각 MD/QMD 파일이 참조하는 이미지가 실제로 attachments/에 존재하는지 확인.

사용: uv run python check_links.py [--verbose]
"""

import re
import sys
import argparse
import unicodedata
from pathlib import Path
from urllib.parse import unquote


class Colors:
    RED = '\033[0;31m'
    GREEN = '\033[0;32m'
    YELLOW = '\033[1;33m'
    BLUE = '\033[0;34m'
    GRAY = '\033[0;90m'
    NC = '\033[0m'


def nfc(s: str) -> str:
    return unicodedata.normalize('NFC', s)


IMAGE_EXT = r'\.(png|jpg|jpeg|gif|webp|svg|PNG|JPG|JPEG|GIF|WEBP|SVG)'
REF_PATTERNS = [
    rf'\!\[[^\]]*\]\(([^\n]+?{IMAGE_EXT})\)',
    rf'<img[^>]+src=["\']([^"\']*{IMAGE_EXT})["\']',
]


def find_refs_ordered(md_path: Path) -> list[str]:
    """MD 파일의 이미지 참조를 등장 순서대로 반환 (NFC, URL-decoded 파일명)."""
    refs = []
    seen = set()
    try:
        content = md_path.read_text(encoding='utf-8', errors='ignore')
    except Exception:
        return refs
    for pat in REF_PATTERNS:
        for m in re.finditer(pat, content, re.IGNORECASE):
            path = m.group(1)
            name = nfc(unquote(Path(path).name))
            if name and name not in seen:
                seen.add(name)
                refs.append(name)
    return refs


def collect_attachments(attachments_dir: Path) -> set[str]:
    """attachments/ 최상위의 모든 이미지 파일명 (NFC). unused_bck/ 제외."""
    if not attachments_dir.exists():
        return set()
    exts = {'.png', '.jpg', '.jpeg', '.gif', '.webp', '.svg'}
    return {
        nfc(f.name)
        for f in attachments_dir.iterdir()
        if f.is_file() and f.suffix.lower() in exts
    }


def check(posts_dir: Path, attachments_dir: Path, verbose: bool = False):
    files_on_disk = collect_attachments(attachments_dir)

    docs = sorted(list(posts_dir.glob('*.md')) + list(posts_dir.glob('*.qmd')))

    total_refs = 0
    total_broken = 0
    broken_docs = []  # [(doc, [broken_refs])]
    ok_docs = []

    for doc in docs:
        refs = find_refs_ordered(doc)
        if not refs:
            continue
        broken = [r for r in refs if r not in files_on_disk]
        total_refs += len(refs)
        total_broken += len(broken)
        if broken:
            broken_docs.append((doc, refs, broken))
        else:
            ok_docs.append((doc, refs))

    print(f"{Colors.BLUE}{'=' * 60}{Colors.NC}")
    print(f"{Colors.BLUE}문서별 이미지 링크 상태 체크{Colors.NC}")
    print(f"{Colors.BLUE}{'=' * 60}{Colors.NC}")
    print()
    print(f"  문서: {len(ok_docs) + len(broken_docs)}개 (이미지 참조 있음)")
    print(f"  {Colors.GREEN}✓ 정상{Colors.NC}: {len(ok_docs)}개")
    print(f"  {Colors.RED}✗ 깨진 링크 있음{Colors.NC}: {len(broken_docs)}개")
    print()
    print(f"  총 참조: {total_refs}개 | {Colors.RED}깨진 참조: {total_broken}개{Colors.NC}")
    print()

    if broken_docs:
        print(f"{Colors.RED}{'-' * 60}{Colors.NC}")
        print(f"{Colors.RED}깨진 링크 상세{Colors.NC}")
        print(f"{Colors.RED}{'-' * 60}{Colors.NC}")
        for doc, refs, broken in broken_docs:
            print()
            print(f"{Colors.YELLOW}{doc.name}{Colors.NC}")
            print(f"  참조 {len(refs)}개 중 {Colors.RED}{len(broken)}개 누락{Colors.NC}:")
            for b in broken:
                print(f"    {Colors.RED}✗{Colors.NC} {b}")

    if verbose and ok_docs:
        print()
        print(f"{Colors.GRAY}{'-' * 60}{Colors.NC}")
        print(f"{Colors.GRAY}정상 문서 ({len(ok_docs)}){Colors.NC}")
        print(f"{Colors.GRAY}{'-' * 60}{Colors.NC}")
        for doc, refs in ok_docs:
            print(f"  {Colors.GREEN}✓{Colors.NC} {doc.name} ({len(refs)} refs)")

    print()
    return 0 if total_broken == 0 else 1


def main():
    ap = argparse.ArgumentParser(description="문서별 이미지 링크 상태 체크")
    ap.add_argument('--verbose', '-v', action='store_true',
                    help='정상 문서 목록도 출력')
    args = ap.parse_args()

    root = Path.cwd()
    posts_dir = root / "Posts"
    attachments_dir = posts_dir / "attachments"

    if not posts_dir.exists():
        print(f"{Colors.RED}에러: Posts 디렉토리를 찾을 수 없습니다.{Colors.NC}")
        return 1

    return check(posts_dir, attachments_dir, verbose=args.verbose)


if __name__ == "__main__":
    sys.exit(main())
