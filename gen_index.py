"""포스트 제목의 ▷ 구분자를 파싱하여 index.qmd를 자동 생성.

분류 규칙 (categories frontmatter 불필요):
- title: "연구 ▷ HST ▷ ..." → main="연구", sub="HST"
- title: "공부 ▷ 토폴로지"   → main="공부", sub=없음
- title: "'메모 ▷ ..."       → main="메모" (따옴표 제거)

메인 탭 순서: 연구 > 리뷰 > 공부 > 책공부 > 나머지 가나다순
하위 카테고리가 있는 메인 탭만 중첩 panel-tabset 생성.
"""
import os, re, unicodedata, yaml
from collections import defaultdict
from pathlib import Path

POSTS_DIR = Path(__file__).parent / "Posts"
INDEX_PATH = Path(__file__).parent / "index.qmd"

MAIN_ORDER = ["연구", "리뷰", "공부", "책공부", "메모"]


def scan_posts():
    """제목에서 ▷ 파싱 → {main: {sub: count}}"""
    tree = defaultdict(lambda: defaultdict(int))

    for f in sorted(POSTS_DIR.iterdir()):
        if f.name.startswith(("_", "-")) or f.suffix not in (".md", ".qmd", ".ipynb"):
            continue
        text = f.read_text(encoding="utf-8")
        m = re.match(r"^---\n(.*?)\n---", text, re.DOTALL)
        if not m:
            continue
        try:
            fm = yaml.safe_load(m.group(1))
        except Exception:
            continue
        title = fm.get("title", "")
        if not title:
            continue
        title = unicodedata.normalize("NFC", str(title).strip().strip("'\""))

        parts = [p.strip() for p in title.split("▷")]
        if len(parts) < 2 or not parts[0]:
            continue

        main = parts[0]
        sub = parts[1] if len(parts) == 3 else None

        tree[main]["__total__"] += 1
        if sub:
            tree[main][sub] += 1

    return tree


def slugify(name):
    return re.sub(r"[^a-zA-Z0-9가-힣]", "", name).lower()


def build_index(tree):
    mains = sorted(tree.keys(), key=lambda k: (
        MAIN_ORDER.index(k) if k in MAIN_ORDER else len(MAIN_ORDER),
        k
    ))

    listings = []
    listings.append({
        "id": "all",
        "contents": "Posts",
        "sort": "date desc",
        "type": "table",
        "sort-ui": True,
        "filter-ui": True,
        "max-items": 999,
    })

    for main in mains:
        slug = slugify(main)
        subs = {k: v for k, v in tree[main].items() if k != "__total__"}

        listings.append({
            "id": f"{slug}-all",
            "contents": "Posts",
            "sort": "date desc",
            "type": "table",
            "include": {"title": f"*{main}*"},
            "max-items": 999,
        })

        for sub in sorted(subs.keys()):
            sub_slug = slugify(sub)
            listings.append({
                "id": f"{slug}-{sub_slug}",
                "contents": "Posts",
                "sort": "date desc",
                "type": "table",
                "include": {"title": f"*{main}*{sub}*"},
                "max-items": 999,
            })

    fm = yaml.dump(
        {"listing": listings, "page-layout": "full", "title-block-banner": False},
        default_flow_style=False,
        allow_unicode=True,
        sort_keys=False,
    )

    lines = []
    lines.append("::: {.panel-tabset}\n")
    lines.append("## 전체\n")
    lines.append("::: {#all}\n:::\n")

    for main in mains:
        slug = slugify(main)
        subs = {k: v for k, v in tree[main].items() if k != "__total__"}

        lines.append(f"## {main}\n")

        if subs:
            lines.append("::: {.panel-tabset}\n")
            lines.append("### 전체\n")
            lines.append(f"::: {{#{slug}-all}}\n:::\n")
            for sub in sorted(subs.keys()):
                sub_slug = slugify(sub)
                lines.append(f"### {sub}\n")
                lines.append(f"::: {{#{slug}-{sub_slug}}}\n:::\n")
            lines.append(":::\n")
        else:
            lines.append(f"::: {{#{slug}-all}}\n:::\n")

    lines.append(":::\n")

    content = f"---\n{fm}---\n\n" + "\n".join(lines)
    return content


if __name__ == "__main__":
    tree = scan_posts()
    print("감지된 카테고리:")
    for main in sorted(tree.keys()):
        subs = {k: v for k, v in tree[main].items() if k != "__total__"}
        total = tree[main]["__total__"]
        if subs:
            print(f"  {main} ({total}) → {', '.join(sorted(subs.keys()))}")
        else:
            print(f"  {main} ({total})")

    content = build_index(tree)
    INDEX_PATH.write_text(content, encoding="utf-8")
    print(f"\n{INDEX_PATH} 생성 완료")
