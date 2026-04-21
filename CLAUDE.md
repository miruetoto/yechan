# Yechan Blog

## 프로젝트 개요
- **프로젝트명**: 신록예찬's Blog
- **기술스택**: Quarto 기반 정적 블로그
- **배포**: GitHub Pages
- **저장소**: https://github.com/miruetoto/yechan

## 프로젝트 구조
```
yechan/
├── _quarto.yml           # Quarto 설정
├── index.qmd             # 메인 페이지
├── Posts/                # 블로그 포스트
│   ├── *.md              # 마크다운 포스트 (YYMMDD_카테고리_제목.md)
│   ├── *.qmd             # Quarto 문서
│   ├── *.ipynb           # Jupyter 노트북
│   └── attachments/      # 이미지 파일 ({MD파일명}_NN.png 형식)
├── docs/                 # 빌드 결과 (GitHub Pages 배포용)
├── cleanup.py            # 블로그 관리 도구 (포스트/이미지 정리)
├── pub.sh                # 배포 스크립트 (render + commit + push)
└── styles.css            # 커스텀 스타일
```

## 파일 명명 규칙

### 포스트 파일
```
YYMMDD_카테고리_제목.md
```
예: `250819_책공부_확률변수의 수렴.md`

### 이미지 파일
```
{YYMMDD}_{hash6}_NN.ext
```
- `YYMMDD`: MD 파일 stem 앞 6자리 (날짜)
- `hash6`: MD 파일 stem의 SHA-256 앞 6자 (hex) — 같은 날 여러 포스트도 충돌 없이 구분
- `NN`: 해당 문서 내 등장 순서 (01, 02, ...)

> **ASCII 전용**. Unicode 정규화(NFC/NFD) 이슈를 원천 차단하기 위해 한글 파일명을 쓰지 않는다.
> 파일명만 봐도 동일 MD 소속 이미지인지 바로 식별 가능 (`hash6` 공유).

예: MD 파일이 `250819_책공부_확률변수의수렴.md` 이면
→ stem 해싱 → `a1b2c3` (예시)
→ `250819_a1b2c3_01.png`, `250819_a1b2c3_02.png`, ...

## 블로그 관리 도구

### cleanup.py
```bash
uv run python cleanup.py
```

**기능:**
1. **포스트 파일명 정리** - frontmatter(date, title) 기반으로 `YYMMDD_카테고리_제목.md` 형식으로 변환
2. **미사용 이미지 삭제** - 어떤 포스트에서도 참조되지 않는 이미지 자동 삭제
3. **이미지 파일명 정리** - 각 MD 파일이 참조하는 이미지를 등장 순서대로 `{MD파일명}_NN.png` 형식으로 리네임 + MD 내 참조 경로도 함께 업데이트

**이미지 참조 방식:**
```markdown
![](attachments/250819_a1b2c3_01.png)
```

**macOS 주의사항:** APFS는 한글 파일명을 NFD(자모 분해형)로 저장하므로, 스크립트는 내부적으로 NFC 정규화를 거쳐 MD 소스 문자열과 매칭한다.

## 배포

### 배포 흐름 (`pub.sh` 사용)
```bash
bash pub.sh
```

내부 동작:
1. **precomposeunicode 가드**: `git config core.precomposeunicode`가 true인지 확인 (아니면 중단)
2. **백업 커밋** → `cleanup.py` → `quarto render` → 최종 커밋 → `git push origin main`
3. **실패 시**: `git reset --hard <BACKUP_COMMIT>` + `git clean -fd`로 자동 복원

### 한글 파일명 NFC 처리
- **macOS APFS는 한글 파일명을 NFD(자모 분해형)로 저장**
- 맥에서 그대로 `git push`하면 NFD 바이트가 GitHub 리포에 박힘
- GitHub Pages(Linux)는 브라우저의 NFC 요청을 NFD 파일과 매칭 못 함 → 이미지 404
- 해결책: 이 리포에 `git config core.precomposeunicode true` 설정. 맥 git이 `git add`/commit 시점에 NFD→NFC 자동 정규화하여 커밋 오브젝트의 바이트가 NFC로 기록됨.

### 주의사항
- `core.precomposeunicode`가 해제되면 다시 NFD 오염 위험. `pub.sh`가 진입 시점에 체크함.
- 맥에서 직접 `git push origin main` 해도 OK (precomposeunicode 덕분). 다만 `pub.sh` 경유해야 백업/정리/렌더까지 일관.
