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
{MD파일명}_NN.png
```
- {MD파일명}: 참조하는 MD 파일의 stem (확장자 제외한 전체 이름)
- NN: 해당 문서 내 등장 순서 (01, 02, ...)

> 같은 날짜에 여러 포스트가 올라갈 수 있으므로 날짜만으로 구분하지 않고 MD 파일명 전체를 prefix로 사용한다.

예: MD 파일이 `250819_책공부_확률변수의수렴.md` 이면
→ `250819_책공부_확률변수의수렴_01.png`, `250819_책공부_확률변수의수렴_02.png`

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
![](attachments/250819_책공부_확률변수의수렴_01.png)
```

**macOS 주의사항:** APFS는 한글 파일명을 NFD(자모 분해형)로 저장하므로, 스크립트는 내부적으로 NFC 정규화를 거쳐 MD 소스 문자열과 매칭한다.

## 배포

### 배포 흐름 (`pub.sh` 사용)
```bash
bash pub.sh
```

내부 동작:
1. **맥 로컬**: 백업 커밋 → `cleanup.py` → `quarto render` → 최종 커밋
2. **Dropbox**: `.git`을 210.117.173.182로 동기화 (맥↔서버 파일명 NFC 정규화는 Dropbox가 수행)
3. **원격 (182)**: `ssh`로 `git push origin main`만 실행
4. **실패 시**: 맥에서 자동으로 `git reset --hard <BACKUP_COMMIT>` + `git clean -fd`로 복원

### 왜 push만 원격에서 하나?
- **macOS APFS는 한글 파일명을 NFD(자모 분해형)로 저장**
- 맥에서 그대로 `git push`하면 NFD 바이트가 GitHub 리포에 박힘
- GitHub Pages(Linux)는 브라우저의 NFC 요청을 NFD 파일과 매칭 못 함 → 이미지 404
- Dropbox가 맥→Linux 서버 동기화 시 파일명을 **NFC로 정규화**해 보내므로, Linux(182)에서 push하면 리포 바이트가 NFC가 되어 정상 배포됨
- 전 과정을 원격에서 돌리지 않고 **push만** 원격에 두는 이유: git 조작을 한 호스트(맥)에서만 해야 Dropbox `.git` 양방향 충돌(`conflicted copy`)이 안 생김

### 주의사항
- 맥에서 직접 `git push origin main`은 위 이유로 금지. 항상 `pub.sh` 경유.
- 충돌 시 `git push --force`는 원격(182)에서 실행.
