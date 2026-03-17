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
│   └── attachments/      # 이미지 파일 (YYMMDD_NN.png 형식)
├── docs/                 # 빌드 결과 (GitHub Pages 배포용)
├── image_cleanup.py      # 이미지 정리 도구
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
YYMMDD_NN.png
```
- YYMMDD: 참조하는 MD 파일의 날짜
- NN: 해당 문서 내 등장 순서 (01, 02, ...)

예: `250819_01.png`, `250819_02.png`

## 이미지 관리 도구

### image_cleanup.py
```bash
uv run python image_cleanup.py
```

**기능:**
1. **미사용 이미지 정리** - 참조되지 않는 이미지 자동 삭제
2. **파일명 정리** - 이미지 파일명을 `YYMMDD_NN.png` 형식으로 변환
3. **모두 실행** - 1 + 2 동시 실행

**이미지 참조 방식:**
```markdown
![](attachments/250819_01.png)
```

## 배포

### 수동 배포
```bash
quarto render
git add -A
git commit -m "Update"
git push origin main
```

### 충돌 시
```bash
git push --force origin main
```
