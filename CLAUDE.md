# 프로젝트 작업 내용

## 프로젝트 개요
- **프로젝트명**: Yechan Blog (신록예찬's Blog)
- **기술스택**: Quarto 기반 정적 블로그
- **배포**: GitHub Pages
- **저장소**: https://github.com/miruetoto/yechan

## 프로젝트 구조
```
yechan/
├── _quarto.yml          # Quarto 설정 파일
├── index.qmd           # 메인 페이지
├── posts/              # 소스 컨텐츠
│   ├── Yechan-md/      # 마크다운 포스트
│   │   ├── 첨부파일들/   # 이미지 파일들
│   │   └── *.md        # 블로그 포스트들
│   └── Yechan-ipynb/   # Jupyter 노트북
├── docs/               # 빌드 결과 (GitHub Pages 배포용)
├── filters/            # Quarto Lua 필터
│   └── fix-image-paths.lua
└── styles.css          # 커스텀 스타일
```

## 해결한 문제들

### 1. Git 충돌 문제
**문제**: 로컬과 리모트 브랜치가 diverged 상태
**해결**: `git push --force origin main`으로 강제 푸시

### 2. 이미지 표시 문제
**문제**: GitHub Pages에서 이미지가 깨져서 표시됨
- GitHub 마크다운에서는 정상 표시
- GitHub Pages에서만 이미지 깨짐

**원인 분석**:
- 마크다운 소스: `![](./첨부파일들/image.png)` (상대경로)
- 로컬 HTML: 절대경로 `file:///...` 생성
- GitHub Pages: 상대경로가 제대로 해석되지 않음

**해결 과정**:
1. **1차 시도**: `quarto render` 재실행
   - 상대경로는 올바르게 생성됨
   - 여전히 GitHub Pages에서 이미지 표시 안됨

2. **2차 시도**: Lua 필터를 이용한 자동 경로 변환
   - `filters/fix-image-paths.lua` 생성
   - 렌더링 시 이미지 경로를 GitHub raw URL로 자동 변환
   - 소스 파일은 변경하지 않고 HTML 출력만 수정

### 3. Lua 필터 구현
**파일**: `filters/fix-image-paths.lua`
```lua
function Image(elem)
  if elem.src:match("^%./첨부파일들/") then
    elem.src = elem.src:gsub("^%./첨부파일들/", "https://raw.githubusercontent.com/miruetoto/yechan/main/posts/Yechan-md/첨부파일들/")
  end
  return elem
end
```

**설정**: `_quarto.yml`에 필터 추가
```yaml
filters:
  - filters/fix-image-paths.lua
```

**결과**:
- 변환 전: `./첨부파일들/image.png`
- 변환 후: `https://raw.githubusercontent.com/miruetoto/yechan/main/posts/Yechan-md/첨부파일들/image.png`

### 4. 경로 수정
**문제**: 처음에 잘못된 경로 사용
- 잘못된 경로: `docs/posts/Yechan-md/첨부파일들/`
- 올바른 경로: `posts/Yechan-md/첨부파일들/`

**수정**: Lua 필터에서 경로 수정하여 GitHub raw URL이 정상 작동하도록 함

### 5. publish.command 스크립트 개선
**개선 사항**: 백업 및 동기화 프로세스 체계화
- 기존: Obsidian 파일 복사 후 바로 렌더링
- 개선: 백업 → 삭제 → 복사 → 렌더링 → 배포 순서로 체계화

**새로운 프로세스**:
1. **백업용 커밋**: 현재 상태를 백업 커밋으로 저장
2. **Yechan-md 내용 삭제**: 기존 파일들만 삭제 (폴더는 유지)
3. **Obsidian 파일 복사**: `cp -R` 명령으로 복사
4. **Quarto 렌더링**: 이미지 경로 자동 변환
5. **Git 커밋 및 푸시**: GitHub Pages에 배포

### 6. 이미지 파일 Git 추적 문제 해결
**문제**: 이미지 파일들이 Git에 추적되지 않아 GitHub Pages에서 404 오류
- 로컬에서는 정상 표시되지만 GitHub Pages에서는 이미지가 보이지 않음
- 120+ 개의 이미지 파일이 Git에 추가되지 않은 상태

**해결**:
- 모든 이미지 파일을 Git에 추가하여 추적 상태로 변경
- GitHub 저장소에 이미지 파일들 푸시 완료
- Lua 필터 URL 형식을 `raw.githubusercontent.com` 형식으로 수정

## 현재 상태
- ✅ 이미지 경로 문제 완전 해결
- ✅ GitHub Pages에서 이미지 정상 표시
- ✅ 소스 파일 변경 없이 자동화된 솔루션
- ✅ 향후 새로운 포스트 작성 시에도 자동으로 적용
- ✅ publish.command 스크립트로 백업/동기화/배포 자동화
- ✅ 모든 이미지 파일이 Git에 추적되어 GitHub Pages에서 접근 가능

## 사용법
1. 마크다운 파일에서 이미지 참조 시 기존대로 상대경로 사용
   ```markdown
   ![](./첨부파일들/image.png)
   ```

2. `quarto render` 실행하면 자동으로 GitHub raw URL로 변환
   ```html
   <img src="https://raw.githubusercontent.com/miruetoto/yechan/main/posts/Yechan-md/첨부파일들/image.png">
   ```

3. GitHub Pages에서 이미지 정상 표시

## 배포 명령어

### 수동 배포
```bash
quarto render
git add .
git commit -m "Update content"
git push origin main
```

### 자동 배포 (권장)
```bash
./publish.command
```
- 백업 커밋 → Obsidian 파일 동기화 → 렌더링 → 배포까지 자동 처리