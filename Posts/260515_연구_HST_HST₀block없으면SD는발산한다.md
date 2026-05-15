---
title: "연구 ▷ HST ▷ HST₀: block 없으면 SD는 발산한다"
author: 신록예찬
date: 05/15/2026
draft: false
output-file: 260515_d3bd76.html
---

HST의 block/reset 메커니즘이 없으면 어떻게 될까? **HST₀**를 정의하자: 순수 랜덤워크로 이동하고, 도착한 노드에 무조건 눈을 적립. block도 reset도 없다.

Wheel $W_{21}$ (hub deg 20, leaf deg 3) 에서 uniform $\mu_0$, $b=0.05$, $\tau=10^8$ 으로 비교.

### 결과

![](attachments/260515_d3bd76_01.gif)

**HST** (Row 1):

`-` 전체 $SD^2 \propto t^{1.0}$ — Thm A 수렴. $SD^2/t \to c$.

`-` leaf 임베딩에서 ring 구조가 선명하게 드러남.

**HST₀** (Row 2):

`-` 전체 $SD^2 \propto t^{3.0}$ — **발산**. hub에 눈이 집중되어 drift regime과 동일한 $t^3$ 스케일링.

`-` leaf끼리만 보면 $SD^2 \propto t^{2.0}$ — 같은 degree인데도 $t$보다 빠르게 발산.

`-` 임베딩에서 구조 정보가 사라짐 (hub/leaf 분리만 남음).

### 해석

block 메커니즘은 단순히 "수렴 속도"를 높이는 게 아니라, **수렴 자체를 가능하게** 하는 핵심이다. block이 없으면:

`-` 워커가 한 지역에 머물면서 눈이 국소적으로 과적립

`-` 높이 차이가 시간에 비례하여 누적 → $SD^2$가 $t$보다 빠르게 증가

`-` 그래프의 세밀한 구조 (ring topology 등)가 drift에 묻혀 사라짐

HST의 "높은 곳에서 낮은 곳으로" flow + "갇히면 reset" 메커니즘이 국소 편향을 해소하고, 그래프-신호 상호작용을 SD에 기록하는 유일한 경로임을 보여준다.
