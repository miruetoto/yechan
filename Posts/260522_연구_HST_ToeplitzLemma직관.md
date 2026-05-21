---
title: 연구 ▷ HST ▷ Toeplitz Lemma 직관 — 가중평균은 무게 따라간다
author: 유진
date: 05/22/2026
draft: false
output-file: 260522_cdd7e3.html
fontsize: 0.85em
---

ABC증명 블로그 Theorem~B 풀이에서 한 번 등장했던 외부 보조정리이다. 본 글에서는 그 정리가 **왜 그렇게 성립하는지** 만 직관 위주로 풀어쓴다.

**Toeplitz Lemma (Knopp 1928, §43).** 실수 배열 $\{a_{ns}\}_{n, s \geq 1} \subset \mathbb{R}$ 이 다음 세 조건을 만족한다 하자.

(T1) 어떤 상수 $c < \infty$ 가 존재해 $\sum_{s=1}^{\infty} |a_{ns}| \leq c$ 가 모든 $n \geq 1$ 에서 성립.
(T2) $\sum_{s=1}^{\infty} a_{ns} \to 1$ as $n \to \infty$.
(T3) 각 고정된 $s \geq 1$ 에 대해 $a_{ns} \to 0$ as $n \to \infty$.

그러면 임의의 수렴 수열 $x_s \to x \in \mathbb{R}$ 에 대해
$$
y_n \;:=\; \sum_{s=1}^{\infty} a_{ns}\, x_s \;\xrightarrow{n \to \infty}\; x.
$$

요약하면 **"수렴 수열의 가중평균은 그 극한으로 수렴"** — 단, 가중치 배열이 (T1)(T2)(T3) 을 만족하면.

한 줄로 직관을 잡으면 **가중평균은 무게가 실린 데로 따라간다.** (T3) 에 의해 무게가 점점 "이미 $x$ 에 가까워진 항들" 쪽으로만 실리니까 평균도 $x$ 로 간다. 세 조건은 각각 그 직관이 깨질 수 있는 시나리오를 막는 안전장치다. (T1) 은 가중치 절댓값 합이 어떤 상수 $c$ 안에 묶이게 해서 "$x$ 근처로 가도 가중치 합이 발산하면 곱해진 평균도 발산" 하는 일을 막는다. (T2) 는 가중치 총합이 1 로 수렴하게 해서 평균이 평균답게 작동하도록 한다 (없으면 $cx$ 같은 배수로 수렴할 수 있다). (T3) 은 어떤 고정된 위치 $s$ 의 가중치도 결국 0 으로 사라지게 해서 "초반에 $x$ 와 동떨어진 $x_1, x_2, \ldots$ 가 끝까지 무게 유지하며 평균을 그쪽으로 끌고 가는" 일을 막는다.

이 직관을 시험 성적 비유로 풀어보자. 학생 한 명의 시험 성적을 가중평균으로 평가한다고 하자. $x_s$ 가 $s$번째 시험 점수, $a_{ns}$ 가 학기 $n$ 시점에 $s$번째 시험에 부여되는 가중치, $y_n$ 이 학기 $n$ 시점의 평가 점수다. 학기말 실력이 $x$ 점으로 수렴 — 초반 시험은 들쭉날쭉했고 후반 시험은 다 $x$ 근처라고 하자. (T3) "옛날 시험 가중치는 학기 진행될수록 0 으로" 는 초반의 들쭉날쭉했던 점수가 학기말 평가에 거의 안 반영되게 한다. (T2) "가중치 총합은 1" 은 평균이 평균답게 작동하도록 보장한다 (1.5배·0.7배 평균으로 변질되지 않는다). (T1) "시험 가중치 절댓값 합이 상수 안에 묶임" 은 특정 시험 또는 그 묶음이 비정상적으로 부풀려지지 않게 한다. 결과적으로 학기말 평가는 후반 시험들 (이미 $x$ 근처) 의 가중평균이라 $x$ 로 수렴한다.

가장 익숙한 특수 케이스는 산술평균 (Cesàro) 이다. 가중치를 단순히
$$
a_{ns} \;=\; \begin{cases} 1/n, & s \leq n, \\ 0, & s > n, \end{cases}
$$
로 잡으면 $y_n = (x_1 + \cdots + x_n)/n$, 즉 그냥 산술평균이다. 세 조건은 모두 통과한다.
$$
\begin{aligned}
\text{(T1)} &\quad \sum_{s=1}^\infty |a_{ns}| = \sum_{s=1}^n \tfrac{1}{n} = 1 \leq 1
   && (\because c = 1) \\
\text{(T2)} &\quad \sum_{s=1}^\infty a_{ns} = 1 \to 1
   && (\because \text{고정 1}) \\
\text{(T3)} &\quad a_{ns} = \tfrac{1}{n} \to 0 \quad \text{각 고정 } s
   && (\because n \to \infty)
\end{aligned}
$$
그러므로 "$x_s \to x$ 이면 산술평균도 $x$" — 누구나 직관적으로 동의하는 명제다 (앞쪽 몇 개 이상한 값들은 $n$ 으로 나누면서 희석된다). Toeplitz Lemma 는 이 직관을 "가중치 모양만 (T1)(T2)(T3) 만족하면 다 OK" 로 일반화한 것이다.

증명도 직관 그대로 따라간다. 임의의 $\epsilon > 0$ 을 고정하자. $x_s \to x$ 이니 어떤 $S$ 가 존재해 $s \geq S$ 면 $|x_s - x| < \epsilon$ 이다. $y_n - x$ 를 다음과 같이 분해한다.
$$
\begin{aligned}
y_n - x
&= \sum_{s=1}^\infty a_{ns} x_s - x
   && (\because y_n \text{ 정의}) \\
&= \sum_{s=1}^\infty a_{ns}(x_s - x) + \Bigl(\sum_{s=1}^\infty a_{ns} - 1\Bigr) x
   && (\because x \sum_s a_{ns} \text{ 더하고 빼기}) \\
&= \underbrace{\sum_{s < S} a_{ns}(x_s - x)}_{=:\,T_1}
 + \underbrace{\sum_{s \geq S} a_{ns}(x_s - x)}_{=:\,T_2}
 + \underbrace{\Bigl(\sum_s a_{ns} - 1\Bigr) x}_{=:\,T_3}
\end{aligned}
$$
세 항을 각각 bound 한다.
$$
\begin{aligned}
|T_1|
&\leq \sum_{s < S} |a_{ns}| \cdot \max_{s < S} |x_s - x|
   && (\because \text{삼각부등식}) \\
\Rightarrow\ |T_1|
&\xrightarrow{n \to \infty} 0
   && (\because \text{(T3)},\ S \text{ 유한}) \\[2pt]
|T_2|
&\leq \epsilon \cdot \sum_{s \geq S} |a_{ns}|
   && (\because |x_s - x| < \epsilon \text{ for } s \geq S) \\
&\leq \epsilon \cdot c
   && (\because \text{(T1)}) \\[2pt]
|T_3|
&= \Bigl|\sum_s a_{ns} - 1\Bigr| \cdot |x|
   && (\because x \text{ 분리}) \\
\Rightarrow\ |T_3|
&\xrightarrow{n \to \infty} 0
   && (\because \text{(T2)})
\end{aligned}
$$
합치면 $\limsup_{n \to \infty} |y_n - x| \leq 0 + \epsilon\,c + 0 = \epsilon\,c$. $\epsilon > 0$ 임의이므로 $y_n \to x$. $\square$

증명에서 각 조건의 역할이 깔끔히 분리된다. (T3) 은 $T_1$ 처리 — 앞쪽 항의 무게가 사라진다. (T1) 은 $T_2$ 처리 — 뒤쪽 항이 $\epsilon$-가까워도 무한 합 가중치로 곱해지면 발산할 수 있는데 (T1) 이 그것을 막는다. (T2) 는 $T_3$ 처리 — 총합이 1 로 가야 진짜 평균이 된다.

마무리하면, **가중평균은 무게가 실린 데로 따라가고, 무게가 점점 $x$ 근처로만 실리니까 평균도 $x$.** 세 조건이 그 직관이 깨질 수 있는 세 시나리오 — (T1) 발산, (T2) 배수 변질, (T3) 옛 항 잔류 — 를 각각 막는 안전장치다. Cesàro 산술평균이 가장 익숙한 특수 케이스이며, Toeplitz Lemma 는 그 직관을 가중치 모양 일반으로 확장한다. ABC증명 Theorem~B 의 $SD^2_{ij}(t)/t^2$ 수렴 증명에서는 가중치 $\alpha_{ts} = 2s/[t(t+1)]$ 를 잡아 위 세 조건 검증 후 Toeplitz Lemma 한 번 적용하면 끝나는 짧은 증명이 된다.
