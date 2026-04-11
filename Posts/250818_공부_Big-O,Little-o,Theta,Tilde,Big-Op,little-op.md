---
title: 공부 ▷ Big-O, Little-o, Theta, Tilde, Big-Op, little-op
author: 신록예찬
date: 08/18/2025
draft: false
---

# $O$, $o$, $\Theta$, $\sim$ 

| 번호 | 문제 | 정답 |
|------|------|------|
| 1 | $n = O(n^2)$ | T |
| 2 | $n^2 = O(n)$ | F |
| 3 | $\log n = O(n)$ | T |
| 4 | $n^{0.5} = O(n)$ | T |
| 5 | $n^2 = O(n^2)$ | T |
| 6 | $n = o(n^2)$ | T |
| 7 | $n^2 = o(n)$ | F |
| 8 | $\log n = o(n)$ | T |
| 9 | $\sqrt{n} = o(n)$ | T |
| 10 | $n^{0.9} = o(n)$ | T |
| 11 | $n^{1.1} = o(n)$ | F |
| 12 | $1 = o(\log n)$ | T |
| 13 | $\log n = o(1)$ | F |
| 14 | $n^3 = \Theta(n^3)$ | T |
| 15 | $n^3 = \Theta(n^4)$ | F |
| 16 | $n\log n = \Theta(n\log n)$ | T |
| 17 | $n\log n = \Theta(n)$ | F |
| 18 | $n^2\log n = \Theta(n^2\log_{10} n)$ | T |
| 19 | $n^2 = \Theta(n^2)$ | T |
| 20 | $n^2 = \Theta(n^3)$ | F |
| 21 | $n^2 + O(n) \sim n^2$ | T |
| 22 | $n^2 + o(n^2) \sim n^2$ | T |
| 23 | $n^2 + n^2 \sim n^2$ | F |
| 24 | $\log_{10} n \sim \log n$ | F |
| 25 | $\sqrt{n^2+n} \sim n + 0.5$ | T |
| 26 | $\sqrt{n^2+n} \sim n$ | T |
| 27 | $n^3 + 2n^2 \sim n^3$ | T |
| 28 | $n^3 + 0.5n^3 \sim n^3$ | F |
| 29 | $n^2\log n \sim n^2\log_{10} n$ | F |
| 30 | $n^2(1 + 1/n) \sim n^2$ | T |
| 31 | $n^2(1 + 1/\log n) \sim n^2$ | T |
| 32 | $n^2(1 + \log n) \sim n^2$ | F |
| 33 | $\log^2 n = O(\log^3 n)$ | T |
| 34 | $\log^3 n = O(\log^2 n)$ | F |
| 35 | $\log^2 n = o(\log^3 n)$ | T |
| 36 | $\log^3 n = o(\log^2 n)$ | F |
| 37 | $n^4 = O(n^5)$ | T |
| 38 | $n^5 = O(n^4)$ | F |
| 39 | $n^4 = o(n^5)$ | T |
| 40 | $n^5 = o(n^4)$ | F |
| 41 | $2^n = O(3^n)$ | T |
| 42 | $3^n = O(2^n)$ | F |
| 43 | $2^n = o(3^n)$ | T |
| 44 | $3^n = o(2^n)$ | F |
| 45 | $n^k = o(n^{k+1})$ | T |
| 46 | $n^{k+1} = o(n^k)$ | F |
| 47 | $n\log^2 n = O(n^{1.1})$ | T |
| 48 | $n\log^2 n = o(n^{1.1})$ | T |
| 49 | $n^{1.1} = o(n\log^2 n)$ | F |
| 50 | $n\log^2 n \sim n\log^2 n$ | T |

# $O_p$, $o_p$ 

`-` 전제조건 

- $X_1, X_2, \dots$ : 독립·동일분포(i.i.d.) 확률변수    
- $\mathbb{E}[X_i] = \mu, \ \mathrm{Var}(X_i) = \sigma^2 < \infty$
- 표본평균: $\overline{X}_n = \frac1n\sum_{i=1}^n X_i$
- 표본합: $S_n = \sum_{i=1}^n X_i$
- 중심극한정리(CLT)와 약한대수의법칙(WLLN) 사용 가능

`-` 참고사항

- $\frac{X_n}{a_n} \overset{p}{\to} 0$ $\Leftrightarrow$ $X_n= o_P(a_n)$
- $X_n \overset{d}{\to} X$ $\Rightarrow$ $X_n=O_P(1)$ 

| 번호  | 문제                                                       | 정답  |
| --- | -------------------------------------------------------- | --- |
| 1   | $\overline{X}_n - \mu = O_p(n^{-1/2})$                   | T   |
| 2   | $\overline{X}_n - \mu = o_p(1)$                          | T   |
| 3   | $\overline{X}_n - \mu = o_p(n^{-1/2})$                   | F   |
| 4   | $\overline{X}_n = O_p(1)$                                | T   |
| 5   | $S_n = O_p(n)$                                           | T   |
| 6   | $S_n = o_p(n)$                                           | F   |
| 7   | $S_n/n = O_p(1)$                                         | T   |
| 8   | $S_n/n = o_p(1)$                                         | F   |
| 9   | $\sqrt{n}(\overline{X}_n - \mu) = O_p(1)$                | T   |
| 10  | $\sqrt{n}(\overline{X}_n - \mu) = o_p(1)$                | F   |
| 11  | $\overline{X}_n = o_p(n^{-1/2})$                         | F   |
| 12  | $\overline{X}_n - \mu = O_p(n^{-1})$                     | F   |
| 13  | $\overline{X}_n - \mu = O_p(1)$                          | T   |
| 14  | $\overline{X}_n - \mu = o_p(n^{-1})$                     | F   |
| 15  | $n(\overline{X}_n - \mu) = O_p(1)$                       | F   |
| 16  | $\sqrt{\log n}(\overline{X}_n - \mu) = O_p(1)$           | T   |
| 17  | $S_n / n^{1.1} = o_p(1)$                                 | T   |
| 18  | $S_n / n^{0.9} = o_p(1)$                                 | F   |
| 19  | $\overline{X}_n^2 = O_p(1)$                              | T   |
| 20  | $\overline{X}_n^2 = o_p(1)$                              | F   |
| 21  | $S_n^2 = O_p(n^2)$                                       | T   |
| 22  | $S_n^2 = o_p(n^2)$                                       | F   |
| 23  | $\overline{X}_n - \mu = O_p(n^{-0.4})$                   | T   |
| 24  | $\overline{X}_n - \mu = o_p(n^{-0.4})$                   | F   |
| 25  | $S_n / \sqrt{n} = O_p(1)$                                | T   |
| 26  | $S_n / \sqrt{n} = o_p(1)$                                | F   |
| 27  | $\overline{X}_n - \mu = o_p(1/\sqrt{\log n})$            | F   |
| 28  | $\overline{X}_n - \mu = O_p(1/\sqrt{\log n})$            | T   |
| 29  | $S_n = O_p(n^{0.9})$                                     | F   |
| 30  | $S_n = o_p(n^{1.1})$                                     | T   |
| 31  | $\overline{X}_n - \mu = O_p(n^{-a}) \ \forall a \le 1/2$ | T   |
| 32  | $\overline{X}_n - \mu = o_p(n^{-a}) \ \forall a < 1/2$   | T   |
| 33  | $\overline{X}_n - \mu = o_p(n^{-1/2})$                   | F   |
| 34  | $S_n / n^{1/2} = O_p(1)$                                 | T   |
| 35  | $S_n / n^{1/2} = o_p(1)$                                 | F   |
| 36  | $S_n / n^{1/2} = O_p(\sqrt{\log n})$                     | T   |
| 37  | $S_n / n^{1/2} = o_p(\sqrt{\log n})$                     | F   |
| 38  | $\overline{X}_n - \mu = O_p(n^{-0.49})$                  | T   |
| 39  | $\overline{X}_n - \mu = o_p(n^{-0.49})$                  | F   |
| 40  | $\overline{X}_n - \mu = O_p(n^{-0.6})$                   | F   |
| 41  | $S_n / n = O_p(1)$                                       | T   |
| 42  | $S_n / n = o_p(1)$                                       | F   |
| 43  | $S_n / n^{1.5} = o_p(1)$                                 | T   |
| 44  | $S_n / n^{1.5} = O_p(1)$                                 | T   |
| 45  | $S_n / n^{0.5} = o_p(1)$                                 | F   |
| 46  | $S_n / n^{0.5} = O_p(1)$                                 | T   |
| 47  | $\overline{X}_n^k = O_p(1)$ for fixed $k$                | T   |
| 48  | $\overline{X}_n^k = o_p(1)$ for fixed $k$                | F   |
| 49  | $S_n^k = O_p(n^k)$                                       | T   |
| 50  | $S_n^k = o_p(n^k)$                                       | F   |
| 51  | $S_n / \sqrt{n\log n} = O_p(1)$                          | T   |
| 52  | $S_n / \sqrt{n\log n} = o_p(1)$                          | F   |
| 53  | $\overline{X}_n - \mu = O_p(\log^{-1/2} n)$              | T   |
| 54  | $\overline{X}_n - \mu = o_p(\log^{-1/2} n)$              | F   |
| 55  | $\overline{X}_n = O_p(1/\sqrt{n})$                       | T   |
| 56  | $\overline{X}_n = o_p(1/\sqrt{n})$                       | F   |
| 57  | $S_n = O_p(n^{3/2})$                                     | F   |
| 58  | $S_n = o_p(n^{3/2})$                                     | T   |
| 59  | $\overline{X}_n - \mu = O_p(1/n)$                        | F   |
| 60  | $\overline{X}_n - \mu = o_p(1/n)$                        | F   |
| 61  | $S_n / n^a = O_p(1)$ for $a \ge 1$                       | T   |
| 62  | $S_n / n^a = o_p(1)$ for $a > 1$                         | T   |
| 63  | $S_n / n^a = o_p(1)$ for $a = 1$                         | F   |
| 64  | $S_n / n^{0.8} = o_p(1)$                                 | F   |
| 65  | $S_n / n^{1.2} = o_p(1)$                                 | T   |
| 66  | $S_n / n^{0.8} = O_p(1)$                                 | F   |
| 67  | $S_n / n^{1.2} = O_p(1)$                                 | T   |
| 68  | $\overline{X}_n = O_p(1/n^a)$ for $a \le 1/2$            | T   |
| 69  | $\overline{X}_n = o_p(1/n^a)$ for $a < 1/2$              | T   |
| 70  | $\overline{X}_n = o_p(1/n^{1/2})$                        | F   |
| 71  | $\sqrt{n}(\overline{X}_n - \mu) = O_p(1)$                | T   |
| 72  | $\sqrt{n}(\overline{X}_n - \mu) = o_p(1)$                | F   |
| 73  | $\sqrt{\log n}(\overline{X}_n - \mu) = O_p(1)$           | T   |
| 74  | $\sqrt{\log n}(\overline{X}_n - \mu) = o_p(1)$           | F   |
| 75  | $S_n / n^{0.5} = O_p(1)$                                 | T   |
| 76  | $S_n / n^{0.5} = o_p(1)$                                 | F   |
| 77  | $S_n / \sqrt{n\log\log n} = O_p(1)$                      | T   |
| 78  | $S_n / \sqrt{n\log\log n} = o_p(1)$                      | F   |
| 79  | $\overline{X}_n - \mu = O_p((\log n)^{-1/2})$            | T   |
| 80  | $\overline{X}_n - \mu = o_p((\log n)^{-1/2})$            | F   |
| 81  | $S_n / \sqrt{n} = O_p(1)$                                | T   |
| 82  | $S_n / \sqrt{n} = o_p(1)$                                | F   |
| 83  | $S_n / n^{1.5} = O_p(1)$                                 | T   |
| 84  | $S_n / n^{1.5} = o_p(1)$                                 | T   |
| 85  | $S_n / n^{0.4} = O_p(1)$                                 | F   |
| 86  | $S_n / n^{0.4} = o_p(1)$                                 | F   |
| 87  | $\overline{X}_n = O_p(1)$                                | T   |
| 88  | $\overline{X}_n = o_p(1)$                                | T   |
| 89  | $\overline{X}_n^2 = O_p(1)$                              | T   |
| 90  | $\overline{X}_n^2 = o_p(1)$                              | F   |
| 91  | $\overline{X}_n - \mu = O_p(1/\sqrt{n})$                 | T   |
| 92  | $\overline{X}_n - \mu = o_p(1/\sqrt{n})$                 | F   |
| 93  | $S_n = O_p(n)$                                           | T   |
| 94  | $S_n = o_p(n)$                                           | F   |
| 95  | $S_n / n = O_p(1)$                                       | T   |
| 96  | $S_n / n = o_p(1)$                                       | F   |
| 97  | $S_n / \sqrt{n} = O_p(1)$                                | T   |
| 98  | $S_n / \sqrt{n} = o_p(1)$                                | F   |
| 99  | $\overline{X}_n - \mu = O_p(1/n^{0.49})$                 | T   |
| 100 | $\overline{X}_n - \mu = o_p(1/n^{0.49})$                 | F   |

# 혼합형태

| 번호 | 문제 | 정답 |
|---|---|---|
| 1 | $S_n/\sqrt{n} + O_p(1) = O_p(1)$ | T |
| 2 | $S_n/n + o_p(1) = \mu + o_p(1)$ | T |
| 3 | $\overline X_n-\mu + O_p(n^{-1/2}) = O_p(n^{-1/2})$ | T |
| 4 | $(\overline X_n-\mu)^2 + o_p(n^{-1}) = O_p(n^{-1})$ | T |
| 5 | $S_n/n + O_p(n^{-1/2}) = \mu + O_p(n^{-1/2})$ | T |
| 6 | $S_n/n + o_p(n^{-1}) = \mu + o_p(n^{-1})$ | T |
| 7 | $\sqrt{n}(\overline X_n-\mu) + O_p(1) = O_p(1)$ | T |
| 8 | $(\overline X_n-\mu)^{-1} + O_p(n^{1/2}) = O_p(n^{1/2})$ | F |
| 9 | $S_n/n^2 + O_p(n^{-1}) = O_p(n^{-1})$ | T |
| 10 | $n(\overline X_n-\mu)^2 + o_p(1) = O_p(1)$ | T |
| 11 | $S_n/\sqrt{n} + o(1) = O_p(1)$ | T |
| 12 | $\overline X_n-\mu + O(n^{-1}) = o_p(1)$ | T |
| 13 | $S_n/n + O_p(1/\log n) = \mu + o_p(1)$ | T |
| 14 | $S_n/\sqrt{n\log n} + O_p(1) = O_p(1)$ | T |
| 15 | $\sqrt{n}(\overline X_n-\mu) + o_p(1) = O_p(1)$ | T |
| 16 | $(\overline X_n-\mu)^3 + O_p(n^{-3/2}) = o_p(n^{-1})$ | T |
| 17 | $S_n/n + o(1) = \mu + o_p(1)$ | T |
| 18 | $S_n/\sqrt{n} + O(n^{-1/2}) = O_p(1)$ | T |
| 19 | $(S_n-n\mu)/n + O_p(n^{-1/2}) = o_p(1)$ | T |
| 20 | $(S_n-n\mu)/\sqrt{n} + o(1) = O_p(1)$ | T |
| 21 | $n(\overline X_n-\mu)^2 + O_p(1/\log n) = O_p(1)$ | T |
| 22 | $n^{1/2}(\overline X_n-\mu)^2 + O_p(n^{-1/2}) = O_p(n^{-1/2})$ | T |
| 23 | $S_n/n + O_p(n^{-1}) = \mu + o_p(1)$ | T |
| 24 | $S_n/n + o(n^{-1}) = \mu + o_p(1)$ | T |
| 25 | $(\overline X_n-\mu)\log n + O_p(n^{-1/2}) = O_p(n^{-1/2}\log n)$ | T |
| 26 | $S_n/n^2 + o_p(n^{-1}) = O_p(n^{-1})$ | T |
| 27 | $(S_n-n\mu)^2/n + O_p(1) = O_p(1)$ | T |
| 28 | $(S_n-n\mu)^2/n + o(1) = O_p(1)$ | T |
| 29 | $S_n^{-1} + O_p(n^{-1}) = O_p(n^{-1})$ | F |
| 30 | $(\overline X_n-\mu)^{-2} + O_p(n) = O_p(n)$ | F |
| 31 | $\overline X_n-\mu + O_p(n^{-1/2}\log^{-1}n) = O_p(n^{-1/2})$ | T |
| 32 | $S_n/\sqrt{n\log\log n} + O_p(1) = O_p(1)$ | T |
| 33 | $S_n/\sqrt{n} + O(1/\log n) = O_p(1)$ | T |
| 34 | $S_n/n + O(1) = \mu + O(1)$ | F |
| 35 | $n(\overline X_n-\mu)^2 + o(1/\log n) = O_p(1)$ | T |
| 36 | $(\overline X_n-\mu)\sqrt{\log n} + O_p(n^{-1/2}) = O_p(n^{-1/2}\sqrt{\log n})$ | T |
| 37 | $S_n/\sqrt{n} + o_p(1) = O_p(1)$ | T |
| 38 | $\overline X_n-\mu + o(n^{-1/2}) = O_p(n^{-1/2})$ | T |
| 39 | $(S_n-n\mu)/n + o_p(n^{-1/2}) = o_p(1)$ | T |
| 40 | $S_n/n^2 + O(n^{-1}) = O_p(n^{-1})$ | T |
| 41 | $(\overline X_n-\mu)^2\log n + O_p(n^{-1}) = O_p(n^{-1}\log n)$ | T |
| 42 | $\sqrt{n}(\overline X_n-\mu) + O(1) = O_p(1)$ | T |
| 43 | $S_n/\sqrt{n} + O_p(n^{-1/2}) = O_p(1)$ | T |
| 44 | $S_n/n + o(1/\log n) = \mu + o_p(1/\log n)$ | T |
| 45 | $(\overline X_n-\mu)^3 + o(n^{-3/2}) = o_p(n^{-1})$ | T |
| 46 | $S_n/\sqrt{n\log n} + o_p(1) = O_p(1)$ | T |
| 47 | $S_n/n + O_p(n^{-3/4}) = \mu + O_p(n^{-3/4})$ | T |
| 48 | $S_n/n + o_p(n^{-3/4}) = \mu + o_p(n^{-3/4})$ | T |
| 49 | $(S_n-n\mu)/\sqrt{n} + O_p(1/\log n) = O_p(1)$ | T |
| 50 | $n^{1/2}(\overline X_n-\mu)^2 + O_p(n^{-1/2}) = O_p(n^{-1/2})$ | T |

