---
title: (공부) Big-O, Little-o, Theta, Tilde, Big-Op, little-op
author: 신록예찬
date: 08/14/2025
draft: false
---
$$X_i \overset{iid}{\sim} (\mu, \sigma^2)$$

| 번호  | 문제                                             | 정답  |
| --- | ---------------------------------------------- | --- |
| 1   | $n = O(n^2)$                                   | T   |
| 2   | $n^2 = O(n)$                                   | F   |
| 3   | $\log n = O(n)$                                | T   |
| 4   | $n^{0.5} = O(n)$                               | T   |
| 5   | $n^2 = O(n^2)$                                 | T   |
| 6   | $n = o(n^2)$                                   | T   |
| 7   | $n^2 = o(n)$                                   | F   |
| 8   | $\log n = o(n)$                                | T   |
| 9   | $\sqrt{n} = o(n)$                              | T   |
| 10  | $n^{0.9} = o(n)$                               | T   |
| 11  | $n^{1.1} = o(n)$                               | F   |
| 12  | $1 = o(\log n)$                                | T   |
| 13  | $\log n = o(1)$                                | F   |
| 14  | $n = O(n)$                                     | T   |
| 15  | $n^3 = O(n^4)$                                 | T   |
| 16  | $n^4 = O(n^3)$                                 | F   |
| 17  | $n^2 = O(n^5)$                                 | T   |
| 18  | $n^5 = O(n^2)$                                 | F   |
| 19  | $n\log n = O(n^2)$                             | T   |
| 20  | $n^2\log n = O(n^2)$                           | F   |
| 21  | $n^2 = o(n^3)$                                 | T   |
| 22  | $n^3 = o(n^2)$                                 | F   |
| 23  | $\log^2 n = o(n^\epsilon)$                     | T   |
| 24  | $n^\epsilon = o(\log n)$                       | F   |
| 25  | $n^{0.5} = o(n^{0.6})$                         | T   |
| 26  | $n^{0.6} = o(n^{0.5})$                         | F   |
| 27  | $n^{0.1} = o(1)$                               | T   |
| 28  | $1 = o(n^{-0.1})$                              | F   |
| 29  | $n^2 + 3n = O(n^3)$                            | T   |
| 30  | $n^3 + 3n^2 = o(n^4)$                          | T   |
| 31  | $n^3 + 3n^2 = o(n^3)$                          | F   |
| 32  | $n^2 = O(n^2\log n)$                           | T   |
| 33  | $n^2\log n = O(n^2)$                           | F   |
| 34  | $\log n = O(\sqrt{n})$                         | T   |
| 35  | $\sqrt{n} = O(\log n)$                         | F   |
| 36  | $n = O(2^n)$                                   | T   |
| 37  | $2^n = O(n^k)$                                 | F   |
| 38  | $n^k = o(2^n)$                                 | T   |
| 39  | $2^n = o(3^n)$                                 | T   |
| 40  | $3^n = o(2^n)$                                 | F   |
| 41  | $n^2 + 3n = \Theta(n^2)$                       | T   |
| 42  | $3n^2 = \Theta(n^2)$                           | T   |
| 43  | $n^2 = \Theta(n)$                              | F   |
| 44  | $\log n = \Theta(\log_{10} n)$                 | T   |
| 45  | $n\log n = \Theta(n)$                          | F   |
| 46  | $n^{0.5} = \Theta(\sqrt{n})$                   | T   |
| 47  | $5n + 3\log n = \Theta(n)$                     | T   |
| 48  | $n = \Theta(n^2)$                              | F   |
| 49  | $\log n = \Theta(\sqrt{\log n})$               | F   |
| 50  | $n^{0.5} = \Theta(n^{0.6})$                    | F   |
| 51  | $n\log n = \Theta(n\log_{10} n)$               | T   |
| 52  | $n\log n = \Theta(n\log^2 n)$                  | F   |
| 53  | $\sqrt{n} + \log n = \Theta(\sqrt{n})$         | T   |
| 54  | $n^3 + n^2 + 5 = \Theta(n^3)$                  | T   |
| 55  | $n^{1.5} + n = \Theta(n^2)$                    | F   |
| 56  | $n^2 = \Theta(n^2 + n)$                        | T   |
| 57  | $n^2 + n = \Theta(n^2 + \sqrt{n})$             | T   |
| 58  | $n^2\log n = \Theta(n^2)$                      | F   |
| 59  | $n^2\log n = \Theta(n^2\log_{10} n)$           | T   |
| 60  | $n^2 + 3n + 1 = \Theta(n^3)$                   | F   |
| 61  | $n^2 + n\log n = \Theta(n^2)$                  | T   |
| 62  | $n^5 + n^4 = \Theta(n^5)$                      | T   |
| 63  | $n^4 + n^5 = \Theta(n^4)$                      | F   |
| 64  | $n^4 + \log n = \Theta(n^4)$                   | T   |
| 65  | $\log^3 n + \log^2 n = \Theta(\log^3 n)$       | T   |
| 66  | $\log^2 n + \log^3 n = \Theta(\log^2 n)$       | F   |
| 67  | $n^2 + \sqrt{n} = \Theta(n^2)$                 | T   |
| 68  | $n^2 + \sqrt{n} = \Theta(\sqrt{n})$            | F   |
| 69  | $n\log^2 n + n\log n = \Theta(n\log^2 n)$      | T   |
| 70  | $n\log^2 n + n\log n = \Theta(n\log n)$        | F   |
| 71  | $n^2 + n \sim n^2$                             | T   |
| 72  | $3n^2 \sim n^2$                                | F   |
| 73  | $n + \log n \sim n$                            | T   |
| 74  | $n^3 + 5n^2 \sim n^3$                          | T   |
| 75  | $n^3 + n \sim n^2$                             | F   |
| 76  | $\log n + 1 \sim \log n$                       | T   |
| 77  | $\sqrt{n} + n \sim n$                          | T   |
| 78  | $n(1 + 1/n) \sim n$                            | T   |
| 79  | $n(1 + 1/n)^2 \sim n$                          | T   |
| 80  | $n^2 + o(n^2) \sim n^2$                        | T   |
| 81  | $n^2(1 + o(1)) \sim n^2$                       | T   |
| 82  | $n^2 + \Theta(n) \sim n^2$                     | F   |
| 83  | $n^2 + O(n) \sim n^2$                          | F   |
| 84  | $n^2 + 2n^2 \sim n^2$                          | F   |
| 85  | $n + n^{0.5} \sim n$                           | T   |
| 86  | $n + 3n \sim n$                                | F   |
| 87  | $n\log n + n \sim n\log n$                     | T   |
| 88  | $n\log n + \log n \sim n\log n$                | T   |
| 89  | $n\log n + n^2 \sim n^2$                       | T   |
| 90  | $n^2 + n\log n \sim n^2$                       | T   |
| 91  | $n^2 + \log n \sim n^2$                        | T   |
| 92  | $\log(n^2) \sim 2\log n$                       | T   |
| 93  | $\log_{10} n \sim \log n$                      | T   |
| 94  | $\sqrt{n^2 + n} \sim n$                        | T   |
| 95  | $\sqrt{n^2 + n} \sim n + 0.5$                  | F   |
| 96  | $\overline{X}_n - \mu = O_p(n^{-1/2})$         | T   |
| 97  | $\overline{X}_n - \mu = o_p(1)$                | T   |
| 98  | $\overline{X}_n - \mu = o_p(n^{-1/2})$         | F   |
| 99  | $\overline{X}_n = O_p(1)$                      | T   |
| 100 | $S_n = \sum_{i=1}^n X_i = O_p(n)$              | T   |
| 101 | $S_n = o_p(n)$                                 | F   |
| 102 | $S_n/n = O_p(1)$                               | T   |
| 103 | $S_n/n = o_p(1)$                               | F   |
| 104 | $\sqrt{n}(\overline{X}_n - \mu) = O_p(1)$      | T   |
| 105 | $\sqrt{n}(\overline{X}_n - \mu) = o_p(1)$      | F   |
| 106 | $\overline{X}_n = o_p(n^{-1/2})$               | F   |
| 107 | $\overline{X}_n - \mu = O_p(n^{-1})$           | F   |
| 108 | $\overline{X}_n - \mu = O_p(1)$                | T   |
| 109 | $\overline{X}_n - \mu = o_p(n^{-1})$           | F   |
| 110 | $n(\overline{X}_n - \mu) = O_p(1)$             | F   |
| 111 | $\sqrt{\log n}(\overline{X}_n - \mu) = O_p(1)$ | T   |
| 112 | $S_n / n^{1.1} = o_p(1)$                       | T   |
| 113 | $S_n / n^{0.9} = o_p(1)$                       | F   |
| 114 | $\overline{X}_n^2 = O_p(1)$                    | T   |
| 115 | $\overline{X}_n^2 = o_p(1)$                    | F   |
| 116 | $S_n^2 = O_p(n^2)$                             | T   |
| 117 | $S_n^2 = o_p(n^2)$                             | F   |
| 118 | $\overline{X}_n - \mu = O_p(n^{-0.4})$         | T   |
| 119 | $\overline{X}_n - \mu = o_p(n^{-0.4})$         | F   |
| 120 | $S_n / \sqrt{n} = O_p(1)$                      | T   |
| 121 | $n^3 = o(n^4)$                                 | T   |
| 122 | $n^3 = O(n^4)$                                 | T   |
| 123 | $n^3 = \Theta(n^4)$                            | F   |
| 124 | $n^3 + 5n = \Theta(n^3)$                       | T   |
| 125 | $n^3 + 5n \sim n^3$                            | T   |
| 126 | $\log n = o(n^\epsilon)$                       | T   |
| 127 | $n^\epsilon = o(\log n)$                       | F   |
| 128 | $n\log n = o(n^{1+\epsilon})$                  | T   |
| 129 | $n^2\log n = O(n^3)$                           | T   |
| 130 | $n^2\log n = o(n^3)$                           | T   |
| 131 | $n^2\log n = \Theta(n^2)$                      | F   |
| 132 | $n^2\log n = \Theta(n^2\log_{10} n)$           | T   |
| 133 | $n^{1.5} = o(n^2)$                             | T   |
| 134 | $n^{1.5} = O(n^2)$                             | T   |
| 135 | $n^{1.5} = \Theta(n^2)$                        | F   |
| 136 | $n(1 + o(1)) = \Theta(n)$                      | T   |
| 137 | $n(1 + o(1)) \sim n$                           | T   |
| 138 | $n(1 + O(1)) = \Theta(n)$                      | T   |
| 139 | $n(1 + O(1)) \sim n$                           | F   |
| 140 | $n + o(n) = \Theta(n)$                         | T   |
| 141 | $n + o(n) \sim n$                              | T   |
| 142 | $n + O(n) = \Theta(n)$                         | T   |
| 143 | $n + O(n) \sim n$                              | F   |
| 144 | $n + \Theta(n) = \Theta(n)$                    | T   |
| 145 | $n + \Theta(n) \sim n$                         | F   |
| 146 | $n^4 = O(n^5)$                                 | T   |
| 147 | $n^5 = O(n^4)$                                 | F   |
| 148 | $n^4 = o(n^5)$                                 | T   |
| 149 | $n^5 = o(n^4)$                                 | F   |
| 150 | $n^k = o(n^{k+1})$                             | T   |
